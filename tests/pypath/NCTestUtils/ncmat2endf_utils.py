
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################


from NCrystalDev.ncmat2endf import ncmat2endf
from NCrystalDev._ncmat2endf_impl import _endf_clean
from NCrystalDev.exceptions import NCBadInput
from NCrystalDev._numpy import _np
import NCrystalDev.cli as nc_cli

from .common import ( print_text_file_with_snipping,
                      require_flteq )

def test_cfg( cfg, check_teff=None,
              ref_parsed=None, ref_bragg_edges=None,
              compare_xsec=False,
             **kwargs ):
    import pprint
    print()
    print()
    print()
    print('-'*80)
    if 'force_save' not in kwargs:
        kwargs['force_save'] = True
    kwargs['ncmat_cfg']=cfg
    pprint.pprint(kwargs)
    temperatures = None
    import NCrystal.core as nc_core
    temperatures = (nc_core.createInfo(cfg).dyninfos[0].temperature,)
    if 'temperatures' in kwargs:
        temperatures = temperatures+tuple(kwargs['temperatures'])
    res = ncmat2endf(**kwargs)
    if compare_xsec:
        E = _np.geomspace(1e-5, 5.0, 1000)
        xs_test = _np.zeros(_np.shape(E))
    for endf_fn, frac in res:
        print(f"Created file {endf_fn} with fraction {frac}")
        with open(endf_fn) as f:
            text = "".join(f.readlines())
            print_text_file_with_snipping(text,
                                          nstart=140,
                                          nend=70,
                                          prefix='endf>')
        parser, endf_dic = None, None
        EndfParser, list_parsed_sections = None, None
        if check_teff or ref_parsed or ref_bragg_edges:
            from endf_parserpy import EndfParser, list_parsed_sections
            parser = EndfParser(cache_dir=False)
            endf_dic = parser.parsefile(endf_fn)

        if check_teff:
            teff = endf_dic[7][4]['teff0_table']['Teff0']
            za = endf_dic[1][451]['ZA']
            from NCrystal import atomDB
            label = atomDB(Z=int(za/1000)).elementName()
            teff_vals = [_endf_clean([compute_teff(cfg+f';temp={T}', label)])[0]
                         for T in temperatures]
            require_flteq(teff, teff_vals)

        if ref_parsed:
            if endf_fn not in ref_parsed.keys():
                raise RuntimeError( 'No reference parsed ENDF sections for '
                                   f'{endf_fn}')
            parsed = " ".join([" ".join(str(x) for x in _)
                               for _ in list_parsed_sections(endf_dic)])
            if parsed != ref_parsed[endf_fn]:
                raise RuntimeError( 'ENDF sections {parsed} expected but '
                                   f'sections {ref_parsed[endf_fn]} found')
        if ref_bragg_edges:
            if endf_fn not in ref_bragg_edges.keys():
                raise RuntimeError( 'No reference Bragg edges for '
                                   f'{endf_fn}')
            ref_Eint, ref_S0 = ref_bragg_edges[endf_fn]
            Eint = tuple(endf_dic[7][2]['S_T0_table']['Eint'][:len(ref_Eint)])
            S0 = tuple(endf_dic[7][2]['S_T0_table']['S'][:len(ref_Eint)])
            require_flteq(Eint, ref_Eint)
            require_flteq(S0, ref_S0[0])
            if temperatures:
                Slist = tuple(tuple(S for k, S in Sdict.items())
                              for tindex, Sdict in
                              endf_dic[7][2]['S'].items())
                for v1,v2 in zip(Slist[:len(ref_Eint)], ref_S0[1]):
                    v1 = tuple(v1)
                    require_flteq(v1, v2)
        if compare_xsec:
            xs_test += frac*get_scatxs_from_endf(endf_fn, E)
    if compare_xsec:
        from NCrystal import load as nc_load
        m = nc_load(cfg+';comp=inelas')
        xs = m.scatter.xsect(E)
        require_flteq(xs, xs_test, tol=0.02)

def compute_teff(cfg, label):
    import NCrystal.core as nc_core
    import NCrystal.vdos as nc_vdos
    info_obj = nc_core.createInfo(cfg)
    di = info_obj.findDynInfo(label)
    emin = di.vdosData()[0][0]
    emax = di.vdosData()[0][1]
    rho = di.vdosData()[1]
    res = nc_vdos.analyseVDOS(emin, emax, rho, di.temperature,
                              di.atomData.averageMassAMU())
    return res['teff']

def test_cfg_fail( e, *args, **kwargs ):
    try:
        test_cfg(*args,**kwargs)
    except NCBadInput as e:
        print("FAILED (as expected): %s"%e)
        return
    raise SystemExit('Did not fail as expected')

def test_cli( args ):
    import shlex
    if isinstance(args,str):
        args = shlex.split(args)
    hr=f"============= CLI >>{shlex.join(args)}<< ===================="
    print(hr)
    nc_cli.run('ncmat2endf',*args)
    print('='*len(hr))

def get_scatxs_from_endf(endf_fn, E=None):
    """
    Computes scattering XS from first temperature in ENDF-6 TSL file
    """
    from endf_parserpy import EndfParser
    from NCrystal import NCMATComposer, atomDB

    parser = EndfParser(cache_dir=False)
    endf_dic = parser.parsefile(endf_fn)
    T0 = 293.6
    emax = endf_dic[1][451]['EMAX']
    za = endf_dic[1][451]['ZA']
    label = atomDB(Z=int(za/1000)).elementName()
    awr = endf_dic[7][4]['AWR']
    lat = endf_dic[7][4]['LAT']
    S_table =  endf_dic[7][4]['S_table']
    S = _np.array([v['S'] for k,v in S_table.items()])
    T = endf_dic[7][4]['teff0_table']['Tint'][0]
    beta = _np.array([ v for k, v in endf_dic[7][4]['beta'].items()])
    alpha = _np.array(S_table[1]['alpha'])*awr
    print(alpha)
    if lat == 1:
        beta = beta*T0/T
        alpha = alpha*T0/T
    c_test = NCMATComposer()
    c_test.set_dyninfo_scatknl(label,
                               alphagrid=alpha,
                               betagrid=beta,
                               temperature=T,
                               sab_scaled=S,
                               fraction=1.0)
    c_test.set_density(1.0,'g/cm3')
    m = c_test.load()
    E = _np.geomspace(1e-5, emax, 1000) if E is None else E
    return m.scatter.xsect(E)

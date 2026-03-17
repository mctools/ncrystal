
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
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

# Utilities for detecting changes in given MiniMC output, and making it easy to
# inspect differences and update reference histograms in case of changes.
import NCrystalDev.minimc as ncminimc
from NCrystalDev.exceptions import NCBadInput
from NCrystalDev.hist import Hist1D
from NCrystalDev.core import load as ncload

#Specialisation for std case with which we might want to test many (or even all)
#materials.:
def minimc_unittest_stdsphere( *,
                               cfgstr,
                               neutron_energy,#string like "1eV" or "1Aa"
                               illuminate_uniformly = False,
                               sphere_diam_meter = None,
                               n = None,
                               **kwargs #passed to minimc_unittest
                               ):
    mat = ncload( cfgstr )
    ekin = _parse_energy( neutron_energy )
    if sphere_diam_meter is None:
        ( sphere_diam_str,
          sphere_diam ) = _approx_mfp_as_length_string( mat, ekin = ekin )
    else:
        sphere_diam_str = '%.14g'%sphere_diam_meter
        sphere_diam = float(sphere_diam_str)
        assert sphere_diam == sphere_diam_meter

    sphere_radius = sphere_diam/2
    if n is None:
        n = 1e4 if mat.scatter.isOriented() else 1e5
    srcz = (-sphere_radius)*(1-1e-13)

    if illuminate_uniformly:
        srccfg = f'circular;r={sphere_radius}'
    else:
        srccfg = 'constant'
    srccfg += f';ekin={ekin};z={srcz};n={n}'
    minimc_unittest( cfgstr = cfgstr,
                     srccfg = srccfg,
                     geomcfg = f'sphere;r={sphere_radius}',
                     **kwargs )

def minimc_unittest_scenariostr( cfgstr, scenariostr, override_n = None,
                                 **kwargs ):
    d = ncminimc.decode_scenario( cfgstr, scenariostr )
    for k,v in sorted(d.items()):
        assert isinstance(k,str)
        assert isinstance(v,str)
        print(f'  decoded {k}="{v}"')
    if override_n is not None:
        ns = str(override_n)
        assert float(ns)==int(float(ns))
        ns = f';n={ns}'
        print(f'  override_n set => Appending "{ns}" to src_cfg')
        d['srccfg'] += ns
    return minimc_unittest( cfgstr = cfgstr,
                            srccfg = d['srccfg'],
                            geomcfg = d['geomcfg'],
                            **kwargs )

def _parse_sysargv():
    import sys
    args = sys.argv[1:]
    if '--help' in args or '-h' in args:
        print("Possible options: --plot or --update")
    return dict( do_plot = '--plot' in args,
                 do_updateref = '--update' in args )

def main_minimc_unittest_stdsphere( *a, **kw ):
    minimc_unittest_stdsphere( *a, **kw, **_parse_sysargv() )

def main_minimc_unittest_scenariostr( *a, **kw ):
    minimc_unittest_scenariostr( *a, **kw, **_parse_sysargv() )

def main_minimc_unittest( *a, **kw ):
    minimc_unittest( *a, **kw, **_parse_sysargv() )

def _detect_caller_filebasename( drop_ext = True ):
    #Find first file in callstack outside the present one.
    import inspect
    import os
    res = None
    for e in inspect.stack():
        if not os.path.samefile(__file__,e.filename):
            res = e.filename
            break
    assert res
    res = os.path.basename( res )
    if drop_ext:
        res = os.path.splitext( res )[0]
    return res

def minimc_unittest( *,
                     key = '<auto>', #<auto> expands to name taken from
                                     #filename of caller.
                     cfgstr,
                     srccfg,
                     geomcfg,
                     do_plot = False,
                     do_updateref = False,
                     tally='theta',
                     tallybins=None,
                     extra_enginecfg='',
                     quiet = False
                    ):
    """Run MiniMC with chosen configuration and compare result with reference
    histogram. Will keep reference histogram in <testsdata>mmcrefs/key.json (if
    key contains a slash, a subdir will be used). If do_plot is True, show
    plot vs. ref., and if do_updateref is True, update the reference.

    Runs with nthreads=2.
    """
    from .dirs import get_named_test_data_dir

    if '/' not in key:
        key = f'mmcref/{key}'
    if '<auto>' in key:
        autokey = _detect_caller_filebasename()
        if autokey.startswith('sb_nctest'):
            autokey = autokey.split('_',2)[2]
            assert autokey.startswith('test')
            autokey=autokey[len('test'):]
        key = key.replace('<auto>',autokey)
    assert key.count('/') <= 1, "no more than one slash allowed in key"
    if '/' in key:
        subdir, key = key.split('/')
        assert subdir and key
    else:
        subdir, key = None, key

    reffile = get_named_test_data_dir(
        subdir,
        for_updates = do_updateref
    ).joinpath( key+'.json' )

    if tallybins is None:
        tallybins=';tallybins=theta:90:0:180' if tally=='theta' else ''
    else:
        tallybins=f';tallybins={tallybins}'

    enginecfg = ( 'nthreads=2'
                  f';tally={tally}{tallybins}'
                  f';{extra_enginecfg}' )

    res = ncminimc.run( cfgstr=cfgstr,
                        geomcfg=geomcfg,
                        srccfg=srccfg,
                        enginecfg=enginecfg )

    if len(res.tally_names)>1:
        tally = tally.split(',')[0].strip()
        assert tally in res.tally_names
        print("WARNING: Multiple tallies enabled. Will only actually"
              f" compare, plot or persistify the '{tally}' tally (since"
              " it was listed first)")

    h = res.tally(tally).hist_total
    if do_plot:
        h.dump(contents=False)
        res.tally(tally).plot(logy=True)

    if do_updateref:
        reffile.write_text(h.to_json())
        print(f"Updated {reffile}")
        return

    if not reffile.is_file():
        raise SystemExit('Reffile not found (run with --update to create):'
                         f'\n\n  {reffile}\n\n')
    h_ref = Hist1D(reffile.read_text())
    pval = h.check_compat( h_ref, return_pval = True )
    if do_plot or not quiet:
        print( "Pvalue for comp. with ref"
               f" (higher is more compatible): {pval:g}" )
    if do_plot:
        h_ref.plot(do_show=False,error_bands=1.0,
                   alpha=0.3,color='blue',label='ref')
        plt=h.plot(do_show=False,color='none',logy=True,label='new')
        plt.legend()
        plt.grid()
        plt.show()
    if pval < 0.001:
        raise SystemExit("ERROR: Possible compatibility issues detected!")

def _approx_mfp_as_length_string( mat, **xsect_kwargs ):
    macroxs_scatter = _macroxs_if_isotropic( mat, **xsect_kwargs )
    if not macroxs_scatter:
        return '1cm', 0.01#fallback
    else:
        mfp_scatter = 1.0 / macroxs_scatter
        return ( _encode_length_to_str(mfp_scatter,round2digits=True),
                 mfp_scatter )

def _macroxs_if_isotropic( mat, **xsect_kwargs ):
    #macroxs_scatter in units of [1/m]
    return ( None if mat.scatter.isOriented()
             else ( mat.info.factor_macroscopic_xs
                    * mat.scatter.xsect( **xsect_kwargs )
                    / _parse_length('1cm') ) )

#fixme: do we need these here??:
_length_units = {'km':1000.0,
                 'm':1.0,
                 'meter':1.0,
                 'cm':0.01,
                 'mm':0.001,
                 'mfp': None,#special
                 'nm':1e-9,
                 'aa':1e-10,
                 'Aa':1e-10,
                 'AA':1e-10,
                 'angstrom':1e-10}

_energy_units = {'eV':1.0,
                 'keV':1e3,
                 'MeV':1e6,
                 'GeV':1e9,
                 'meV':0.001,
                 'neV':1e-9,
                 'aa':None,
                 'Aa':None,
                 'AA':None,
                 'angstrom':None}

def _tofloat(s):
    try:
        return float(s)
    except ValueError:
        return None

def _parse_unit(valstr,unitmap):
    v = _tofloat(valstr)
    if v is not None:
        return v, None, None
    valstr=valstr.strip()
    for unit,unitvalue in sorted(unitmap.items(),key=lambda x : (-len(x),x)):
        if valstr.endswith(unit):
            v = _tofloat(valstr[:-len(unit)])
            if v is not None:
                return v, unit, unitvalue
    return None,None,None

def _parse_energy( valstr ):
    v,u,uv = _parse_unit( valstr, _energy_units )
    if v is not None and u is None and uv is None:
        raise NCBadInput('Invalid energy specification (missing unit'
                         f' like Aa or eV): "{valstr}"')
    if v is None:
        raise NCBadInput(f'Invalid energy specification: "{valstr}"')
    if u is not None and u.lower() in ('aa','angstrom'):
        from NCrystalDev.constants import wl2ekin
        return wl2ekin(v)
    v *= uv
    return v

def _parse_length( valstr, mfp = None ):
    v,u,uv = _parse_unit( valstr, _length_units )
    if v is not None and u is None and uv is None:
        _ex0="mfp" if mfp is not None else "mm"
        raise NCBadInput('Invalid length specification (missing unit like '
                         f'{_ex0} or cm): "{valstr}"')
    if v is None:
        raise NCBadInput(f'Invalid length specification: "{valstr}"')
    if u=='mfp':
        if mfp is None:
            raise ValueError('Invalid length specification ("mfp" '
                             f'not supported for this parameter): "{valstr}"')
        v *= mfp
    else:
        v *= uv
    return v

def _encode_length_to_str( length_meters, round2digits = False ):
    assert length_meters>=0.0
    if not length_meters:
        return '0mm'
    unit_mm = _parse_length('1mm')
    unit_cm = _parse_length('1cm')
    if round2digits:
        def roundval(x):
            return float('%.2g'%x)
    else:
        def roundval(x):
            return x
    if length_meters <= unit_cm:
        return f'{roundval(length_meters/unit_mm):.14g}mm'
    unit_m = _parse_length('1m')
    assert unit_m == 1.0
    if length_meters <= unit_m:
        return f'{roundval(length_meters/unit_cm):.14g}cm'
    return f'{roundval(length_meters/unit_m):.14g}m'

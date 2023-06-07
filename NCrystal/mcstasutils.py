"""

Module containing various McStas-related utilities, such as functions helping
with automatic setup of McStas-Union materials from NCrystal cfg-strings, or
production of .laz/.lau files for McStas components such as PowderN.comp or
Single_crystal.comp (with reduced physics capabilities of course).

"""

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2023 NCrystal developers                                   ##
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

from . import api as _NC
from ._common import warn as _nc_warn

def cfgstr_detect_components( cfgstr ):
    """
    Internal helper function which can detect which physics components are
    present once a given cfg-string is loaded, and provide derived cfg-strings
    suitable for picking out those components. The return value is a list of
    component names (e.g. "inelas", "cohelas", "sans", etc.) and the associated
    cfg-strings: [(cfgstr1,compname1), (cfgstr2,compname2),...].
    """
    #Normalise original cfg (so syntax errors will refer to the cfgstr as
    #originally specified), and add flags which cause initialisation speedup
    #without affecting which components are present:
    probecfgstr = _NC.normaliseCfg(cfgstr) + ';vdoslux=0'
    res=[]
    from . import misc as nc_misc
    for ct in nc_misc.standard_comp_types:
        extracfg = f';comp={ct}'
        if not _NC.createScatter( probecfgstr + extracfg ).isNull():
            res.append( ( _NC.normaliseCfg(cfgstr+extracfg), ct.replace('_','') ) )
    return res

def cfgstr_2_unioncfg( *, cfgstr, split_by_physics = False ):
    """Analyse cfg-string and return data needed in order to set up Union
    components accordingly. Specifically the function returns
    (inv_pen_depth_2200,proclist) where inv_pen_depth is the inverse penetration
    depth in inverse meters for a neutron at 2200m/s, and proclist is a list
    [(cfgstr1,compname1), (cfgstr2,compname2),...] with the relevant
    components. If split_by_physics is True, this process list will be as
    provided by cfgstr_detect_components(cfgstr), and otherwise it will simply
    be [(normaliseCfg(cfgstr),'total')].
    """
    info=_NC.createInfo(cfgstr)
    inv_pen_depth_2200 = 100.0 * info.getNumberDensity() * info.getXSectAbsorption()# in inverse meters
    if split_by_physics:
        return inv_pen_depth_2200, cfgstr_detect_components( cfgstr )
    else:
        return inv_pen_depth_2200, [(_NC.normaliseCfg(cfgstr),'total')]

def cfgstr_2_union_instrument_code( *, cfgstr, name, split_by_physics = False ):
    """
    Analyse cfg-string (via a call to cfgstr_2_unioncfg(cfgstr,split_by_physics)),
    and use the results to produce (and return) McStas-instrument code suitable
    for defining a McStas-Union material with the given name as specified via
    the "name" parameter. Depending on the value of split_by_physics, this Union
    material will either be set up with either just a single NCrystal_process,
    or a whole list of NCrystal_processes. Using split_by_physics=False is
    likely slightly more efficient, but split_by_physics=True will allow insight
    into the contributions of the various physics components ("inelas",
    "coh_elas", "sans", "incoh_elas", ...) at the McStas-Union level.
    """
    for ch in ';.:':
        if ch in name:
            raise ValueError(f'union process name "{name}" contains a "{ch}"'
                             +' character which is likely not what was intended.')
    assert not '/*' in cfgstr
    assert not '*/' in cfgstr
    assert not '/*' in name
    assert not '*/' in name

    out_absorption, physlist = cfgstr_2_unioncfg( cfgstr=cfgstr, split_by_physics=split_by_physics )
    res = f"""
/*
   The following code was auto generated by NCrystal v{_NC.get_version()} via Python:

     NCrystal.mcstasutils.cfgstr_2_union_instrument_code(
         cfgstr = {repr(str(cfgstr))},
         name = {repr(str(name))}"""
    if split_by_physics:
        res += """,
         split_by_physics = True )
"""
    else:
        res += ' )\n'
    res+="""
   Please rerun in case of major changes to input data or NCrystal.
*/
"""

    procnames=[]
    mcstas_string_max = 256
    for proccfgstr, proctypename in physlist:
        if len(physlist)==1 and proctypename=='total':
            procname=f'{name}_ncrystal_proc'
        else:
            procname=f'{name}_ncrystal_{proctypename}_proc'
        procnames.append( procname )
        if len(proccfgstr) > mcstas_string_max:
            _nc_warn(f'cfg_string might be too long for a McStas string: "{proccfgstr}"')

        res += f"""
COMPONENT {procname} = NCrystal_process(
    cfg = "{proccfgstr}" )
AT (0,0,0) ABSOLUTE
"""
    procnamestr = ','.join(procnames)
    if len(procnamestr) > mcstas_string_max:
        _nc_warn(f'process_string might be too long for a McStas string: "{procnamestr}" '
                 +'(it might work if you use a shorter material name)')
    res += f"""
COMPONENT {name} = Union_make_material(
    process_string = "{procnamestr}",
    my_absorption = {out_absorption:.15g} )
AT (0,0,0) ABSOLUTE

/* End of auto generated code from NCrystal v{_NC.get_version()}. */
"""



    return res

def cfgstr_2_hkl(*, cfgstr, tgtformat, verbose=True, fp_format = '%.14g', prefer_systemexit_exceptions = False ):
    """Function which can be used to create input files with reflections for
    McStas crystalline sample components like PowderN and Single_crystal, based
    on NCrystal cfg-strings (usually referring to NCMAT files with crystalline
    single-phase materials). The tgtformat must be either 'laz' or 'lau'. If
    verbose is True, the files might contain strictly unneccessary content
    (e.g. white-space for adjusting columns, a dspacing column in .lau files,
    ...). Finally the fp_format parameter can be used to change the precision of
    floating point numbers in the file.

    This function returns an iterable, yielding one line of the output file at a
    time.

    """
    import numbers
    import functools
    import math

    assert tgtformat in ('laz','lau')
    doPowder = (tgtformat=='laz')
    doDsp = (doPowder or verbose)
    doMult = True#Always needed, even in lau

    cfgstr = _NC.normaliseCfg(cfgstr)
    info = _NC.createInfo(cfgstr)
    def errmsg(msg):
        raise (SystemExit if prefer_systemexit_exceptions else RuntimeError)('Error: %s'%msg)

    if info.isMultiPhase():
        errmsg('this script does not handle multiphase materials')
    if not info.hasHKLInfo():
        errmsg('this script does not handle non-crystalline materials')
    if not info.hasStructureInfo():
        errmsg('this script does not handle crystalline materials without unit cell structure')
    si = info.getStructureInfo()
    if not si.get('spacegroup',None):
        errmsg('this script does not handle crystalline materials without spacegroup number')
    if not info.hasAtomInfo():
        errmsg('this script does not handle crystalline materials without info of atoms placement in the unit cell')
    if info.hklInfoType() != _NC.HKLInfoType.SymEqvGroup:
        errmsg('this script does not handle crystalline materials without symmetry equivalent HKL groupings')

    orig_header = []
    decoded_cfg = _NC.decodeCfg(cfgstr)
    data_name = decoded_cfg.get('data_name',None)
    cfgstr_nodataname = cfgstr.replace(data_name,'')
    if cfgstr_nodataname.startswith(';'):
        cfgstr_nodataname = cfgstr_nodataname[1:]
    cfgstr_nodataname = cfgstr_nodataname.strip()

    if ( decoded_cfg.get('density',{}).get('type',None) != 'scalefactor'
         or decoded_cfg.get('density',{}).get('value',None) != 1.0 ):
        errmsg('this script does not handle configurations with density overrides'
               +' (since resulting files might not load correctly by all clients)')

    if data_name:
        for i,l in enumerate(_NC.createTextData(data_name)):
            if l.startswith('#') or not l.strip():
                orig_header.append(l.rstrip())
            else:
                if not (i==0 and l.startswith('NCMAT')):
                    break
    assert info.hasStructureInfo()
    assert info.hasHKLInfo()
    assert info.hasTemperature()
    assert info.hasAtomInfo()

    si = info.getStructureInfo()
    fmtfp=lambda x : fp_format%x if isinstance(x, numbers.Real) else str(x)

    yield f'# File created by NCrystal v{_NC.get_version()}'
    yield '#'
    yield f'# ncrystal_cfgstr "{cfgstr}"'
    yield '#'
    if tgtformat=='laz':
        yield '# Format: "laz" (suitable for McStas PowderN component)'
    else:
        assert tgtformat=='lau'
        yield '# Format: "lau" (suitable for McStas Single_crystal component)'
    yield '#'

    if orig_header:
        had_embedded_cfg = any( 'NCRYSTALMATCFG[' in e for e in orig_header)
        guard=' orighdr :'#something that hopefully will cause McStas code to ignore the line
        yield f'#{guard} Original file had the following comments at the top:'
        if had_embedded_cfg:
            yield f'#{guard}'
            yield f'#{guard} (note the ncrystalmatcfg statement is disabled in this'
            yield f'#{guard}  reproduction by adding "<disable>" around it)'
        yield f'#{guard}'
        for l in orig_header:
            if 'NCRYSTALMATCFG[' in l:
                l=l.replace('NCRYSTALMATCFG[','NCRYSTAL<disable>MATCFG[')
            yield f'#{guard} {l}'
        yield '#'
    if cfgstr_nodataname:
        yield '#'
        yield f'# ncrystal_embedded_cfg : NCRYSTALMATCFG[{cfgstr_nodataname}]'
        yield '#'
    yield f'# temperature {fmtfp(info.getTemperature())} [kelvin]'
    yield '#'
    yield f'# spacegroup {si["spacegroup"]}'
    yield f'# lattice_a {fmtfp(si["a"])} [Aa]'
    yield f'# lattice_b {fmtfp(si["b"])} [Aa]'
    yield f'# lattice_c {fmtfp(si["c"])} [Aa]'
    yield f'# lattice_aa {fmtfp(si["alpha"])} [degrees]'
    yield f'# lattice_bb {fmtfp(si["beta"])} [degrees]'
    yield f'# lattice_cc {fmtfp(si["gamma"])} [degrees]'
    yield f'# Vc {fmtfp(si["volume"])} [Aa^3]'
    yield '#'
    n_atoms = sum( ai.count for ai in info.atominfos )
    assert si["n_atoms"] == n_atoms
    yield f'# multiplicity {n_atoms} [atoms/unit cell]'
    yield f'# density {fmtfp(info.density)} [g/cm^3]'
    yield f'# number_density {fmtfp(info.numberdensity)} [atoms/Aa^3]'
    daltons_per_unitcell = sum(ai.atomData.averageMassAMU() * ai.count for ai in info.atominfos)
    yield f'# weight {fmtfp(daltons_per_unitcell)} [g/mol of entire unit cell]'
    yield f'# average_mass {fmtfp(daltons_per_unitcell/n_atoms)} [average atomic g/mol]'

    sigma_coh = sum(ai.atomData.coherentXS() * ai.count for ai in info.atominfos)
    sigma_inc = sum(ai.atomData.incoherentXS() * ai.count for ai in info.atominfos)
    sigma_abs = sum(ai.atomData.captureXS() * ai.count for ai in info.atominfos)
    yield '#'
    yield f'# sigma_coh {fmtfp(sigma_coh)} [barn/unitcell]'
    yield f'# sigma_inc {fmtfp(sigma_inc)} [barn/unitcell]'
    yield f'# sigma_abs {fmtfp(sigma_abs)} [barn/unitcell]'
    yield '#'

    has_debye_temp = True
    debye_temp_sum = 0.0
    debye_temp_sumw = 0.0

    if all(ai.atomData.isElement() for ai in info.atominfos):
        d = {}
        for ai in info.atominfos:
            if has_debye_temp:
                if ai.debyeTemperature is not None:
                    _dt = ai.debyeTemperature
                elif ai.msd is not None:
                    _dt = _NC.debyeTempFromIsotropicMSD( msd = ai.msd,
                                                        temperature = info.getTemperature(),
                                                        mass = ai.atomData.averageMassAMU() )
                else:
                    has_debye_temp = False
                    _dt = 0.0
                _dtw = ai.count * ai.atomData.scatteringXS()
                debye_temp_sum += _dt*_dtw
                debye_temp_sumw += _dtw

            d[ai.atomData.elementName()] = d.get(ai.atomData.elementName(),0) + ai.count
            nformula_per_unitcell = functools.reduce(math.gcd, list(c for _,c in d.items()))
        formula = ''.join(( '%s%i'%(k,v/nformula_per_unitcell) if v!=nformula_per_unitcell else k) for k,v in sorted(d.items()))
        yield f'# formula {formula}'
        yield f'# nformula_per_unitcell {nformula_per_unitcell}'
        if debye_temp_sum:
            yield f'# debye_temperature {fmtfp(debye_temp_sum/debye_temp_sumw)} [kelvin, weighted by sigma_b]'

    cols = ['h','k','l']
    if doMult:
        cols.append('j')
    if doDsp:
        cols.append('d')
    cols.append('F2')

    colwidths={'h':3,'k':3,'l':3,'j':2,'d':16,'F2':16}
    header = '#'
    for i,c in enumerate(cols):
        header += ' '
        if verbose:
            header += c.rjust(colwidths[c])
        else:
            header += c

    col_comments={
        'd':"[d-spacing in Aa]",
        'j': "multiplicity",
        'F2':"[norm of scattering factor |F|^2 in barn]"
    }

    yield '#'
    maxcw=max(len(c) for c in col_comments.keys())
    for i,c in enumerate(cols):
        cc = col_comments.get(c,'')
        if cc:
            cc = '  '+cc
        yield f'# column_{c.ljust(maxcw)} {i+1}{cc}'
    yield '# unit_F2 barn'
    yield '#'

    yield header
    def format_line(data):
        s=' ' if verbose else ''
        for c in cols:
            if s:
                s+=' '
            s+=fmtfp(data[c]).rjust(colwidths[c])
        return s

    if doPowder:
        for h,k,l,mult,dsp,F2_barn in info.hklList():
            yield format_line({'h':h,'k':k,'l':l,'d':dsp,'j':mult,'F2':F2_barn})
    else:
        for h,k,l,mult,dsp,F2_barn in info.hklList(all_indices=True):
            for i in range(len(h)):
                #Yield both hkl and -hkl and put mult=1 (since we put all planes explicitly)
                yield format_line({'h':h[i],'k':k[i],'l':l[i],'d':dsp,'j':1,'F2':F2_barn})
                yield format_line({'h':-h[i],'k':-k[i],'l':-l[i],'d':dsp,'j':1,'F2':F2_barn})

def main(argv):
    args = argv[1:]
    if args and isinstance(args[0],bytes):
        args = list(e.decode() for e in args)
    def usage(*,err):
        if err:
            print("ERROR - wrong usage!")
            print()
        import os
        pn=os.path.basename(argv[0])
        if pn.endswith('.py'):
            pn = pn[0:-3]
        print("""Usage:

  --help|-h                         Show these instructions

  --union [--split] NAME CFGSTR     Output McStas-Union code to define material
                                    with NAME based on specified NCrystal CFGSTR.
                                    Providing --split will split NCrystal processes
                                    into physics types at the McStas-Union level,
                                    rather than internally in NCrystal.

  --laz CFGSTR                      Convert NCrystal CFGSTR to McStas .laz format
                                    for PowderN.

  --lau CFGSTR                      Convert NCrystal CFGSTR to McStas .lau format
                                    for Single_crystal.

""")
        raise SystemExit(1 if err else 0)
    err = lambda : usage(err=True)
    if '-h' in args or '--help' in args:
        usage(err=False)
    if not args or len(args)<2:
        err()
    if args[0] == '--union':
        dosplit = False
        if '--split' in args:
            args.remove('--split')
            dosplit=True
        if len(args)<3:
            err()
        print(cfgstr_2_union_instrument_code( cfgstr=';'.join(args[2:]), name=args[1], split_by_physics = dosplit ))
        return
    if args[0] in ('--laz','--lau'):
        for e in cfgstr_2_hkl( cfgstr=';'.join(args[1:]),tgtformat=args[0][2:], prefer_systemexit_exceptions = True):
            print(e)
        return
    err()

if __name__ == '__main__':
    import sys
    main(sys.argv)

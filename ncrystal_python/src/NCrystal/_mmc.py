
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

# For backwards compatibility we keep __all__ here in this internal file. In
# NCrystal 4.0.0 we actually renamed the public API file from _mmc.py to mmc.py.

__all__ = ['runsim_diffraction_pattern',
           'quick_diffraction_pattern',
           'quick_diffraction_pattern_autoparams']

from ._mmc_utils import MMCResults#fixme: MMCResults to main mmc.py
from .exceptions import NCBadInput

def run_minimc( cfgstr, *, geomcfg, srccfg, enginecfg ):
    """Run the NCrystal MiniMC with the provided cfg-strings for material,
    geometry, source and engine. Returns result wrapped in MMCResults object.
    """
    res = raw_minimc_invocation( cfgstr = cfgstr,
                                 geomcfg = geomcfg,
                                 srccfg = srccfg,
                                 enginecfg = enginecfg,
                                 unpack = True )
    return MMCResults(res)

def minimc_decode_scenario( cfgstr, scenario ):
    """fixme"""
    from ._chooks import _get_raw_cfcts
    return _get_raw_cfcts()['minimc_scenario']( cfgstr, scenario )

def raw_minimc_invocation( cfgstr, *,
                           geomcfg, srccfg, enginecfg, unpack=False ):
    """Raw invocation of the ncrystal_minimc C function, accepting 4
    configuration strings and returning a JSON string with the results.

    Set unpack=True to decode the JSON data and replace any Hist1D data inside
    it with Hist1D objects. Set unpack='json' to simply decode the JSON data.
    """
    assert isinstance(cfgstr, str)
    assert isinstance(geomcfg, str)
    assert isinstance(srccfg, str)
    assert isinstance(enginecfg, str)
    from ._chooks import _get_raw_cfcts
    _rawfct = _get_raw_cfcts()
    res = _rawfct['minimc']( cfgstr, geomcfg, srccfg, enginecfg )
    if unpack:
        import json
        res = json.loads(res)
        if unpack != 'json':
            from ._hist import Hist1D
            res = Hist1D.objectify_data(res)
    return res

_cache_availtallies=[None]
def available_tallies():
    """Returns dictionary with information about available MiniMC talies that
    can be used in the enginecfg string like enginecfg="tally=cosmu,de").
"""
    if _cache_availtallies[0] is None:
        from types import MappingProxyType
        import json
        res = MappingProxyType(
            dict( (k,tuple(v)) for k,v in
                  json.loads(raw_minimc_invocation( "[[show_tally_lists]]",
                                                    geomcfg="", srccfg="",
                                                    enginecfg="",
                                                    unpack=False)).items())
        )
        _cache_availtallies[0] = res
    return _cache_availtallies[0]

#def runsim_diffraction_pattern_NEW( cfgstr, *,
#                                    geomcfg,
#                                    srccfg,
#                                    enginecfg = ''):
#    res = raw_minimc_invocation( cfgstr = cfgstr,
#                                 geomcfg = geomcfg,
#                                 srccfg = srccfg,
#                                 enginecfg = enginecfg,
#                                 unpack = True )
#    return res

__cache_rsdwarn=[True]
def runsim_diffraction_pattern( cfgstr, *,
                                geomcfg,
                                srccfg,
                                nthreads = 'auto',
                                tally_detail_lvl = 2,
                                suppress_obsoletion_warning = False ):

    """Run embedded "mini-MonteCarlo" (MiniMC) to produce diffraction pattern,
    including both effects of multiple scattering and absorption physics. This
    for now results in exit angle (angle of emitted neutrons wrt. the Z-axis)
    histograms, with a perfect 4pi detector.

    This is highly experimental, and not yet fully documented.

    Example geomcfg: 'sphere;r=0.1'. This describes a 0.1m=10cm radius sphere
    centered at (0,0,0).

    Example srccfg: 'constant;ekin=0.025;n=1e6;z=-0.1'. This starts 1e6 0.025eV
    neutrons at (0,0,-10cm) with a direction (0,0,1).

    tally_detail_lvl can be reduced to 1 or 0, if only the exit_angle histogram
    is needed. tally_detail_lvl=2 provides more details, including histograms
    for specific components (e.g. separating the contributions from single and
    multiple scattering, and inelastic/elastic scatterings.

    nthreads can be 'auto' or a specific integral value of threads to use.

    """
    #fixme: update docs and emit warning.
    if __cache_rsdwarn[0] and not suppress_obsoletion_warning:
        __cache_rsdwarn[0] = False
        from ._common import warn
        warn('The runsim_diffraction_pattern(..) function is obsolete. Please'
             ' migrate your code to use the run_minimc(..) function instead')

    assert tally_detail_lvl in (0,1,2)
    ec = f'nthreads={nthreads};tally=mu'
    if tally_detail_lvl < 2:
        ec += ',nobreakdown'
    return run_minimc( cfgstr,
                       geomcfg=geomcfg,
                       srccfg=srccfg,
                       enginecfg=ec)
#
#
#    #Fixme: make this fct still work, but emit warning and redirect to new
#    #function!
#    from ._numpy import _ensure_numpy
#    _ensure_numpy()
#
#    nthreads = 9999 if nthreads=='auto' else min(9999,max(1,int(nthreads)))
#
#    setup_info = dict( nthreads = int(nthreads),
#                       tally_detail_lvl = int(tally_detail_lvl),
#                       material_cfgstr = str(cfgstr),
#                       geomcfg = str(geomcfg),
#                       srccfg = str(srccfg) )
#
#
#
#    from ._chooks import _get_raw_cfcts
#    _rawfct = _get_raw_cfcts()
#    cfct = _rawfct['runmmcsim_stdengine']
#
#    ct,errsq,res_json = cfct( nthreads = nthreads,
#                              tally_detail_lvl = int(tally_detail_lvl),
#                              mat_cfgstr = str(cfgstr),
#                              mmc_geomcfg = str(geomcfg),
#                              mmc_srccfg = str(srccfg) )
#
#    from ._mmc_results import MMCDiffractionResults
#    return MMCDiffractionResults._from_C( main_hist_content = ct,
#                                          main_hist_errsq = errsq,
#                                          json_details = res_json,
#                                          cfgstr = str(cfgstr),
#                                          setup_info = setup_info)

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
        from .constants import wl2ekin
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

def _macroxs_if_isotropic( mat, **xsect_kwargs ):
    #macroxs_scatter in units of [1/m]
    return ( None if mat.scatter.isOriented()
             else ( mat.info.factor_macroscopic_xs
                    * mat.scatter.xsect( **xsect_kwargs )
                    / _parse_length('1cm') ) )

def _approx_mfp_as_length_string( mat, **xsect_kwargs ):
    macroxs_scatter = _macroxs_if_isotropic( mat, **xsect_kwargs )
    if not macroxs_scatter:
        return '1cm', 0.01#fallback
    else:
        mfp_scatter = 1.0 / macroxs_scatter
        return ( _encode_length_to_str(mfp_scatter,round2digits=True),
                 mfp_scatter )

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

#fixme: docstrings?
def quick_diffraction_pattern_autoparams( cfgstr ):
    from .core import load as ncload
    mat = ncload( cfgstr )
    neutron_wl = 1.8 # todo: depend on e.g. Bragg threshold?
    material_thickness = _approx_mfp_as_length_string( mat,
                                                       wl = neutron_wl )[0]
    return dict( neutron_energy_str = f'{neutron_wl}Aa',
                 material_thickness_str = material_thickness )

def minimc_quick_tally( cfgstr, tally='mu', nstat='auto' ):
    """fixme"""
    scenario_cfg = '0.8BT pencil on 1mfp sphere'
    cfg = minimc_decode_scenario( cfgstr, scenario_cfg )
    cfg['enginecfg'] += ';tally=mu'
    def simfct( nn, cfgstr ):
        import time
        t0 = time.time()
        cfg2 = dict(cfg.items())
        cfg2['srccfg'] += f';n={nn}'
        res = run_minimc( **cfg2 )
        t1 = time.time()
        return t1-t0, res

    if nstat is None or nstat=='auto':
        for nstat in [1e4,1e5,1e6,1e7]:
            t,res = simfct(nstat,cfgstr)
            #Usually, end within a second in total, but in worst cases, up to
            #10seconds:
            if t>0.1 and nstat >= 1e5:
                break
            if t>1.0:
                break
    else:
        t,res=simfct(nstat,cfgstr)
    return res

__cache_qdpwarn=[True]
def quick_diffraction_pattern( cfgstr, *,
                               neutron_energy,
                               material_thickness,
                               nstat = 'auto',
                               nthreads = 'auto',
                               suppress_obsoletion_warning = False ):

    if __cache_qdpwarn[0] and not suppress_obsoletion_warning:
        __cache_qdpwarn[0] = False
        from ._common import warn
        warn('The quick_diffraction_pattern(..) function is obsolete. Please'
             ' migrate your code to use the run_minimc(..) function instead')

    scenario_cfg = f'{neutron_energy} pencil on {material_thickness} sphere'
    cfg = minimc_decode_scenario( cfgstr, scenario_cfg )
    cfg['enginecfg'] += f';nthreads={nthreads};tally=mu'

    def simfct( n, cfgstr ):
        import time
        t0 = time.time()
        cfg2 = dict(cfg.items())
        cfg2['srccfg'] += f';n={n}'
        res = run_minimc( **cfg2 )
        t1 = time.time()
        return t1-t0, res

    if nstat is None or nstat=='auto':
        for nstat in [1e4,1e5,1e6,1e7]:
            t,res = simfct(nstat,cfgstr)
            #Usually, end within a second in total, but in worst cases, up to
            #10seconds:
            if t>0.1 and nstat >= 1e6:
                break
            if t>1.0:
                break
    else:
        t,res=simfct(nstat,cfgstr)
    return res

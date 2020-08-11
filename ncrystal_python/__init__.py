#!/usr/bin/env python
"""Python module for using the NCrystal library for thermal neutron transport in crystals

Please find more information about NCrystal at the website:

   https://mctools.github.io/ncrystal/

In particular, a small example using the NCrystal python module can be found at:

   https://github.com/mctools/ncrystal/blob/master/examples/ncrystal_example_py

A substantial effort went into developing NCrystal. If you use it for your work,
we would appreciate it if you would use the following reference in your work:

  X.-X. Cai and T. Kittelmann, NCrystal: A library for thermal neutron
  transport, Computer Physics Communications 246 (2020) 106851,
  https://doi.org/10.1016/j.cpc.2019.07.015

For detailed usage conditions and licensing of this open source project, see:

   https://github.com/mctools/ncrystal/blob/master/NOTICE
   https://github.com/mctools/ncrystal/blob/master/LICENSE
   https://github.com/mctools/ncrystal/blob/master/ncrystal_extra/LICENSE

"""

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2020 NCrystal developers                                   ##
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

#enable py3 behaviour in py2.6+ (avoid unicode_literals on purpose!)
from __future__ import division, print_function, absolute_import

__license__ = "Apache 2.0, http://www.apache.org/licenses/LICENSE-2.0"
__version__ = '1.99.1'
__status__ = "Production"
__author__ = "NCrystal developers (Thomas Kittelmann, Xiao Xiao Cai)"
__copyright__ = "Copyright 2015-2020 %s"%__author__
__maintainer__ = __author__
__email__ = "ncrystal-developers@cern.ch"
#Only put the few most important items in __all__, to prevent cluttering on
#wildcard imports. Specifically this is the exceptions, the most important API
#classes, the factory functions, and the constants:
__all__ = [ 'NCException','NCFileNotFound','NCDataLoadError',
            'NCCalcError','NCLogicError','NCBadInput',
            'RCBase','Info','CalcBase','Process','Absorption','Scatter',
            'createInfo','createScatter','createAbsorption',
            'constant_c','constant_dalton2kg','constant_dalton2eVc2',
            'constant_avogadro','constant_boltzmann',
            'const_neutron_mass_amu','constant_planck']

import sys

###################################
#For python2/python3 support:

__metaclass__ = type  #classes are new-style in py2 without inheriting from "object"

try:
    xrange
except NameError:
    pass#in py3, range is py2's xrange and xrange is absent
else:
    range,xrange = xrange,None#emulate py3 range in py2

def _str2cstr(s):
    #converts any string (str,bytes,unicode,path) to bytes
    if hasattr(s,'__fspath__'):
        s=str(s)
    try:
        return s if isinstance(s,bytes) else s.encode('ascii')
    except UnicodeEncodeError:
        #Attempt with file-system encoding, in case of non-ASCII path names:
        return s.encode(sys.getfilesystemencoding())

def _cstr2str(s):
    #converts bytes object to str (unicode in py3, bytes in py2)
    try:
        return s if isinstance(s,str) else s.decode('ascii')
    except UnicodeDecodeError:
        return s.decode(sys.getfilesystemencoding())

###################################
#Same as NCRYSTAL_VERSION macro:
version_num = sum(int(i)*j for i,j in zip(__version__.split('.'),(1000000,1000,1)))

class NCException(RuntimeError):
    """Base class for all exceptions raised by NCrystal code"""
    pass
class NCFileNotFound(NCException):
    pass
class NCDataLoadError(NCException):
    pass
class NCMissingInfo(NCException):
    pass
class NCCalcError(NCException):
    pass
class NCLogicError(NCException):
    pass
class NCBadInput(NCException):
    pass


#some constants (NB: Copied here from NCMath.hh - must keep synchronized!! Also,
#remember to include in __all__ list above):
constant_c  = 299792458e10#  speed of light in Aa/s
constant_dalton2kg =  1.660539040e-27#  amu to kg
constant_dalton2eVc2 =  931494095.17#  amu to eV/c^2
constant_avogadro = 6.022140857e23#  mol^-1
constant_boltzmann = 8.6173303e-5#  eV/K
const_neutron_mass_amu = 1.00866491588#  [amu]
constant_planck = 4.135667662e-15 # [eV*s]

def _find_nclib():

    import glob, os

    #If NCRYSTAL_LIB env var is set, we try that and only that:
    override=os.getenv('NCRYSTAL_LIB')
    if override:
        if not os.path.exists(override):
            raise NCFileNotFound('NCRYSTAL_LIB environment variable is set but does not point to an actual file.')
        return override

    #Try to locate the lib via various library paths:
    pathenvs=['LD_LIBRARY_PATH','DYLD_LIBRARY_PATH','SHLIB_PATH','LIBPATH']
    searchpaths = []
    for libpath in pathenvs:
        for d in [_f for _f in os.environ.get(libpath,'').split(':') if _f]:
            if os.path.isdir(d):
                searchpaths += [d]

    #Try to locate the lib by looking relatively to this python file (makes
    #standard NCrystal installation work out of the box + makes it possible to
    #simply put the shared library next to the python file):
    searchpaths += [ os.path.dirname(__file__),
                     os.path.join(os.path.dirname(__file__),'../lib'),
                     os.path.join(os.path.dirname(__file__),'../lib64'),
                     os.path.join(os.path.dirname(__file__),'../../lib'),
                     os.path.join(os.path.dirname(__file__),'../../lib64'),
                     os.path.join(os.path.dirname(__file__),'../bin'),#on Windows, seems NCrystal.dll ends up here?
                     '/usr/lib' ]

    #Look inside searched directories for libNCrystal:
    libname_patterns = ['libNCrystal.so','libNCrystal.dylib','NCrystal.dll','libNCrystal.*']
    for libnamepat in libname_patterns:
        for d in searchpaths:
            for f in sorted(glob.glob(os.path.join(d,libnamepat))):
                if os.path.exists(f) and not os.path.isdir(f):
                    #Note os.path.exists returns False for broken symlinks, True for
                    #both actual files and unbroken symlinks, so the check here is
                    #sufficient.
                    return f

    raise NCFileNotFound('Could not find NCrystal shared library (specify its location with NCRYSTAL_LIB env var)')

try:
    import numpy as _np
except ImportError:
    _np = None
def _ensure_numpy():
    if not _np:
        raise NCException("Numpy not available - array based functionality is unavailable")

_globalstates = {}

def _load(nclib_filename):
    import ctypes

    _nclib = ctypes.CDLL(nclib_filename)
    _int,_intp,_uint,_uintp,_dbl,_dblp,_cstr,_voidp = (ctypes.c_int, ctypes.POINTER(ctypes.c_int),
                                                       ctypes.c_uint,ctypes.POINTER(ctypes.c_uint), ctypes.c_double,
                                                       ctypes.POINTER(ctypes.c_double), ctypes.c_char_p, ctypes.c_void_p)
    _cstrp = ctypes.POINTER(_cstr)
    _dblpp = ctypes.POINTER(_dblp)
    ndarray_to_dblp = lambda a : a.ctypes.data_as(_dblp)

    _npempty = [None]
    def _create_numpy_array(n):
        _ensure_numpy()
        a=_np.empty(n)
        return a,ndarray_to_dblp(a)

    class ncrystal_info_t(ctypes.Structure):
        _fields_ = [('internal', _voidp)]
    class ncrystal_process_t(ctypes.Structure):
        _fields_ = [('internal', _voidp)]
    class ncrystal_scatter_t(ctypes.Structure):
        _fields_ = [('internal', _voidp)]
    class ncrystal_absorption_t(ctypes.Structure):
        _fields_ = [('internal', _voidp)]

    functions = {}

    #Exceptions:
    _errmap = {'FileNotFound':NCFileNotFound,
               'DataLoadError':NCDataLoadError,
               'MissingInfo':NCMissingInfo,
               'CalcError':NCCalcError,
               'LogicError':NCLogicError,
               'BadInput':NCBadInput}

    def _raise_err():
        assert _ncerror()#checks there was an error
        tm=(_cstr2str(_ncerror_type()),_cstr2str(_ncerror_msg()))
        _ncerror_clear()
        #TODO for NC2: Provide line number / file as well?
        e=_errmap.get(tm[0],NCException)(tm[1])
        e.message = tm[1]#to avoid warnings in py 2.6
        raise e

    #helper class for exporting the functions:

    def _wrap(fct_name,restype,argtypes,take_ref = False, hide=False, error_check=True):
        assert isinstance(argtypes,tuple)
        raw=getattr(_nclib,fct_name)
        raw.argtypes=argtypes
        raw.restype=restype

        if take_ref:
            assert len(argtypes)==1
            fct = lambda arg : raw(ctypes.byref(arg))
        else:
            fct = lambda *args : raw(*args)
        if error_check:
            #NB: we should read about return types in the ctypes tutorial. Apparently one
            #can just set an error checking function as the restype.
            raw_fct = fct
            def fcte(*aaa):
                r = raw_fct(*aaa)
                if _ncerror():
                    _raise_err()
                return r
            fct=fcte
        if not hide:
            functions[fct_name] = fct
        return fct

    lib_version = _cstr2str(_wrap('ncrystal_version_str',_cstr,tuple(),hide=True,error_check=False)())
    if lib_version != __version__:
        print("WARNING: Version mismatch detected between NCrystal python code (v%s)"
              " and loaded binary"" library (v%s). Control which NCrystal library"
              " to load with the NCRYSTAL_LIB env var."%(__version__,lib_version))

    _wrap('ncrystal_sethaltonerror',_int,(_int,),hide=True,error_check=False)(False)
    _wrap('ncrystal_setquietonerror',_int,(_int,),hide=True,error_check=False)(True)
    _ncerror       = _wrap('ncrystal_error',_int,tuple(),hide=True,error_check=False)
    _ncerror_msg   = _wrap('ncrystal_lasterror',_cstr,tuple(),hide=True,error_check=False)
    _ncerror_type  = _wrap('ncrystal_lasterrortype',_cstr,tuple(),hide=True,error_check=False)
    _ncerror_clear = _wrap('ncrystal_clearerror',None,tuple(),hide=True,error_check=False)

    _wrap('ncrystal_refcount',_int,(_voidp,),take_ref=True)
    _wrap('ncrystal_unref',None,(_voidp,),take_ref=True)

    _wrap('ncrystal_cast_scat2proc',ncrystal_process_t,(ncrystal_scatter_t,))
    _wrap('ncrystal_cast_abs2proc',ncrystal_process_t,(ncrystal_absorption_t,))


    _wrap('ncrystal_dump',None,(ncrystal_info_t,))
    _wrap('ncrystal_ekin2wl',_dbl,(_dbl,))
    _wrap('ncrystal_wl2ekin',_dbl,(_dbl,))
    _raw_natelemdata=_wrap('ncrystal_natelemdata',None,( _uint, _cstrp,_dblp,_dblp,_dblp,_dblp),hide=True)
    def ncrystal_natelemdata(z):
        mass_amu,sigma_inc,sigma_coh,sigma_abs=_dbl(),_dbl(),_dbl(),_dbl()
        name=_cstr()
        _raw_natelemdata(z,ctypes.byref(name),mass_amu,sigma_inc,sigma_coh,sigma_abs)
        if mass_amu.value==0.0:
            return {}
        return dict(name=_cstr2str(name.value),
                    z=z,
                    mass_amu=mass_amu.value,
                    sigma_inc=sigma_inc.value,
                    sigma_coh=sigma_coh.value,
                    sigma_abs=sigma_abs.value)
    functions['ncrystal_natelemdata'] = ncrystal_natelemdata
    _wrap('ncrystal_isnonoriented',_int,(ncrystal_process_t,))
    _wrap('ncrystal_name',_cstr,(ncrystal_process_t,))

    for s in ('temperature','xsectabsorption','xsectfree','globaldebyetemp','density','numberdensity'):
        _wrap('ncrystal_info_get%s'%s,_dbl,(ncrystal_info_t,))
    _raw_info_getstruct = _wrap('ncrystal_info_getstructure',_int,(ncrystal_info_t,_uintp,_dblp,_dblp,_dblp,_dblp,_dblp,_dblp,_dblp,_uintp))
    def ncrystal_info_getstructure(nfo):
        sg,natom=_uint(),_uint()
        a,b,c,alpha,beta,gamma,vol = _dbl(),_dbl(),_dbl(),_dbl(),_dbl(),_dbl(),_dbl(),
        if _raw_info_getstruct(nfo,sg,a,b,c,alpha,beta,gamma,vol,natom) == 0:
            return {}
        return dict(spacegroup=int(sg.value),a=a.value,b=b.value,c=c.value,alpha=alpha.value,
                    beta=beta.value,gamma=gamma.value,volume=vol.value,n_atoms=int(natom.value))
    functions['ncrystal_info_getstructure'] = ncrystal_info_getstructure

    _wrap('ncrystal_info_nhkl',_int,(ncrystal_info_t,))
    _wrap('ncrystal_info_hkl_dlower',_dbl,(ncrystal_info_t,))
    _wrap('ncrystal_info_hkl_dupper',_dbl,(ncrystal_info_t,))
    _wrap('ncrystal_info_gethkl',None,(ncrystal_info_t,_int,_intp,_intp,_intp,_intp,_dblp,_dblp))
    _wrap('ncrystal_info_dspacing_from_hkl',_dbl,(ncrystal_info_t,_int,_int,_int))
    functions['ncrystal_info_gethkl_setuppars'] = lambda : (_int(),_int(),_int(),_int(),_dbl(),_dbl())

    _wrap('ncrystal_info_ndyninfo',_uint,(ncrystal_info_t,))
    _raw_di_base = _wrap('ncrystal_dyninfo_base',None,(ncrystal_info_t,_uint,_dblp,_cstrp,_dblp,_uintp),hide=True)
    _raw_di_scatknl = _wrap('ncrystal_dyninfo_extract_scatknl',None,(ncrystal_info_t,_uint,_uint,_dblp,_uintp,_uintp,_uintp,
                                                                     _dblpp,_dblpp,_dblpp,_dblpp),hide=True)
    _raw_di_vdos = _wrap('ncrystal_dyninfo_extract_vdos',None,(ncrystal_info_t,_uint,_dblp,_dblp,_uintp,_dblpp),hide=True)
    _raw_di_vdosdebye = _wrap('ncrystal_dyninfo_extract_vdosdebye',None,(ncrystal_info_t,_uint,_dblp),hide=True)
    _raw_di_vdos_input = _wrap('ncrystal_dyninfo_extract_vdos_input',None,(ncrystal_info_t,_uint,_uintp,_dblpp,_uintp,_dblpp),hide=True)

    def ncrystal_dyninfo_base(key):
        infoobj,dynidx = key
        fr,tt,en,ditype=_dbl(),_dbl(),_cstr(),_uint()
        _raw_di_base(infoobj,dynidx,fr,ctypes.byref(en),tt,ditype)
        return (fr.value,tt.value,_cstr2str(en.value),ditype.value)
    def ncrystal_dyninfo_extract_scatknl(key,vdoslux):
        infoobj,dynidx = key
        sugEmax,ne,na,nb,e,a,b,sab = _dbl(),_uint(),_uint(),_uint(),_dblp(),_dblp(),_dblp(),_dblp()
        _raw_di_scatknl(infoobj,dynidx,vdoslux,sugEmax,ne,na,nb,
                        ctypes.byref(e),ctypes.byref(a),ctypes.byref(b),ctypes.byref(sab))
        return (sugEmax.value,ne.value,na.value,nb.value,e,a,b,sab)
    def ncrystal_dyninfo_extract_vdos(key):
        infoobj,dynidx = key
        egrid_min,egrid_max,ndensity,densityptr = _dbl(),_dbl(),_uint(),_dblp()
        _raw_di_vdos(infoobj,dynidx,egrid_min,egrid_max,ndensity,ctypes.byref(densityptr))
        return (egrid_min.value,egrid_max.value,ndensity.value,densityptr)
    def ncrystal_dyninfo_extract_vdosdebye(key):
        infoobj,dynidx = key
        td=_dbl()
        _raw_di_vdosdebye(infoobj,dynidx,td)
        return td.value
    def ncrystal_dyninfo_extract_vdos_input(key):
        infoobj,dynidx = key
        negrid,egridptr,ndensity,densityptr = _uint(),_dblp(),_uint(),_dblp()
        _raw_di_vdos_input(infoobj,dynidx,negrid,ctypes.byref(egridptr),ndensity,ctypes.byref(densityptr));
        return (negrid.value,egridptr,ndensity.value,densityptr)
    functions['ncrystal_dyninfo_base'] = ncrystal_dyninfo_base
    functions['ncrystal_dyninfo_extract_scatknl'] = ncrystal_dyninfo_extract_scatknl
    functions['ncrystal_dyninfo_extract_vdos'] = ncrystal_dyninfo_extract_vdos
    functions['ncrystal_dyninfo_extract_vdosdebye'] = ncrystal_dyninfo_extract_vdosdebye
    functions['ncrystal_dyninfo_extract_vdos_input'] = ncrystal_dyninfo_extract_vdos_input

    def _prepare_many(ekin,repeat):
        if _np is None and not repeat is None:
            raise NCBadInput('Can not use "repeat" parameter when Numpy is absent on the system')
        if repeat is None and not hasattr(ekin,'__len__'):
            return None#scalar case, array interface not triggered
        repeat = 1 if repeat is None else repeat
        ekin = (ekin if hasattr(ekin,'ctypes') else _np.asfarray(ekin) ) if hasattr(ekin,'__len__') else _np.ones(1)*ekin
        #NB: returning then ekin object itself is important in order to keep a reference to it after the call:
        return ndarray_to_dblp(ekin),len(ekin),repeat,ekin

    _raw_xs_no = _wrap('ncrystal_crosssection_nonoriented',None,(ncrystal_process_t,_dbl,_dblp),hide=True)
    _raw_xs_no_many = _wrap('ncrystal_crosssection_nonoriented_many',None,(ncrystal_process_t,_dblp,ctypes.c_ulong,
                                                                           ctypes.c_ulong,_dblp),hide=True)
    def ncrystal_crosssection_nonoriented(scat,ekin,repeat=None):
        many = _prepare_many(ekin,repeat)
        if many is None:
            res = _dbl()
            _raw_xs_no(scat,ekin,res)
            return res.value
        else:
            ekin_ct,n_ekin,repeat,ekin_nparr = many
            xs, xs_ct = _create_numpy_array(n_ekin*repeat)
            _raw_xs_no_many(scat,ekin_ct,n_ekin,repeat,xs_ct)
            return xs
    functions['ncrystal_crosssection_nonoriented'] = ncrystal_crosssection_nonoriented

    _raw_domain = _wrap('ncrystal_domain',None,(ncrystal_process_t,_dblp,_dblp),hide=True)
    def ncrystal_domain(proc):
        a,b = _dbl(),_dbl()
        _raw_domain(proc,a,b)
        return (a.value,b.value)
    functions['ncrystal_domain'] = ncrystal_domain

    _raw_gs_no = _wrap('ncrystal_genscatter_nonoriented',None,(ncrystal_scatter_t,_dbl,_dblp,_dblp),hide=True)
    _raw_gs_no_many = _wrap('ncrystal_genscatter_nonoriented_many',None,(ncrystal_scatter_t,_dblp,ctypes.c_ulong,
                                                                         ctypes.c_ulong,_dblp,_dblp),hide=True)
    def ncrystal_genscatter_nonoriented(scat,ekin,repeat=None):
        many = _prepare_many(ekin,repeat)
        if many is None:
            angle,de = _dbl(),_dbl()
            _raw_gs_no(scat,ekin,angle,de)
            return angle.value,de.value
        else:
            ekin_ct,n_ekin,repeat,ekin_nparr = many
            angle, angle_ct = _create_numpy_array(n_ekin*repeat)
            de, de_ct = _create_numpy_array(n_ekin*repeat)
            _raw_gs_no_many(scat,ekin_ct,n_ekin,repeat,angle_ct,de_ct)
            return angle,de
    functions['ncrystal_genscatter_nonoriented'] = ncrystal_genscatter_nonoriented

    _raw_xs = _wrap('ncrystal_crosssection',None,(ncrystal_process_t,_dbl,_dbl*3,_dblp),hide=True)
    def ncrystal_crosssection( proc, ekin, direction):
        res = _dbl()
        cdir = (_dbl * 3)(*direction)
        _raw_xs(proc,ekin,cdir,res)
        return res.value
    functions['ncrystal_crosssection'] = ncrystal_crosssection

    _raw_gs = _wrap('ncrystal_genscatter',None,(ncrystal_scatter_t,_dbl,_dbl*3,_dbl*3,_dblp),hide=True)
    def ncrystal_genscatter(scat, ekin, direction):
        cdir = (_dbl * 3)(*direction)
        res_dir = (_dbl * 3)(0,0,0)
        res_de = _dbl()
        _raw_gs(scat,ekin,cdir,res_dir,res_de)
        return (res_dir[0],res_dir[1],res_dir[2]),res_de.value
    functions['ncrystal_genscatter']=ncrystal_genscatter

    _wrap('ncrystal_create_info',ncrystal_info_t,(_cstr,))
    _wrap('ncrystal_create_scatter',ncrystal_scatter_t,(_cstr,))
    _wrap('ncrystal_create_absorption',ncrystal_absorption_t,(_cstr,))
    _raw_save_rng = _wrap('ncrystal_save_randgen',None,tuple(),hide=True)
    _raw_restore_rng = _wrap('ncrystal_restore_randgen',None,tuple(),hide=True)
    _wrap('ncrystal_setbuiltinrandgen',None,tuple())

    _RANDGENFCTTYPE = ctypes.CFUNCTYPE( _dbl )
    _raw_setrand    = _wrap('ncrystal_setrandgen',None,(_RANDGENFCTTYPE,),hide=True)
    def ncrystal_setrandgen(randfct):
        #Set random function, keeping references as needed and casting None to a null-ptr.
        if not randfct:
            keepalive=(None,ctypes.cast(None, _RANDGENFCTTYPE))
        else:
            keepalive=(randfct,_RANDGENFCTTYPE(randfct))#keep refs!
        _globalstates['current_rng']=keepalive
        _raw_setrand(keepalive[1])
    functions['ncrystal_setrandgen'] = ncrystal_setrandgen
    def ncrystal_save_randgen():
        import copy
        _globalstates['saved_rng']=copy.copy(_globalstates.get('current_rng',None))
        _raw_save_rng()
    def ncrystal_restore_randgen():
        import copy
        _globalstates['current_rng']=copy.copy(_globalstates.get('saved_rng',None))
        _raw_restore_rng()
    functions['ncrystal_save_randgen'] = ncrystal_save_randgen
    functions['ncrystal_restore_randgen'] = ncrystal_restore_randgen

    _wrap('ncrystal_decodecfg_packfact',_dbl,(_cstr,))
    _wrap('ncrystal_decodecfg_vdoslux',_uint,(_cstr,))
    _wrap('ncrystal_clear_info_caches',None,tuple())
    _wrap('ncrystal_disable_caching',None,tuple())
    _wrap('ncrystal_enable_caching',None,tuple())
    _wrap('ncrystal_clear_factory_registry',None,tuple())
    _wrap('ncrystal_has_factory',_int,(_cstr,))
    _wrap('ncrystal_clear_caches',None,tuple())

    return functions

_rawfct = _load(_find_nclib())

def decodecfg_packfact(cfgstr):
    return float(_rawfct['ncrystal_decodecfg_packfact'](_str2cstr(cfgstr)))

def decodecfg_vdoslux(cfgstr):
    return int(_rawfct['ncrystal_decodecfg_vdoslux'](_str2cstr(cfgstr)))

def createVDOSDebye(debye_temperature):
    _ensure_numpy()
    #NB: Must keep function exactly synchronised with createVDOSDebye function
    #in .cc src (although leaving out temperature,boundXS,elementMassAMU args
    #here):
    debye_energy = constant_boltzmann*debye_temperature;
    vdos_egrid = _np.linspace(0.5*debye_energy,debye_energy,20);
    scale = 1.0 / (debye_energy*debye_energy);
    vdos_density = scale * (vdos_egrid**2)
    #Actual returned egrid should contain only first and last value:
    return (_np.asarray([vdos_egrid[0],vdos_egrid[-1]]) ,vdos_density)

class RCBase:
    """Base class for all NCrystal objects"""
    def __init__(self, rawobj):
        """internal usage only"""
        self._rawobj = rawobj
        #do not ref here, since ncrystal_create_xxx functions in C-interface already did so.
        self._rawunref = _rawfct['ncrystal_unref']#keep fct reference
    def __del__(self):
        if hasattr(self,'_rawunref') and self._rawunref:
            self._rawunref(self._rawobj)
    def refCount(self):
        """Access reference count of wrapped C++ object"""
        return _rawfct['ncrystal_refcount'](self._rawobj)

def nc_assert(b,msg=""):
    """Assertion which throws NCLogicError on failure"""
    if not bool(b):
        raise NCLogicError(msg if msg else 'assertion failed')

class Info(RCBase):
    """Class representing information about a given material"""
    def __init__(self, cfgstr):
        """create Info object based on cfg-string (same as using createInfo(cfgstr))"""
        rawobj = _rawfct['ncrystal_create_info'](_str2cstr(cfgstr))
        super(Info, self).__init__(rawobj)
        self.__dyninfo=None

    def dump(self):
        """Dump contained information to standard output"""
        sys.stdout.flush()
        sys.stderr.flush()
        _rawfct['ncrystal_dump'](self._rawobj)

    def hasTemperature(self):
        return _rawfct['ncrystal_info_gettemperature'](self._rawobj)>-1
    def getTemperature(self):
        t=_rawfct['ncrystal_info_gettemperature'](self._rawobj)
        nc_assert(t>-1)
        return t

    def hasGlobalDebyeTemperature(self):
        return _rawfct['ncrystal_info_getglobaldebyetemp'](self._rawobj)>-1
    def getGlobalDebyeTemperature(self):
        t=_rawfct['ncrystal_info_getglobaldebyetemp'](self._rawobj)
        nc_assert(t>-1)
        return t

    def hasDensity(self):
        return _rawfct['ncrystal_info_getdensity'](self._rawobj)>-1
    def getDensity(self):
        t=_rawfct['ncrystal_info_getdensity'](self._rawobj)
        nc_assert(t>-1)
        return t

    def hasNumberDensity(self):
        return _rawfct['ncrystal_info_getnumberdensity'](self._rawobj)>-1
    def getNumberDensity(self):
        t=_rawfct['ncrystal_info_getnumberdensity'](self._rawobj)
        nc_assert(t>-1)
        return t

    def hasXSectAbsorption(self):
        return _rawfct['ncrystal_info_getxsectabsorption'](self._rawobj)>-1
    def getXSectAbsorption(self):
        t=_rawfct['ncrystal_info_getxsectabsorption'](self._rawobj)
        nc_assert(t>-1)
        return t

    def hasXSectFree(self):
        return _rawfct['ncrystal_info_getxsectfree'](self._rawobj)>-1
    def getXSectFree(self):
        t=_rawfct['ncrystal_info_getxsectfree'](self._rawobj)
        nc_assert(t>-1)
        return t

    def hasStructureInfo(self):
        return bool(_rawfct['ncrystal_info_getstructure'](self._rawobj))
    def getStructureInfo(self):
        d=_rawfct['ncrystal_info_getstructure'](self._rawobj)
        nc_assert(d)
        return d

    def hasHKLInfo(self):
        return bool(_rawfct['ncrystal_info_nhkl'](self._rawobj)>-1)
    def nHKL(self):
        return int(_rawfct['ncrystal_info_nhkl'](self._rawobj))
    def hklDLower(self):
        return float(_rawfct['ncrystal_info_hkl_dlower'](self._rawobj))
    def hklDUpper(self):
        return float(_rawfct['ncrystal_info_hkl_dupper'](self._rawobj))
    def hklList(self):
        """Iterator over HKL info, yielding tuples in the format
        (h,k,l,multiplicity,dspacing,fsquared)"""
        nc_assert(self.hasHKLInfo())
        h,k,l,mult,dsp,fsq = _rawfct['ncrystal_info_gethkl_setuppars']()
        for idx in range(self.nHKL()):
            _rawfct['ncrystal_info_gethkl'](self._rawobj,idx,h,k,l,mult,dsp,fsq)
            yield h.value,k.value,l.value,mult.value,dsp.value,fsq.value
    def dspacingFromHKL(self, h, k, l):
        return float(_rawfct['ncrystal_info_dspacing_from_hkl'](self._rawobj,h,k,l))

    class DynamicInfo:
        """Class representing dynamic information (related to inelastic scattering)
           about a given element"""
        def __init__(self,fr,en,tt,key):
            self.__fraction, self.__elementName, self._key, self.__tt = fr,en,key,tt
        @property
        def fraction(self):
            """Element fraction in material (all fractions must add up to unity)"""
            return self.__fraction
        @property
        def temperature(self):
            """Material temperature (same value as on associated Info object)"""
            return self.__tt
        @property
        def elementName(self):
            """Name of element"""
            return self.__elementName
        def _np(self):
            _ensure_numpy()
            return _np
        def _copy_cptr_2_nparray(self,cptr,n):
            np = self._np()
            return np.copy(np.ctypeslib.as_array(cptr, shape=(n,)))

    class DI_Sterile(DynamicInfo):
        """Class indicating elements for which inelastic neutron scattering is absent
           or disabled."""
        pass

    class DI_FreeGas(DynamicInfo):
        """Class indicating elements for which inelastic neutron scattering should be
           modelled as scattering on a free gas."""
        pass

    class DI_ScatKnl(DynamicInfo):
        """Base class indicating elements for which inelastic neutron scattering will
           be, directly or indirectly, described by a scattering kernel,
           S(alpha,beta). This is an abstract class, and derived classes provide
           actual access to the kernels.
        """

        def __init__(self,fr,en,tt,key):
            super(Info.DI_ScatKnl, self).__init__(fr,en,tt,key)
            self.__lastknl,self.__lastvdoslux = None,None

        def _loadKernel( self, vdoslux = 3 ):
            assert isinstance(vdoslux,int) and 0<=vdoslux<=5
            if self.__lastvdoslux != vdoslux:
                sugEmax,ne,na,nb,eptr,aptr,bptr,sabptr = _rawfct['ncrystal_dyninfo_extract_scatknl'](self._key,vdoslux)
                self.__lastvdoslux = vdoslux
                res={}
                assert ne>=0
                res['suggestedEmax'] = float(sugEmax)
                res['egrid'] = self._copy_cptr_2_nparray(eptr,ne) if ne > 0 else self._np().zeros(0)
                assert na>1 and nb>1
                res['alpha'] = self._copy_cptr_2_nparray(aptr,na)
                res['beta']  = self._copy_cptr_2_nparray(bptr,nb)
                res['sab']   = self._copy_cptr_2_nparray(sabptr,na*nb)
                self.__lastknl = res
            assert self.__lastknl is not None
            return self.__lastknl

    class DI_ScatKnlDirect(DI_ScatKnl):
        """Pre-calculated scattering kernel which at most needs a (hidden) conversion to
           S(alpha,beta) format before it is available."""
        def loadKernel( self ):
            """Prepares and returns the scattering kernel in S(alpha,beta) format"""
            return self._loadKernel(vdoslux=3)#vdoslux value not actually used

    class DI_VDOS(DI_ScatKnl):
        """Solid state material with a phonon spectrum in the form of a Vibrational
        Density Of State (VDOS) parameterisation. This can be expanded into a
        full scattering kernel. How luxurious this expansion will be is
        controlled by an optional vdoslux parameter in the loadKernel call (must
        be integer from 0 to 5)
        """
        def __init__(self,fr,en,tt,key):
            super(Info.DI_VDOS, self).__init__(fr,en,tt,key)
            self.__vdosdata = None
            self.__vdosegrid_expanded = None
            self.__vdosorig = None

        def vdosData(self):
            """Access the VDOS as ([egrid_min,egrid_max],vdos_density)"""
            if self.__vdosdata is None:
                emin,emax,nd,dptr = _rawfct['ncrystal_dyninfo_extract_vdos'](self._key)
                vdos_egrid = (emin,emax)
                vdos_density = self._copy_cptr_2_nparray(dptr,nd)
                self.__vdosdata = (vdos_egrid,vdos_density)
            return self.__vdosdata

        def __loadVDOSOrig(self):
            if self.__vdosorig is None:
                neg,egptr,nds,dsptr = _rawfct['ncrystal_dyninfo_extract_vdos_input'](self._key)
                self.__vdosorig = ( self._copy_cptr_2_nparray(egptr,neg),
                                    self._copy_cptr_2_nparray(dsptr,nds) )
            return self.__vdosorig

        def vdosOrigEgrid(self):
            """Access the original un-regularised VDOS energy grid"""
            return self.__loadVDOSOrig()[0]

        def vdosOrigDensity(self):
            """Access the original un-regularised VDOS energy grid"""
            return self.__loadVDOSOrig()[1]

        @property
        def vdos_egrid(self):
            """Access the VDOS energy grid as [egrid_min,egrid_max]"""
            return self.vdosData()[0]

        @property
        def vdos_egrid_expanded(self):
            """Access the egrid expanded into all actual egrid points"""
            if self.__vdosegrid_expanded is None:
                _ = self.vdosData()
                self.__vdosegrid_expanded = self._np().linspace(_[0][0],_[0][1],len(_[1]))
            return self.__vdosegrid_expanded

        @property
        def vdos_density(self):
            """Access the VDOS density array"""
            return self.vdosData()[1]

        def loadKernel( self, vdoslux = 3 ):
            """Converts VDOS to S(alpha,beta) kernel with a luxury level given by the vdoslux parameter."""
            return self._loadKernel(vdoslux=vdoslux)

    class DI_VDOSDebye(DI_ScatKnl):
        """Similarly to DI_VDOS, but instead of using a phonon VDOS spectrum provided
           externally, an idealised spectrum is used for lack of better
           options. This spectrum is based on the Debye Model, in which the
           spectrum rises quadratically with phonon energy below a cutoff value,
           kT, where T is the Debye temperature
        """

        def __init__(self,fr,en,tt,key):
            super(Info.DI_VDOSDebye, self).__init__(fr,en,tt,key)
            self.__vdosdata = None
            self.__debyetemp = None
            self.__vdosegrid_expanded = None

        def vdosData(self):
            """Access the idealised VDOS as ([egrid_min,egrid_max],vdos_density)"""
            if self.__vdosdata is None:
                self.__vdosdata = createVDOSDebye(self.debyeTemperature())
            return self.__vdosdata

        def debyeTemperature(self):
            #The Debye temperature of the element.
            if self.__debyetemp is None:
                self.__debyetemp = _rawfct['ncrystal_dyninfo_extract_vdosdebye'](self._key)
            return self.__debyetemp

        @property
        def vdos_egrid(self):
            """Access the VDOS energy grid as [egrid_min,egrid_max]"""
            return self.vdosData()[0]

        @property
        def vdos_egrid_expanded(self):
            """Access the egrid expanded into all actual egrid points"""
            if self.__vdosegrid_expanded is None:
                _ = self.vdosData()
                self.__vdosegrid_expanded = self._np().linspace(_[0][0],_[0][1],len(_[1]))
            return self.__vdosegrid_expanded

        @property
        def vdos_density(self):
            """Access the VDOS density array"""
            return self.vdosData()[1]

        def loadKernel( self, vdoslux = 3 ):
            """Converts VDOS to S(alpha,beta) kernel with a luxury level given by the
               vdoslux parameter, which is similar to the vdoslux parameter used
               in DI_VDOS. Notice that the vdoslux parameter specified here on
               DI_VDOSDebye will be reduced internally by 3 (but not less than
               0), since the Debye model is anyway only a crude approximation
               and it accordingly does not need the same level of precise
               treatment as a full externally specified VDOS.
            """
            return self._loadKernel(vdoslux=vdoslux)

    def hasDynamicInfo(self):
        """Whether or not dynamic information for each element is present"""
        return int(_rawfct['ncrystal_info_ndyninfo'](self._rawobj))>0 if self.__dyninfo is None else bool(self.__dyninfo)

    def getDynamicInfoList(self):
        """Get list of DynamicInfo objects (if available). One for each element"""
        if self.__dyninfo is None:
            self.__dyninfo = []
            for idx in range(int(_rawfct['ncrystal_info_ndyninfo'](self._rawobj))):
                key = (self._rawobj,idx)
                fr,tt,en,ditype = _rawfct['ncrystal_dyninfo_base'](key)
                args=(fr,en,tt,key)
                if ditype==0:
                    di = Info.DI_Sterile(*args)
                elif ditype==1:
                    di = Info.DI_FreeGas(*args)
                elif ditype==2:
                    di = Info.DI_ScatKnlDirect(*args)
                elif ditype==3:
                    di = Info.DI_VDOS(*args)
                elif ditype==4:
                    di = Info.DI_VDOSDebye(*args)
                else:
                    raise AssertionError('Unknown DynInfo type id (%i)'%ditype.value)
                self.__dyninfo += [ di ]
        return self.__dyninfo

    @property
    def dyninfos(self):
        """List of DynamicInfo objects (empty if not available). One for each element"""
        return self.getDynamicInfoList()

class CalcBase(RCBase):
    """Base class for all calculators"""
    def getCalcName(self):
        """Calculator name"""
        return _cstr2str(_rawfct['ncrystal_name'](self._rawobj))
    @property
    def name(self):
        """Calculator name as property"""
        return self.getCalcName()

class Process(CalcBase):
    """Base class for calculations of processes in materials.

    Note that kinetic energies are in electronvolt and direction vectors are
    tuples of 3 numbers.

    """
    def domain(self):
        """Domain where process has non-vanishing cross-section.

        Returns the domain as (ekin_low,ekin_high). Outside this range of
        neutron kinetic energy, the process can be assumed to have vanishing
        cross-sections. Thus, processes present at all energies will return
        (0.0,infinity).

        """
        return _rawfct['ncrystal_domain'](self._rawobj)
    def isNonOriented(self):
        """opposite of isOriented()"""
        return bool(_rawfct['ncrystal_isnonoriented'](self._rawobj))
    def isOriented(self):
        """Check if process is oriented and results depend on the incident direction of the neutron"""
        return not self.isNonOriented()
    def crossSection( self, ekin, direction ):
        """Access cross-sections."""
        return _rawfct['ncrystal_crosssection'](self._rawobj,ekin, direction)
    def crossSectionNonOriented( self, ekin, repeat = None ):
        """Access cross-sections (should not be called for oriented processes).

        For efficiency it is possible to provide the ekin parameter as a numpy
        array of numbers and get a corresponding array of cross-sections
        back. Likewise, the repeat parameter can be set to a positive number,
        causing the ekin value(s) to be reused that many times and a numpy array
        with results returned.

        """
        return _rawfct['ncrystal_crosssection_nonoriented'](self._rawobj,ekin,repeat)

class Absorption(Process):
    """Base class for calculations of absorption in materials"""

    def __init__(self, cfgstr):
        """create Absorption object based on cfg-string (same as using createAbsorption(cfgstr))"""
        rawobj_abs = _rawfct['ncrystal_create_absorption'](_str2cstr(cfgstr))
        rawobj_proc = _rawfct['ncrystal_cast_abs2proc'](rawobj_abs)
        super(Absorption, self).__init__(rawobj_proc)

class Scatter(Process):
    """Base class for calculations of scattering in materials.

    Note that kinetic energies are in electronvolt and direction vectors are
    tuples of 3 numbers.

    """

    def __init__(self, cfgstr):
        """create Scatter object based on cfg-string (same as using createScatter(cfgstr))"""
        self._rawobj_scat = _rawfct['ncrystal_create_scatter'](_str2cstr(cfgstr))
        rawobj_proc = _rawfct['ncrystal_cast_scat2proc'](self._rawobj_scat)
        super(Scatter, self).__init__(rawobj_proc)
    def generateScattering( self, ekin, direction ):
        """Randomly generate scatterings.

        Assuming a scattering took place, generate energy transfer (delta_ekin)
        and new neutron direction based on current kinetic energy and direction
        and return tuple(new_direction,delta_ekin).

        """
        return _rawfct['ncrystal_genscatter'](self._rawobj_scat,ekin,direction)
    def generateScatteringNonOriented( self, ekin, repeat = None ):
        """Randomly generate scatterings (should not be called for oriented processes).

        Assuming a scattering took place, generate energy transfer (delta_ekin)
        and scatter angle in radians and return tuple(scatter_angle,delta_ekin)
        (this method should not be invoked on oriented processes).  For
        efficiency it is possible to provide the ekin parameter as a numpy array
        of numbers and get corresponding arrays of angles and energy transfers
        back. Likewise, the repeat parameter can be set to a positive number,
        causing the ekin value(s) to be reused that many times and numpy arrays
        with results returned.

        """
        return _rawfct['ncrystal_genscatter_nonoriented'](self._rawobj_scat,ekin,repeat)

def createInfo(cfgstr):
    """Construct Info object based on provided configuration (using available factories)"""
    return Info(cfgstr)

def createScatter(cfgstr):
    """Construct Scatter object based on provided configuration (using available factories)"""
    return Scatter(cfgstr)

def createAbsorption(cfgstr):
    """Construct Absorption object based on provided configuration (using available factories)"""
    return Absorption(cfgstr)

#numpy compatible wl2ekin and ekin2wl
_c_wl2ekin = float(_rawfct['ncrystal_wl2ekin'](1.0))
_c_ekin2wl = float(_rawfct['ncrystal_ekin2wl'](1.0))

def wl2ekin(wl):
    """Convert neutron wavelength in Angstrom to kinetic energy in electronvolt"""
    if _np and hasattr(wl,'__len__'):
        #reciprocals without zero division:
        wlnonzero = wl != 0.0
        wlinv = 1.0 / _np.where( wlnonzero, wl, 1.0)#fallback 1.0 wont be used
        return _c_wl2ekin * _np.square(_np.where( wlnonzero, wlinv, _np.inf))
    else:
        return _rawfct['ncrystal_wl2ekin'](wl)

def ekin2wl(ekin):
    """Convert neutron kinetic energy in electronvolt to wavelength in Angstrom"""
    if _np and hasattr(ekin,'__len__'):
        #reciprocals without zero division:
        ekinnonzero = ekin != 0.0
        ekininv = 1.0 / _np.where( ekinnonzero, ekin, 1.0)#fallback 1.0 wont be used
        return _c_ekin2wl * _np.sqrt(_np.where( ekinnonzero, ekininv, _np.inf))
    else:
        return _rawfct['ncrystal_ekin2wl'](ekin)

def clearInfoCaches():
    """Clear all cached Info objects held by factory infrastructure. NB: included in clearCaches()"""
    _rawfct['ncrystal_clear_info_caches']()
def clearCaches():
    """Clear various caches"""
    _rawfct['ncrystal_clear_caches']()
def disableCaching():
    """Disable caching of Info objects in factory infrastructure"""
    _rawfct['ncrystal_disable_caching']()
def enableCaching():
    """Enable caching of Info objects in factory infrastructure"""
    _rawfct['ncrystal_enable_caching']()
def clearFactoryRegistry():
    """Clear all registered factories"""
    _rawfct['ncrystal_clear_factory_registry']()
def hasFactory(name):
    """Check if a factory of a given name exists"""
    return bool(_rawfct['ncrystal_has_factory'](_str2cstr(name)))

#Helper function, for scripts creating ncmat files:
def formatVectorForNCMAT(name,values):
    """Utility function for help in python scripts composing .ncmat files,
       transforming an array of of values into a properly formatted text string,
       with word-wrapping, usage of <val>r<n> syntax, etc. Returns list of lines
       (strings) for .ncmat files.
    """

    v,res,l,indent,efmt_prev,efmt_nrepeat=values.flatten(),'  %s'%name,'','',None,0
    if not len(v):
        return res
    ilast,nv=len(v)-1,len(v)
    def _fmtnum(num):
        _ = '%g'%num if num else '0'#avoid 0.0, -0, etc.
        if _.startswith('0.'):
            _=_[1:]
        return _
    i=0
    while i<nv:
        fmt_vi=_fmtnum(v[i])
        #check if is repeated:
        irepeat=i
        while irepeat+1<nv:
            if _fmtnum(v[irepeat+1])==fmt_vi:
                irepeat+=1
            else:
                break
        #Write:
        s = ' %sr%i'%(fmt_vi,1+irepeat-i) if irepeat>i else ' %s'%fmt_vi
        l+=(s if l else (indent+s))
        i=irepeat+1#advance
        if i>=nv or len(l)>80:
            #Flush line
            res += '%s\n'%l
            l,indent='',' '
    return res

#Accept custom random generator:
def setDefaultRandomGenerator(rg):
    """Set the default random generator for CalcBase classes.

    Note that this can only changes the random generator for those CalcBase
    instances that did not already use random numbers). Default generator when
    using the NCrystal python interface is the scientifically sound
    random.random stream from the python standard library (a Mersenne Twister).
    """
    _rawfct['ncrystal_setrandgen'](rg)

#From python we default to the usual python random stream. It can be cleared
#with setDefaultRandomGenerator(None) and its seed controlled with
#random.seed(somevalue)
try:
    import random
    setDefaultRandomGenerator(random.random)
except ImportError:
    pass

class RandomCtxMgr:
    """Context manager which can be used to temporarily change the random stream used by NCrystal.

    This is mainly intended for internal usage in the test() function, and does
    not support multiple instances:

    """
    def __init__(self,randfct):
        self._r = randfct
    def __enter__(self):
        _rawfct['ncrystal_save_randgen']()
        setDefaultRandomGenerator(self._r)
    def __exit__(self,*args):
        _rawfct['ncrystal_restore_randgen']()

__natelemdb={}
def naturalElementDB(z):
    """Access internal database with data for natural element with z
       protons. Returns dictionary with data (empty dictionary if not available)."""
    global __natelemdb
    z=int(z)
    _=__natelemdb.get(z,None)
    if _ is not None:
        return _
    _=_rawfct['ncrystal_natelemdata'](z)
    __natelemdb[z] = _
    return _

__natelemdb_n2z={}
def naturalElementDBByName(atom_name):
    """Same as naturalElementDB but access via name (H, He, Li, ...)."""
    #First check if already loaded:
    for e in __natelemdb:
        _=e.get('name',None)
        if _==atom_name:
            return _
    #The hard way, populate DB until found:
    for z in range(1,121):
        _=naturalElementDB(z)
        if _.get('mass_amu',None)==atom_name:
            return _
    #Not found:
    return {}

def test():
    """Quick test that NCrystal works as expected in the current installation."""
    with RandomCtxMgr(None):
        _actualtest()
    print("Tests completed succesfully")

def _actualtest():
    def require(b):
        if not b:
            raise RuntimeError('check failed')
    def flteq(a,b,rtol=1.0e-6,atol=1.0e-6):
        return abs(a-b) <= 0.5 * rtol * (abs(a) + abs(b)) + atol
    def require_flteq(a,b):
        if not flteq(a,b):
            raise RuntimeError('check failed (%.16g != %.16g, diff %g)'%(a,b,a-b))
        return True

    al = createInfo('Al_sg225.ncmat;dcutoff=1.4')
    require(hasFactory('stdncmat'))
    require(al.hasTemperature() and require_flteq(al.getTemperature(),293.15))
    require(al.hasXSectFree() and require_flteq(al.getXSectFree(),1.39667))
    require(al.hasXSectAbsorption() and require_flteq(al.getXSectAbsorption(),0.231))
    require(al.hasDensity() and require_flteq(al.getDensity(),2.69864547673))
    require(al.hasNumberDensity() and require_flteq(al.getNumberDensity(),0.06023238256131625))
    require(not al.hasGlobalDebyeTemperature())

    require(al.hasStructureInfo())
    si=al.getStructureInfo()
    require( si['spacegroup'] == 225 )
    require_flteq(si['a'],4.04958)
    require_flteq(si['b'],4.04958)
    require_flteq(si['c'],4.04958)
    require( si['alpha'] == 90.0 )
    require( si['beta'] == 90.0 )
    require( si['gamma'] == 90.0 )
    require( si['n_atoms'] == 4 )
    require_flteq(si['volume'],66.4094599932)
    require( al.hasHKLInfo() )
    require( al.nHKL() == 3 )
    require_flteq(al.hklDLower(),1.4)
    require( al.hklDUpper() > 1e36 )
    expected_hkl = { 0  : (1, -1, -1, 8, 2.3380261031049243, 1.77210688280773),
                     1  : (0, 0, 2, 6, 2.02479, 1.730418273564928),
                     2  : (0, 2, -2, 12, 1.4317427394787094, 1.573243021457801) }
    for idx,hkl in enumerate(al.hklList()):
        h,k,l,mult,dsp,fsq = hkl
        require(idx<len(expected_hkl))
        e = expected_hkl[idx]
        require( list(e)[0:4] == [h,k,l,mult] )
        require_flteq(dsp, e[4])
        require_flteq(fsq, e[5])
    #TODO for NC2: Missing is atominfo, bkgd cross-sections, ...

    alpc = createScatter('Al_sg225.ncmat;dcutoff=1.4;incoh_elas=0;inelas=0')
    require( alpc.name == 'PCBragg' )
    require( isinstance(alpc.name,str) )
    require( alpc.refCount() in (1,2) and type(alpc.refCount()) == int )
    require( alpc.isNonOriented() )
    require_flteq(1.631341646154576,alpc.crossSectionNonOriented(wl2ekin(4.0)) )
    require_flteq(1.631341646154576,alpc.crossSection(wl2ekin(4.0),(1,0,0)))
    require( alpc.crossSectionNonOriented(wl2ekin(5.0)) == 0.0 )

    #use NCrystal's own builtin rng for guaranteed reproducibility (set explicitly to avoid warning):
    _rawfct['ncrystal_setbuiltinrandgen']()
    for i in range(60):
        alpc.generateScatteringNonOriented(wl2ekin(4.0))

    alpc = createScatter('Al_sg225.ncmat;dcutoff=1.4')
    require( alpc.name == 'ScatterComp' )
    expected = [ ( 2.0527318521221005, 0.0 ),
                 ( 2.8283092712311073, 0.0 ),
                 ( 2.8283092712311073, 0.0 ),
                 ( 2.0527318521221, 0.0 ),
                 ( 2.0527318521221, 0.0 ),
                 ( 2.0527318521221005, 0.0 ),
                 ( 1.877750222153273, 0.02107623981594055 ),
                 ( 2.8283092712311073, 0.0 ),
                 ( 2.0527318521221, 0.0 ),
                 ( 2.8283092712311073, 0.0 )]

    if _np is None:
        ang,de=[],[]
        for i in range(10):
            _ang,_de=alpc.generateScatteringNonOriented(wl2ekin(4.0))
            ang += [_ang]
            de += [_de]
    else:
        ang,de = alpc.generateScatteringNonOriented(wl2ekin(4.0),repeat=10)

    for i in range(10):
        #print ( f'    ( {ang[i]}, {de[i]} ),');continue
        require_flteq(ang[i],expected[i][0])
        require_flteq(de[i],expected[i][1])

    expected = [ ( (0.6341262223410173, 0.42508486262350836, 0.6458999873880349), 0.0 ),
                 ( (0.6341262223410173, 0.6030479945760753, 0.48395976111375677), 0.0 ),
                 ( (0.6341262223410173, -0.5323426933952337, 0.5607987080300906), 0.0 ),
                 ( (0.5121682964546899, -0.8208697888708928, 0.25269829011245154), 0.0 ),
                 ( (0.6341262223410173, 0.5169272374773588, -0.5750392728271148), 0.0 ),
                 ( (0.6341262223410173, -0.4591420535083147, 0.6221515159827856), 0.0 ),
                 ( (-0.699307619791641, -0.7068851945247946, 0.10621758170375115), 0.0051568432579783925 ),
                 ( (0.02433659290937973, 0.8727814258082862, 0.4875041671715413), 0.0 ),
                 ( (0.6341262223410173, 0.3115365900742872, -0.7076926502263509), 0.0 ),
                 ( (0.6341262223410173, 0.14449943153842737, -0.7596076937634202), 0.0 )]

    for i in range(10):
        outdir,de = alpc.generateScattering(wl2ekin(2.0),(1.0,0.0,0.0))
        #print ( f'    ( {outdir}, {de} ),');continue
        require_flteq(de,expected[i][1])
        require_flteq(outdir[0],expected[i][0][0])
        require_flteq(outdir[1],expected[i][0][1])
        require_flteq(outdir[2],expected[i][0][2])
    gesc = createScatter("""Ge_sg227.ncmat;dcutoff=0.5;mos=40.0arcsec
                            ;dir1=@crys_hkl:5,1,1@lab:0,0,1
                            ;dir2=@crys_hkl:0,-1,1@lab:0,1,0""")
    require_flteq(587.7399848353572,gesc.crossSection(wl2ekin(1.540),( 0., 1., 1. )))
    require_flteq(1.662542957156909,gesc.crossSection(wl2ekin(1.540),( 1., 1., 0. )))

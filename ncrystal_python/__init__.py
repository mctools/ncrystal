#!/usr/bin/env python3
"""Python module for using the NCrystal library for thermal neutron transport in crystals and other materials.

Please find more information about NCrystal at the website:

   https://mctools.github.io/ncrystal/

In particular, a small example using the NCrystal python module can be found at:

   https://github.com/mctools/ncrystal/blob/master/examples/ncrystal_example_py

A substantial effort went into developing NCrystal. If you use it for your work,
we would appreciate it if you would use the following reference in your work:

  X.-X. Cai and T. Kittelmann, NCrystal: A library for thermal neutron
  transport, Computer Physics Communications 246 (2020) 106851,
  https://doi.org/10.1016/j.cpc.2019.07.015

For work benefitting from our inelastic physics, we furthermore request that you
additionally also use the following reference in your work:

  X.-X. Cai, T. Kittelmann, et. al., "Rejection-based sampling of inelastic
  neutron scattering", Journal of Computational Physics 380 (2019) 400-407,
  https://doi.org/10.1016/j.jcp.2018.11.043

For detailed usage conditions and licensing of this open source project, see:

   https://github.com/mctools/ncrystal/blob/master/NOTICE
   https://github.com/mctools/ncrystal/blob/master/LICENSE

"""

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2022 NCrystal developers                                   ##
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

#NB: Synchronize meta-data below with fields in setup.py.in meta data:
__license__ = "Apache 2.0, http://www.apache.org/licenses/LICENSE-2.0"
__version__ = '3.5.1'
__status__ = "Production"
__author__ = "NCrystal developers (Thomas Kittelmann, Xiao Xiao Cai)"
__copyright__ = "Copyright 2015-2022 %s"%__author__
__maintainer__ = __author__
__email__ = "ncrystal-developers@cern.ch"
#Only put the few most important items in __all__, to prevent cluttering on
#wildcard imports. Specifically this is the exceptions, the most important API
#classes, the factory functions, and the constants:
__all__ = [ 'NCException','NCFileNotFound','NCDataLoadError','NCMissingInfo','NCCalcError',
            'NCLogicError','NCBadInput','RCBase','TextData','Info','Process',
            'Absorption','Scatter','AtomData','FileListEntry','createTextData',
            'createInfo','createScatter','createScatterIndependentRNG','createAbsorption',
            'constant_c','constant_dalton2kg','constant_dalton2eVc2','constant_avogadro',
            'constant_boltzmann','const_neutron_mass_amu','constant_planck']

#Place f-string here to catch python <3.6 in a more obvious way than a syntax error below:
f'NCrystal does not work with Python2 (or Python3 < v3.6)'
_minpyversion=(3,6,0)

import sys
pyversion = sys.version_info[0:3]
if pyversion < _minpyversion:
    raise SystemExit('Unsupported python version %i.%i.%i detected (needs %i.%i.%i or later).'%(pyversion+_minpyversion))

import numbers
import pathlib
import os
import copy
import ctypes
import weakref
import enum
import json
import collections

###################################
#Convert cstr<->str:

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
_kPi        = 3.1415926535897932384626433832795028841971694
_k2Pi       = 6.2831853071795864769252867665590057683943388
_k4Pidiv100 = 0.125663706143591729538505735331180115367886776
_k4PiSq     = 39.4784176043574344753379639995046045412547976
standard_comp_types = ('coh_elas','incoh_elas','inelas','sans')#for client code checking all types of components.

def _find_nclib():

    #If NCRYSTAL_LIB env var is set, we try that and only that:
    override=os.environ.get('NCRYSTAL_LIB',None)
    if override:
        override = pathlib.Path(override)
        if not override.exists() or override.is_dir():
            raise NCFileNotFound('NCRYSTAL_LIB environment variable is set but does not point to an actual file.')
        return override.absolute().resolve()

    try:
        if __name__ != '__main__':
            #normal import
            from . import _nclibpath
        else:
            #work if running as script:
            sys.path.insert(0,str(pathlib.Path(__file__).absolute().parent))
            import _nclibpath
            sys.path.pop(0)
    except ImportError:
        raise NCFileNotFound('Autogenerated _nclibpath.py module not found (it should have been generated'
                             +' during installation). In this case you must set the environment variable'
                             +' NCRYSTAL_LIB to point at the compiled NCrystal library.')
    _ = pathlib.Path(_nclibpath.liblocation)
    if not _.is_absolute():
        _ = (pathlib.Path(__file__).absolute().parent / _)
    if not _.exists() or _.is_dir():
        raise NCFileNotFound('Autogenerated _nclibpath.py module was found but no file exists in the indicated'
                             +' library location (%s). Either reinstall NCrystal or try to use the environment variable'%_
                             +' NCRYSTAL_LIB to point at the compiled NCrystal library.')
    return _.resolve()

try:
    import numpy as _np
except ImportError:
    _np = None
def _ensure_numpy():
    if not _np:
        raise NCException("Numpy not available - array based functionality is unavailable")

_keepalive = []

def _np_linspace(start,stop,num=50):
    """linspace with reproducible endpoint value"""
    _ensure_numpy()
    assert num >= 2
    l = _np.linspace(start,stop,num)
    l[0] = start
    l[-1] = stop
    return l

def _np_geomspace(start,stop,num=50):
    """geomspace with reproducible endpoint value"""
    _ensure_numpy()
    assert num >= 2
    l = _np.geomspace(start,stop,num)
    l[0] = start
    l[-1] = stop
    return l

def _np_logspace(start,stop,num=50):
    """logspace with reproducible endpoint value"""
    _ensure_numpy()
    assert num >= 2
    l = _np.logspace(start,stop,num)
    l[0] = start
    l[-1] = stop
    return l

def _load(nclib_filename):

    _nclib = ctypes.CDLL(nclib_filename)
    _int,_intp,_uint,_uintp,_dbl,_dblp,_cstr,_voidp = (ctypes.c_int, ctypes.POINTER(ctypes.c_int),
                                                       ctypes.c_uint,ctypes.POINTER(ctypes.c_uint), ctypes.c_double,
                                                       ctypes.POINTER(ctypes.c_double), ctypes.c_char_p, ctypes.c_void_p)
    _ulong = ctypes.c_ulong
    _charptr = ctypes.POINTER(ctypes.c_char)

    _cstrp = ctypes.POINTER(_cstr)
    _cstrpp = ctypes.POINTER(_cstrp)
    _dblpp = ctypes.POINTER(_dblp)
    ndarray_to_dblp = lambda a : a.ctypes.data_as(_dblp)
    ndarray_to_uintp = lambda a : a.ctypes.data_as(_uintp)
    ndarray_to_intp = lambda a : a.ctypes.data_as(_intp)

    def _create_numpy_double_array(n):
        _ensure_numpy()
        a=_np.empty(n,dtype=_dbl)
        return a,ndarray_to_dblp(a)

    def _create_numpy_unsigned_array(n):
        _ensure_numpy()
        a=_np.empty(n,dtype=_uint)
        return a,ndarray_to_uintp(a)

    def _create_numpy_int_array(n):
        _ensure_numpy()
        a=_np.empty(n,dtype=_int)
        return a,ndarray_to_intp(a)

    class ncrystal_info_t(ctypes.Structure):
        _fields_ = [('internal', _voidp)]
    class ncrystal_process_t(ctypes.Structure):
        _fields_ = [('internal', _voidp)]
    class ncrystal_scatter_t(ctypes.Structure):
        _fields_ = [('internal', _voidp)]
    class ncrystal_absorption_t(ctypes.Structure):
        _fields_ = [('internal', _voidp)]
    class ncrystal_atomdata_t(ctypes.Structure):
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
        #TODO: Provide line number / file as well?
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
        raise RuntimeError("ERROR: Version mismatch detected between NCrystal python code (v%s)"
                           " and loaded binary"" library (v%s). Control which NCrystal library"
                           " to load with the NCRYSTAL_LIB env var."%(__version__,lib_version))

    _wrap('ncrystal_sethaltonerror',_int,(_int,),hide=True,error_check=False)(False)
    _wrap('ncrystal_setquietonerror',_int,(_int,),hide=True,error_check=False)(True)
    _ncerror       = _wrap('ncrystal_error',_int,tuple(),hide=True,error_check=False)
    _ncerror_msg   = _wrap('ncrystal_lasterror',_cstr,tuple(),hide=True,error_check=False)
    _ncerror_type  = _wrap('ncrystal_lasterrortype',_cstr,tuple(),hide=True,error_check=False)
    _ncerror_clear = _wrap('ncrystal_clearerror',None,tuple(),hide=True,error_check=False)

    _wrap('ncrystal_refcount',_int,(_voidp,),take_ref=True)
    _wrap('ncrystal_valid',_int,(_voidp,),take_ref=True)

    #NB: For ncrystal_unref we use take_ref=False, so RCBase.__del__ can cache
    #the result of ctypes.byref(rawobj). This is needed since the ctypes module
    #might have been unloaded before RCBase.__del__ is called:
    _wrap('ncrystal_unref',None,(_voidp,),take_ref=False)

    _wrap('ncrystal_cast_scat2proc',ncrystal_process_t,(ncrystal_scatter_t,))
    _wrap('ncrystal_cast_abs2proc',ncrystal_process_t,(ncrystal_absorption_t,))

    _wrap('ncrystal_dump',None,(ncrystal_info_t,))
    _wrap('ncrystal_dump_verbose',None,(ncrystal_info_t,_uint))
    _wrap('ncrystal_ekin2wl',_dbl,(_dbl,))
    _wrap('ncrystal_wl2ekin',_dbl,(_dbl,))
    _wrap('ncrystal_isnonoriented',_int,(ncrystal_process_t,))
    _wrap('ncrystal_name',_cstr,(ncrystal_process_t,))

    _wrap('ncrystal_debyetemp2msd',_dbl,(_dbl,_dbl,_dbl))
    _wrap('ncrystal_msd2debyetemp',_dbl,(_dbl,_dbl,_dbl))

    _wrap('ncrystal_create_atomdata_fromdb',ncrystal_atomdata_t,(_uint,_uint))
    _wrap('ncrystal_create_atomdata_fromdbstr',ncrystal_atomdata_t,(_cstr,))

    _raw_atomdb_getn = _wrap('ncrystal_atomdatadb_getnentries',_uint,tuple(), hide=True )
    _raw_atomdb_getall = _wrap('ncrystal_atomdatadb_getallentries',_uint,(_uintp,_uintp), hide=True )
    def atomdb_getall_za():
        n = _raw_atomdb_getn()
        zvals,zvalsptr = _create_numpy_unsigned_array(n)
        avals,avalsptr = _create_numpy_unsigned_array(n)
        _raw_atomdb_getall(zvalsptr,avalsptr)
        za=_np.stack((zvals,avals)).T
        return za
    functions['atomdb_getall_za']=atomdb_getall_za

    _wrap('ncrystal_info_natominfo',_uint,(ncrystal_info_t,))
    _wrap('ncrystal_info_hasatommsd',_int,(ncrystal_info_t,))
    _raw_info_getatominfo = _wrap('ncrystal_info_getatominfo',None,(ncrystal_info_t,_uint,_uintp,_uintp,_dblp,_dblp),hide=True)
    def ncrystal_info_getatominfo(nfo,iatom):
        atomidx,n,dt,msd=_uint(),_uint(),_dbl(),_dbl()
        _raw_info_getatominfo(nfo,iatom,atomidx,n,dt,msd)
        return (atomidx.value,n.value,dt.value,msd.value)
    functions['ncrystal_info_getatominfo'] = ncrystal_info_getatominfo
    _raw_info_getatompos = _wrap('ncrystal_info_getatompos',None,(ncrystal_info_t,_uint,_uint,_dblp,_dblp,_dblp),hide=True)
    def ncrystal_info_getatompos(nfo,iatom,ipos):
        x,y,z=_dbl(),_dbl(),_dbl()
        _raw_info_getatompos(nfo,iatom,ipos,x,y,z)
        return x.value, y.value, z.value
    functions['ncrystal_info_getatompos'] = ncrystal_info_getatompos

    for s in ('temperature','xsectabsorption','xsectfree','density','numberdensity','sld'):
        _wrap('ncrystal_info_get%s'%s,_dbl,(ncrystal_info_t,))
    _wrap('ncrystal_info_getstateofmatter',_int,( ncrystal_info_t,))
    _raw_info_getstruct = _wrap('ncrystal_info_getstructure',_int,(ncrystal_info_t,_uintp,_dblp,_dblp,_dblp,_dblp,_dblp,_dblp,_dblp,_uintp))
    def ncrystal_info_getstructure(nfo):
        sg,natom=_uint(),_uint()
        a,b,c,alpha,beta,gamma,vol = _dbl(),_dbl(),_dbl(),_dbl(),_dbl(),_dbl(),_dbl(),
        if _raw_info_getstruct(nfo,sg,a,b,c,alpha,beta,gamma,vol,natom) == 0:
            return {}
        return dict(spacegroup=int(sg.value),a=a.value,b=b.value,c=c.value,alpha=alpha.value,
                    beta=beta.value,gamma=gamma.value,volume=vol.value,n_atoms=int(natom.value))
    functions['ncrystal_info_getstructure'] = ncrystal_info_getstructure

    _wrap('ncrystal_info_nphases',_int,(ncrystal_info_t,))
    _wrap('ncrystal_info_getphase',ncrystal_info_t,(ncrystal_info_t,_int,_dblp))

    _wrap('ncrystal_info_nhkl',_int,(ncrystal_info_t,))
    _wrap('ncrystal_info_hkl_dlower',_dbl,(ncrystal_info_t,))
    _wrap('ncrystal_info_hkl_dupper',_dbl,(ncrystal_info_t,))
    _wrap('ncrystal_info_braggthreshold',_dbl,(ncrystal_info_t,))
    _wrap('ncrystal_info_hklinfotype',_int,(ncrystal_info_t,))
    _wrap('ncrystal_info_gethkl',None,(ncrystal_info_t,_int,_intp,_intp,_intp,_intp,_dblp,_dblp))
    _wrap('ncrystal_info_dspacing_from_hkl',_dbl,(ncrystal_info_t,_int,_int,_int))
    functions['ncrystal_info_gethkl_setuppars'] = lambda : (_int(),_int(),_int(),_int(),_dbl(),_dbl())

    _raw_gethkl_allindices = _wrap('ncrystal_info_gethkl_allindices',None,(ncrystal_info_t,_int,_intp,_intp,_intp), hide=True )
    def iter_hkllist(nfo,all_indices=False):
        h,k,l,mult,dsp,fsq = _int(),_int(),_int(),_int(),_dbl(),_dbl()
        nhkl = int(functions['ncrystal_info_nhkl'](nfo))
        for idx in range(nhkl):
            _rawfct['ncrystal_info_gethkl'](nfo,idx,h,k,l,mult,dsp,fsq)
            if not all_indices:
                yield h.value,k.value,l.value,mult.value,dsp.value,fsq.value
            else:
                nc_assert( mult.value % 2 == 0 )
                n = mult.value // 2
                hvals, hvalsptr = _create_numpy_int_array( n )
                kvals, kvalsptr = _create_numpy_int_array( n )
                lvals, lvalsptr = _create_numpy_int_array( n )
                _raw_gethkl_allindices(nfo,idx,hvalsptr,kvalsptr,lvalsptr)
                yield hvals,kvals,lvals,mult.value,dsp.value,fsq.value
    functions['iter_hkllist']=iter_hkllist

    _wrap('ncrystal_info_ndyninfo',_uint,(ncrystal_info_t,))
    _raw_di_base = _wrap('ncrystal_dyninfo_base',None,(ncrystal_info_t,_uint,_dblp,_uintp,_dblp,_uintp),hide=True)
    _raw_di_scatknl = _wrap('ncrystal_dyninfo_extract_scatknl',None,(ncrystal_info_t,_uint,_uint,_dblp,_uintp,_uintp,_uintp,
                                                                     _dblpp,_dblpp,_dblpp,_dblpp),hide=True)
    _raw_di_vdos = _wrap('ncrystal_dyninfo_extract_vdos',None,(ncrystal_info_t,_uint,_dblp,_dblp,_uintp,_dblpp),hide=True)
    _raw_di_vdosdebye = _wrap('ncrystal_dyninfo_extract_vdosdebye',None,(ncrystal_info_t,_uint,_dblp),hide=True)
    _raw_di_vdos_input = _wrap('ncrystal_dyninfo_extract_vdos_input',None,(ncrystal_info_t,_uint,_uintp,_dblpp,_uintp,_dblpp),hide=True)

    def ncrystal_dyninfo_base(key):
        infoobj,dynidx = key
        fr,tt,atomindex,ditype=_dbl(),_dbl(),_uint(),_uint()
        _raw_di_base(infoobj,dynidx,fr,atomindex,tt,ditype)
        return (fr.value,tt.value,atomindex.value,ditype.value)
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
        _raw_di_vdos_input(infoobj,dynidx,negrid,ctypes.byref(egridptr),ndensity,ctypes.byref(densityptr))
        return (negrid.value,egridptr,ndensity.value,densityptr)
    functions['ncrystal_dyninfo_base'] = ncrystal_dyninfo_base
    functions['ncrystal_dyninfo_extract_scatknl'] = ncrystal_dyninfo_extract_scatknl
    functions['ncrystal_dyninfo_extract_vdos'] = ncrystal_dyninfo_extract_vdos
    functions['ncrystal_dyninfo_extract_vdosdebye'] = ncrystal_dyninfo_extract_vdosdebye
    functions['ncrystal_dyninfo_extract_vdos_input'] = ncrystal_dyninfo_extract_vdos_input

    _raw_vdoseval = _wrap('ncrystal_vdoseval',None,(_dbl,_dbl,_uint,_dblp,_dbl,_dbl,_dblp,_dblp,_dblp,_dblp,_dblp),hide=True)
    def nc_vdoseval(emin,emax,density,temp,mass_amu):
        msd,dt,g0,teff,oint=_dbl(),_dbl(),_dbl(),_dbl(),_dbl()
        _raw_vdoseval(emin,emax,len(density),ndarray_to_dblp(density),temp,mass_amu,
                      msd,dt,g0,teff,oint)
        return dict(msd=msd.value,debye_temp=dt.value,gamma0=g0.value,teff=teff.value,integral=oint.value)
    functions['nc_vdoseval']=nc_vdoseval

    _wrap('ncrystal_info_ncomponents',_uint,(ncrystal_info_t,))
    _raw_info_getcomp=_wrap('ncrystal_info_getcomponent',None,(ncrystal_info_t,_uint,_uintp,_dblp),hide=True)
    def ncrystal_info_getcomp(nfo,icomp):
        aidx,fraction=_uint(),_dbl()
        _raw_info_getcomp(nfo,icomp,aidx,fraction)
        return aidx.value,fraction.value
    functions['ncrystal_info_getcomp']=ncrystal_info_getcomp

    _NEP_FCTTYPE = ctypes.CFUNCTYPE( _uint,_uint,_uintp,_dblp )# -> unsigned (*natelemprovider)(unsigned,unsigned*,double*)
    _raw_info_getflatcompos = _wrap('ncrystal_get_flatcompos',_charptr,(ncrystal_info_t,_int,_NEP_FCTTYPE),hide=True)

    def nc_info_getflatcompos( nfo, natelemprovider = None, prefernatelem = True ):
        _prefnatelem = _int(1 if prefernatelem else 0)
        if not natelemprovider:
            _nullfct = ctypes.cast(None, _NEP_FCTTYPE)
            raw_str = _raw_info_getflatcompos(nfo,_prefnatelem,_nullfct)
        else:
            #natelemprovider is function which takes integer Z and returns a list
            #[(A,fraction),...]. Create a wrapper compatible with _NEP_FCTTYPE:
            #unsigned (*natelemprovider)(unsigned,unsigned*,double*)
            def nep_wrap(Z,bufA,bufFrac):
                _ = natelemprovider(Z)
                if not _:
                    #returned None or [], i.e. there is no info about the
                    #element:
                    return 0
                assert len(_) < 128
                for i,(A,f) in enumerate(_):
                    bufA[i] = A
                    bufFrac[i] = f
                return len(_)
            _c_nep_wrap = _NEP_FCTTYPE(nep_wrap)
            raw_str = _raw_info_getflatcompos(nfo,_prefnatelem,_c_nep_wrap)
        if raw_str is None:
            return 'null'#None in json
        res=_cstr2str(ctypes.cast(raw_str,_cstr).value)
        _raw_deallocstr(raw_str)
        return res
    functions['nc_info_getflatcompos']=nc_info_getflatcompos

    _wrap('ncrystal_create_atomdata',ncrystal_atomdata_t,(ncrystal_info_t,_uint))
    _raw_atomdata_subcomp = _wrap('ncrystal_create_atomdata_subcomp',ncrystal_atomdata_t,
                                  (ncrystal_atomdata_t,_uint,_dblp),hide=True)
    _raw_atomdata_getfields=_wrap('ncrystal_atomdata_getfields',None,(ncrystal_atomdata_t,_cstrp,_cstrp,
                                                                      _dblp,_dblp,_dblp,_dblp,
                                                                      _uintp,_uintp,_uintp),hide=True)
    _wrap('ncrystal_create_component_atomdata',ncrystal_atomdata_t,(ncrystal_info_t,_uint))

    def ncrystal_atomdata_createsubcomp(ad,icomp):
        fraction = _dbl()
        comp_ad = _raw_atomdata_subcomp(ad,icomp,fraction)
        return (comp_ad,fraction.value)
    functions['ncrystal_atomdata_createsubcomp']=ncrystal_atomdata_createsubcomp
    def ncrystal_atomdata_getfields(ad):
        mass_amu,sigma_inc,scatlen_coh,sigma_abs=_dbl(),_dbl(),_dbl(),_dbl()
        dl,descr=_cstr(),_cstr()
        ncomp,zval,aval = _uint(),_uint(),_uint()
        _raw_atomdata_getfields(ad,ctypes.byref(dl),ctypes.byref(descr),
                                mass_amu,sigma_inc,scatlen_coh,sigma_abs,
                                ncomp,zval,aval)
        return dict(m=mass_amu.value,incxs=sigma_inc.value,cohsl_fm=scatlen_coh.value,absxs=sigma_abs.value,
                    dl=_cstr2str(dl.value),descr=_cstr2str(descr.value),
                    ncomp=ncomp.value,z=zval.value,a=aval.value)
    functions['ncrystal_atomdata_getfields'] = ncrystal_atomdata_getfields

    _raw_ncustom = _wrap('ncrystal_info_ncustomsections',_uint,(ncrystal_info_t,),hide=True)
    _raw_csec_name = _wrap('ncrystal_info_customsec_name',_cstr,(ncrystal_info_t,_uint),hide=True)
    _raw_csec_nlines = _wrap('ncrystal_info_customsec_nlines',_uint,(ncrystal_info_t,_uint),hide=True)
    _raw_csec_nparts = _wrap('ncrystal_info_customline_nparts',_uint,(ncrystal_info_t,_uint,_uint),hide=True)
    _raw_csec_part = _wrap('ncrystal_info_customline_getpart',_cstr,(ncrystal_info_t,_uint,_uint,_uint),hide=True)
    def ncrystal_info_getcustomsections(nfo):
        n=_raw_ncustom(nfo)
        if n==0:
            return tuple()
        out=[]
        for isec in range(n):
            lines=[]
            secname = _cstr2str(_raw_csec_name(nfo,isec))
            nlines = _raw_csec_nlines(nfo,isec)
            for iline in range(nlines):
                nparts=_raw_csec_nparts(nfo,isec,iline)
                parts=[]
                for ipart in range(nparts):
                    parts.append(_cstr2str(_raw_csec_part(nfo,isec,iline,ipart)))
                lines.append(tuple(parts))
            out.append((secname,tuple(lines)))
        return tuple(out)
    functions['ncrystal_info_getcustomsections'] = ncrystal_info_getcustomsections

    _raw_reginmemfd = _wrap('ncrystal_register_in_mem_file_data',None,(_cstr,_cstr),hide=True)
    def ncrystal_register_in_mem_file_data(virtual_filename,data):
        _raw_reginmemfd(_str2cstr(virtual_filename),
                        _str2cstr(data))
    functions['ncrystal_register_in_mem_file_data']=ncrystal_register_in_mem_file_data

    def _prepare_many(ekin,repeat):
        if _np is None and not repeat is None:
            raise NCBadInput('Can not use "repeat" parameter when Numpy is absent on the system')
        if repeat is None and not hasattr(ekin,'__len__'):
            return None#scalar case, array interface not triggered
        repeat = 1 if repeat is None else repeat
        ekin = (ekin if hasattr(ekin,'ctypes') else _np.asfarray(ekin) ) if hasattr(ekin,'__len__') else _np.ones(1)*ekin
        #NB: returning the ekin object itself is important in order to keep a reference to it after the call:
        return ndarray_to_dblp(ekin),len(ekin),repeat,ekin

    _raw_xs_no = _wrap('ncrystal_crosssection_nonoriented',None,(ncrystal_process_t,_dbl,_dblp),hide=True)
    _raw_xs_no_many = _wrap('ncrystal_crosssection_nonoriented_many',None,(ncrystal_process_t,_dblp,_ulong,
                                                                           _ulong,_dblp),hide=True)
    def ncrystal_crosssection_nonoriented(scat,ekin,repeat=None):
        many = _prepare_many(ekin,repeat)
        if many is None:
            res = _dbl()
            _raw_xs_no(scat,ekin,res)
            return res.value
        else:
            ekin_ct,n_ekin,repeat,ekin_nparr = many
            xs, xs_ct = _create_numpy_double_array(n_ekin*repeat)
            _raw_xs_no_many(scat,ekin_ct,n_ekin,repeat,xs_ct)
            return xs
    functions['ncrystal_crosssection_nonoriented'] = ncrystal_crosssection_nonoriented

    _raw_domain = _wrap('ncrystal_domain',None,(ncrystal_process_t,_dblp,_dblp),hide=True)
    def ncrystal_domain(proc):
        a,b = _dbl(),_dbl()
        _raw_domain(proc,a,b)
        return (a.value,b.value)
    functions['ncrystal_domain'] = ncrystal_domain

    _raw_samplesct_iso =_wrap('ncrystal_samplescatterisotropic',None,(ncrystal_scatter_t,_dbl,_dblp,_dblp),hide=True)
    _raw_samplesct_iso_many =_wrap('ncrystal_samplescatterisotropic_many',None,
                                   (ncrystal_scatter_t,_dblp,_ulong,_ulong,_dblp,_dblp),hide=True)
    _raw_samplescat = _wrap('ncrystal_samplescatter',None,( ncrystal_scatter_t, _dbl,_dbl*3,_dblp,_dbl*3),hide=True)
    _raw_samplescat_many = _wrap('ncrystal_samplescatter_many',None,( ncrystal_scatter_t,_dbl,_dbl*3,_ulong,
                                                                      _dblp,_dblp,_dblp,_dblp),hide=True)
    def ncrystal_samplesct_iso(scat,ekin,repeat=None):
        many = _prepare_many(ekin,repeat)
        if many is None:
            ekin_final,mu = _dbl(),_dbl()
            _raw_samplesct_iso(scat,ekin,ekin_final,mu)
            return ekin_final.value,mu.value
        else:
            ekin_ct,n_ekin,repeat,ekin_nparr = many
            ekin_final, ekin_final_ct = _create_numpy_double_array(n_ekin*repeat)
            mu, mu_ct = _create_numpy_double_array(n_ekin*repeat)
            _raw_samplesct_iso_many(scat,ekin_ct,n_ekin,repeat,ekin_final_ct,mu_ct)
            return ekin_final,mu
    functions['ncrystal_samplesct_iso'] = ncrystal_samplesct_iso

    def ncrystal_samplesct(scat, ekin, direction, repeat):
        cdir = (_dbl * 3)(*direction)
        if not repeat:
            res_dir = (_dbl * 3)(0,0,0)
            res_ekin = _dbl()
            _raw_samplescat(scat,ekin,cdir,res_ekin,res_dir)
            return res_ekin.value,(res_dir[0],res_dir[1],res_dir[2])
        else:
            assert repeat>=1
            res_ekin, res_ekin_ct = _create_numpy_double_array(repeat)
            res_ux, res_ux_ct = _create_numpy_double_array(repeat)
            res_uy, res_uy_ct = _create_numpy_double_array(repeat)
            res_uz, res_uz_ct = _create_numpy_double_array(repeat)
            _raw_samplescat_many(scat,ekin,cdir,repeat,res_ekin_ct,res_ux_ct,res_uy_ct,res_uz_ct)
            return res_ekin,(res_ux,res_uy,res_uz)
    functions['ncrystal_samplesct']=ncrystal_samplesct

    _raw_xs = _wrap('ncrystal_crosssection',None,(ncrystal_process_t,_dbl,_dbl*3,_dblp),hide=True)
    def ncrystal_crosssection( proc, ekin, direction):
        res = _dbl()
        cdir = (_dbl * 3)(*direction)
        if hasattr(ekin,'__len__'):
            #Todo: this vectorises on the Python side which is slow!
            _ensure_numpy()
            def e2xs(e):
                _raw_xs(proc,e,cdir,res)
                return res.value
            return _np.vectorize(e2xs)(ekin)
        else:
            _raw_xs(proc,ekin,cdir,res)
            return res.value
    functions['ncrystal_crosssection'] = ncrystal_crosssection

    #Obsolete:
    _raw_gs_no = _wrap('ncrystal_genscatter_nonoriented',None,(ncrystal_scatter_t,_dbl,_dblp,_dblp),hide=True)
    _raw_gs_no_many = _wrap('ncrystal_genscatter_nonoriented_many',None,(ncrystal_scatter_t,_dblp,_ulong,
                                                                         _ulong,_dblp,_dblp),hide=True)
    def ncrystal_genscatter_nonoriented(scat,ekin,repeat=None):
        many = _prepare_many(ekin,repeat)
        if many is None:
            angle,de = _dbl(),_dbl()
            _raw_gs_no(scat,ekin,angle,de)
            return angle.value,de.value
        else:
            ekin_ct,n_ekin,repeat,ekin_nparr = many
            angle, angle_ct = _create_numpy_double_array(n_ekin*repeat)
            de, de_ct = _create_numpy_double_array(n_ekin*repeat)
            _raw_gs_no_many(scat,ekin_ct,n_ekin,repeat,angle_ct,de_ct)
            return angle,de
    functions['ncrystal_genscatter_nonoriented'] = ncrystal_genscatter_nonoriented
    _raw_gs = _wrap('ncrystal_genscatter',None,(ncrystal_scatter_t,_dbl,_dbl*3,_dbl*3,_dblp),hide=True)
    _raw_gs_many = _wrap('ncrystal_genscatter_many',None,(ncrystal_scatter_t,_dbl,_dbl*3,
                                                          _ulong,_dblp,_dblp,_dblp,_dblp),hide=True)
    def ncrystal_genscatter(scat, ekin, direction, repeat):
        cdir = (_dbl * 3)(*direction)
        if not repeat:
            res_dir = (_dbl * 3)(0,0,0)
            res_de = _dbl()
            _raw_gs(scat,ekin,cdir,res_dir,res_de)
            return (res_dir[0],res_dir[1],res_dir[2]),res_de.value
        else:
            assert repeat>=1
            res_ux, res_ux_ct = _create_numpy_double_array(repeat)
            res_uy, res_uy_ct = _create_numpy_double_array(repeat)
            res_uz, res_uz_ct = _create_numpy_double_array(repeat)
            res_de, res_de_ct = _create_numpy_double_array(repeat)
            _raw_gs_many(scat,ekin,cdir,repeat,res_ux_ct,res_uy_ct,res_uz_ct,res_de_ct)
            return (res_ux,res_uy,res_uz),res_de
    functions['ncrystal_genscatter']=ncrystal_genscatter

    _wrap('ncrystal_create_info',ncrystal_info_t,(_cstr,))
    _wrap('ncrystal_create_scatter',ncrystal_scatter_t,(_cstr,))
    _wrap('ncrystal_create_scatter_builtinrng',ncrystal_scatter_t,(_cstr,_ulong))
    _wrap('ncrystal_create_absorption',ncrystal_absorption_t,(_cstr,))

    _raw_multicreate_direct = _wrap('ncrystal_multicreate_direct',None,
                                    ( _cstr, _cstr, _cstr,
                                      ctypes.POINTER(ncrystal_info_t),
                                      ctypes.POINTER(ncrystal_scatter_t),
                                      ctypes.POINTER(ncrystal_absorption_t) ),hide=True)
    nullptr_ncrystal_info_t = ctypes.cast(None, ctypes.POINTER(ncrystal_info_t))
    nullptr_ncrystal_scatter_t = ctypes.cast(None, ctypes.POINTER(ncrystal_scatter_t))
    nullptr_ncrystal_absorption_t = ctypes.cast(None, ctypes.POINTER(ncrystal_absorption_t))

    def multicreate_direct(data,dataType,cfg_params,doI,doS,doA):
        rawi = ncrystal_info_t() if doI else None
        raws = ncrystal_scatter_t() if doS else None
        rawa = ncrystal_absorption_t() if doA else None
        _raw_multicreate_direct( _str2cstr(data),_str2cstr(dataType or "" ),_str2cstr(cfg_params or ""),
                                 ctypes.byref(rawi) if rawi else nullptr_ncrystal_info_t,
                                 ctypes.byref(raws) if raws else nullptr_ncrystal_scatter_t,
                                 ctypes.byref(rawa) if rawa else nullptr_ncrystal_absorption_t )
        return rawi,raws,rawa
    functions['multicreate_direct'] = multicreate_direct

    _wrap('ncrystal_setbuiltinrandgen',None,tuple())

    _RANDGENFCTTYPE = ctypes.CFUNCTYPE( _dbl )
    _raw_setrand    = _wrap('ncrystal_setrandgen',None,(_RANDGENFCTTYPE,),hide=True)
    def ncrystal_setrandgen(randfct):
        #Set random function, keeping references as needed (otherwise fct ptrs
        #kept on C++ side will suddenly stop working!) and casting None to a null-ptr.
        if not randfct:
            keepalive=(None,ctypes.cast(None, _RANDGENFCTTYPE))
        else:
            keepalive=(randfct,_RANDGENFCTTYPE(randfct))#keep refs!
        _keepalive.append(keepalive)
        _raw_setrand(keepalive[1])
    functions['ncrystal_setrandgen'] = ncrystal_setrandgen

    _wrap('ncrystal_clone_absorption',ncrystal_absorption_t,(ncrystal_absorption_t,))
    _wrap('ncrystal_clone_scatter',ncrystal_scatter_t,(ncrystal_scatter_t,))
    _wrap('ncrystal_clone_scatter_rngbyidx',ncrystal_scatter_t,(ncrystal_scatter_t,_ulong))
    _wrap('ncrystal_clone_scatter_rngforcurrentthread',ncrystal_scatter_t,(ncrystal_scatter_t,))
    _wrap('ncrystal_decodecfg_vdoslux',_uint,(_cstr,))
    _wrap('ncrystal_has_factory',_int,(_cstr,))
    _wrap('ncrystal_clear_caches',None,tuple())

    _wrap('ncrystal_rngsupportsstatemanip_ofscatter',_int,( ncrystal_scatter_t, ))
    _wrap('ncrystal_setrngstate_ofscatter',None,(ncrystal_scatter_t, _cstr))
    _raw_getrngstate_scat = _wrap('ncrystal_getrngstate_ofscatter',_charptr,( ncrystal_scatter_t,),hide=True)

    def nc_getrngstate_scat(rawscatobj):
        rawstate = _raw_getrngstate_scat(rawscatobj)
        if not rawstate:
            #null ptr, i.e. state manipulation is not supported
            return None
        state=_cstr2str(ctypes.cast(rawstate,_cstr).value)
        _raw_deallocstr(rawstate)
        return state
    functions['nc_getrngstate_scat']=nc_getrngstate_scat

    def _decode_and_dealloc_raw_str(raw_cstr):
        if not raw_cstr:
            return None
        res=_cstr2str(ctypes.cast(raw_cstr,_cstr).value)
        _raw_deallocstr(raw_cstr)
        return res

    _raw_dbg_process = _wrap('ncrystal_dbg_process',_charptr,(ncrystal_process_t,),hide=True)
    def nc_dbg_proc(rawprocobj):
        return json.loads( _decode_and_dealloc_raw_str( _raw_dbg_process( rawprocobj ) ) )
    functions['nc_dbg_proc']=nc_dbg_proc

    _raw_decodecfg_json = _wrap('ncrystal_decodecfg_json',_charptr,(_cstr,),hide=True)
    def nc_cfgstr2json(cfgstr):
        return _decode_and_dealloc_raw_str( _raw_decodecfg_json(_str2cstr(cfgstr) ) )
    functions['nc_cfgstr2json']=nc_cfgstr2json

    _raw_ncmat2json = _wrap('ncrystal_ncmat2json',_charptr,(_cstr,),hide=True)
    def nc_ncmat2json(tdname):
        return _decode_and_dealloc_raw_str( _raw_ncmat2json(_str2cstr(tdname) ) )
    functions['nc_ncmat2json']=nc_ncmat2json

    raw_proc_uid = _wrap('ncrystal_process_uid',_charptr,(ncrystal_process_t,),hide=True)
    functions['procuid'] = lambda rawproc : int(_decode_and_dealloc_raw_str(raw_proc_uid(rawproc)))
    raw_info_uid = _wrap('ncrystal_info_uid',_charptr,(ncrystal_info_t,),hide=True)
    functions['infouid'] = lambda rawinfo : int(_decode_and_dealloc_raw_str(raw_info_uid(rawinfo)))
    raw_info_underlyinguid = _wrap('ncrystal_info_underlyinguid',_charptr,(ncrystal_info_t,),hide=True)
    functions['infouid_underlying'] = lambda rawinfo : int(_decode_and_dealloc_raw_str(raw_info_underlyinguid(rawinfo)))

    _raw_normcfgstr = _wrap('ncrystal_normalisecfg',_charptr,(_cstr,),hide=True)
    def nc_normcfgstr(cfgstr):
        raw_str = _raw_normcfgstr(_str2cstr(cfgstr))
        if not raw_str:
            return None
        res=_cstr2str(ctypes.cast(raw_str,_cstr).value)
        _raw_deallocstr(raw_str)
        return res
    functions['nc_normcfgstr']=nc_normcfgstr

    _raw_gencfgdoc = _wrap('ncrystal_gencfgstr_doc',_charptr,(_int,),hide=True)
    def nc_gencfgdoc(mode):
        raw_str = _raw_gencfgdoc(mode)
        if not raw_str:
            return None
        res=_cstr2str(ctypes.cast(raw_str,_cstr).value)
        _raw_deallocstr(raw_str)
        return res
    functions['nc_gencfgdoc']=nc_gencfgdoc

    _raw_gettextdata = _wrap('ncrystal_get_text_data',_cstrp,(_cstr,),hide=True)
    _raw_deallocstr = _wrap('ncrystal_dealloc_string',None,(_charptr,),hide=True)

    def nc_gettextdata(name):
        l = _raw_gettextdata(_str2cstr(str(name)))
        assert l is not None
        n = 5
        res = [l[i].decode() for i in range(n)]
        assert isinstance(res[0],str)
        _raw_deallocstrlist(n,l)
        return res
    functions['nc_gettextdata'] = nc_gettextdata

    _raw_getfilelist = _wrap('ncrystal_get_file_list',None,(_uintp,_cstrpp),hide=True)
    _raw_deallocstrlist = _wrap('ncrystal_dealloc_stringlist',None,(_uint,_cstrp),hide=True)
    def ncrystal_get_filelist():
        n,l = _uint(),_cstrp()
        _raw_getfilelist(n,ctypes.byref(l))
        assert n.value%4==0
        res=[]
        for i in range(n.value//4):
            res += [ (l[i*4].decode(),l[i*4+1].decode(),l[i*4+2].decode(),l[i*4+3].decode()) ]
        _raw_deallocstrlist(n,l)
        return res
    functions['ncrystal_get_filelist'] = ncrystal_get_filelist

    _raw_getpluginlist = _wrap('ncrystal_get_plugin_list',None,(_uintp,_cstrpp),hide=True)
    def ncrystal_get_pluginlist():
        n,l = _uint(),_cstrp()
        _raw_getpluginlist(n,ctypes.byref(l))
        assert n.value%3==0
        res=[]
        for i in range(n.value//3):
            pluginname,filename,plugintype=l[i*3].decode(),l[i*3+1].decode(),l[i*3+2].decode()
            res+=[(pluginname,filename,plugintype)]
        _raw_deallocstrlist(n,l)
        return res
    functions['ncrystal_get_pluginlist'] = ncrystal_get_pluginlist

    _wrap('ncrystal_add_custom_search_dir',None,(_cstr,))
    _wrap('ncrystal_remove_custom_search_dirs',None,tuple())
    _wrap('ncrystal_enable_abspaths',None,(_int,))
    _wrap('ncrystal_enable_relpaths',None,(_int,))
    _wrap('ncrystal_enable_stddatalib',None,(_int,_cstr))
    _wrap('ncrystal_enable_stdsearchpath',None,(_int,))
    _wrap('ncrystal_remove_all_data_sources',None,tuple())
    return functions

_rawfct = _load(_find_nclib())

def decodecfg_packfact(cfgstr):
    """OBSOLETE FUNCTION (always returns 1.0 now)."""
    return 1.0

def decodecfg_vdoslux(cfgstr):
    """Extract vdoslux value from cfgstr"""
    return int(_rawfct['ncrystal_decodecfg_vdoslux'](_str2cstr(cfgstr)))

def createVDOSDebye(debye_temperature):
    """Create simplified VDOS according to the Debye model"""
    _ensure_numpy()
    #NB: Must keep function exactly synchronised with createVDOSDebye function
    #in .cc src (although leaving out temperature,boundXS,elementMassAMU args
    #here):
    debye_energy = constant_boltzmann*debye_temperature
    vdos_egrid = _np_linspace(0.5*debye_energy,debye_energy,20)
    scale = 1.0 / (debye_energy*debye_energy)
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
        self.__rawobj_byref = ctypes.byref(rawobj)#keep byref(rawobj), since ctypes might
                                                  #disappear before __del__ is called.
    def __del__(self):
        if hasattr(self,'_rawunref') and self._rawunref:
            self._rawunref(self.__rawobj_byref)
    def refCount(self):
        """Access reference count of wrapped C++ object"""
        return _rawfct['ncrystal_refcount'](self._rawobj)

def nc_assert(b,msg=""):
    """Assertion which throws NCLogicError on failure"""
    if not bool(b):
        raise NCLogicError(msg if msg else 'assertion failed')

class AtomData(RCBase):
    """Class providing physical constants related to a particular mix of
    isotopes. This can be used to represent elements (i.e. all isotopes having
    same Z) in either natural or enriched form, but can also be used to
    represent atoms in doped crystals. E.g. if a small fraction (0.1%) of
    Cr-ions replace some Al-ions in a Al2O3 lattice, the AtomData could
    represent a mix of 0.1% Cr and 99.9% Al.
    """
    def __init__(self,rawobj):
        """internal usage only"""
        super(AtomData, self).__init__(rawobj)
        f=_rawfct['ncrystal_atomdata_getfields'](rawobj)
        self.__m = f['m']
        self.__incxs = f['incxs']
        self.__cohsl_fm = f['cohsl_fm']
        self.__absxs = f['absxs']
        self.__dl = f['dl']
        self.__descr = f['descr']
        self.__ncomp = f['ncomp']
        self.__z = f['z']
        self.__a = f['a']
        self.__b2f = (self.__m/(self.__m+const_neutron_mass_amu))**2
        self.__comp = [None]*self.__ncomp
        self.__compalldone = (self.__ncomp==0)

    def averageMassAMU(self):
        """Atomic mass in Daltons (averaged appropriately over constituents)"""
        return self.__m
    def coherentScatLen(self):
        """Coherent scattering length in sqrt(barn)=10fm"""
        return self.__cohsl_fm*0.1#0.1 is fm/sqrt(barn)
    def coherentScatLenFM(self):
        """Coherent scattering length in fm"""
        return self.__cohsl_fm
    def coherentXS(self):
        """Bound coherent cross section in barn. Same as 4*pi*coherentScatLen()**2"""
        return _k4Pidiv100*self.__cohsl_fm**2
    def incoherentXS(self):
        """Bound incoherent cross section in barn"""
        return self.__incxs
    def scatteringXS(self):
        """Bound scattering cross section in barn (same as coherentXS()+incoherentXS())"""
        return self.__incxs+self.coherentXS()
    def captureXS(self):
        """Absorption cross section in barn"""
        return self.__absxs

    def freeScatteringXS(self):
        """Free scattering cross section in barn (same as freeCoherentXS()+freeIncoherentXS())"""
        return self.__b2f * self.scatteringXS()
    def freeCoherentXS(self):
        """Free coherent cross section in barn."""
        return self.__b2f * self.coherentXS()
    def freeIncoherentXS(self):
        """Free incoherent cross section in barn."""
        return self.__b2f * self.incoherentXS()

    def isNaturalElement(self):
        """Natural element with no composition."""
        return self.__z!=0 and self.__ncomp==0 and self.__a==0

    def isSingleIsotope(self):
        """Single isotope with no composition."""
        return self.__a!=0

    def isComposite(self):
        """Composite definition. See nComponents(), getComponent() and components property"""
        return self.__ncomp!=0

    def isElement(self):
        """If number of protons per nuclei is well defined. This is true for natural
           elements, single isotopes, and composites where all components
           have the same number of protons per nuclei."""
        return self.__z!=0

    def Z(self):
        """Number of protons per nuclei (0 if not well defined)."""
        return self.__z

    def elementName(self):
        """If Z()!=0, this returns the corresponding element name ('H', 'He', ...).
           Returns empty string when Z() is 0."""
        if not self.__z:
            return ''
        #NB: We are relying on natural elements to return their element names in
        #description(false). This is promised by a comment in NCAtomData.hh!
        if self.isNaturalElement():
            return self.__descr
        return atomDB(self.__z).description(False)

    def A(self):
        """Number of nucleons per nuclei (0 if not well defined or natural element)."""
        return self.__a

    class Component:
        def __init__(self,fr,ad):
            """internal usage only"""
            self.__fr = fr
            self.__ad = ad
            assert not ad.isTopLevel()
        @property
        def fraction(self):
            """Fraction (by count) of component in mixture"""
            return self.__fr
        @property
        def data(self):
            """AtomData of component"""
            return self.__ad
        def __str__(self):
            return '%g*AtomData(%s)'%(self.__fr,self.__ad.description(True))

    def nComponents(self):
        """Number of sub-components in a mixture"""
        return self.__ncomp
    def getComponent(self,icomponent):
        """Get component in a mixture"""
        c=self.__comp[icomponent]
        if c:
            return c
        rawobj_subc,fraction=_rawfct['ncrystal_atomdata_createsubcomp'](self._rawobj,icomponent)
        ad = AtomData(rawobj_subc)
        c = AtomData.Component(fraction,ad)
        self.__comp[icomponent] = c
        return c
    def getAllComponents(self):
        """Get list of all components"""
        if self.__compalldone:
            return self.__comp
        for i,c in enumerate(self.__comp):
            if not c:
                self.getComponent(i)
        self.__compalldone=True
        return self.__comp
    components = property(getAllComponents)

    def displayLabel(self):
        """Short label which unique identifies an atom role within a particular material."""
        return self.__dl

    def isTopLevel(self):
        """Whether or not AtomData appears directly on an Info object. If not,
         it will most likely either be a component (direct or indirect) of a top
         level AtomData object, or be taken from the composition list of a
         multi-phase object.
        """
        return bool(self.__dl)

    def description(self,includeValues=True):
        """Returns description of material as a string, with or without values."""
        if includeValues:
            zstr=' Z=%i'%self.__z if self.__z else ''
            astr=' A=%i'%self.__a if self.__a else ''
            _=(self.__descr,self.__cohsl_fm,self.coherentXS(),self.__incxs,
               self.__absxs,self.__m,zstr,astr)
            return'%s(cohSL=%gfm cohXS=%gbarn incXS=%gbarn absXS=%gbarn mass=%gamu%s%s)'%_
        return self.__descr

    def __str__(self):
        descr=self.description()
        return '%s=%s'%(self.__dl,descr) if self.__dl else descr

class StateOfMatter(enum.Enum):
    """State of matter. Note that Solid's might be either amorphous or crystalline."""
    #NB: List here must be synchronized with list and values in NCInfo.hh:
    Unknown = 0
    Solid = 1
    Gas = 2
    Liquid = 3

class HKLInfoType(enum.Enum):
    """Describes the kind of information about plane normals and Miller (hkl)
       indices available on each entry in the HKLList."""
    #NB: List here must be synchronized with list and values in NCInfoTypes.hh:
    SymEqvGroup = 0
    ExplicitHKLs = 1
    ExplicitNormals = 2
    Minimal = 3


class Info(RCBase):
    """Class representing information about a given material.

     Objects might represent either multi- or single phase
     materials. Multi-phase objects contain a list of phases (which might
     themselves be either single or multi-phase objects). Most other fields
     (structure, hkl lists, dynamics, etc.) are single-phase specific and will
     be unavailable on multiphase-phase objects. Exceptions are phase-structure
     information (todo) which is only available on multi-phase objects, and
     fields which are available for both multi- and single-phase objects such as
     density, composition, temperature, and state of matter (where such are well
     defined).

    """
    def __init__(self, cfgstr):
        """create Info object based on cfg-string (same as using createInfo(cfgstr))"""
        if isinstance(cfgstr,tuple) and len(cfgstr)==2 and cfgstr[0]=='_rawobj_':
            #Already got an ncrystal_info_t object:
            rawobj = cfgstr[1]
        else:
            rawobj = _rawfct['ncrystal_create_info'](_str2cstr(cfgstr))
        super(Info, self).__init__(rawobj)
        self.__dyninfo=None
        self.__atominfo=None
        self.__custom=None
        self.__atomdatas=[]
        self.__comp=None
        self._som=None
        self._nphases = int(_rawfct['ncrystal_info_nphases'](rawobj))
        assert self._nphases == 0 or self._nphases >= 2
        self._phases = tuple() if self._nphases == 0 else None

    def getUniqueID(self):
        """Unique identifier of object (UID)."""
        return _rawfct['infouid'](self._rawobj)
    uid = property(getUniqueID)

    def _getUnderlyingUniqueID(self):
        """Unique identifier of underlying object, which does not change on simple
        density or cfg-data overrides (expert usage only!)."""
        return _rawfct['infouid_underlying'](self._rawobj)

    def isSinglePhase(self):
        """Single phase object."""
        return self._nphases == 0

    def isMultiPhase(self):
        """Multi phase object."""
        return self._nphases != 0

    def __initPhases(self):
        assert self._phases is None and self._nphases > 1
        l=[]
        for i in range(self._nphases):
            fraction = ctypes.c_double()
            ph_info_raw = _rawfct['ncrystal_info_getphase'](self._rawobj,i,fraction)
            ph_info = Info( ('_rawobj_',ph_info_raw) )
            l.append( ( float(fraction.value), ph_info ) )
        self._phases = tuple(l)
        return self._phases

    def getPhases(self):
        """Daughter phases in a multi-phase object. Returns a list of fractions and Info
        objects of the daughter phases, in the format [(fraction_i,daughter_i),...]
        """
        return self.__initPhases() if self._phases is None else self._phases
    phases=property(getPhases)

    def _initComp(self):
        assert self.__comp is None
        nc = _rawfct['ncrystal_info_ncomponents'](self._rawobj)
        l = []
        for icomp in range(nc):
            atomidx,fraction = _rawfct['ncrystal_info_getcomp'](self._rawobj,icomp)
            #NB: atomidx will be invalid in case of multiphase objects!
            if atomidx < 65535:
                #Most likely a single-phase object with valid atomidx, we can
                #use self._provideAtomData and share the AtomData objects also here on the python side:
                l += [(fraction,self._provideAtomData(atomidx))]
            else:
                #Most likely a multi-phase object with invalid atomidx, we must
                #create new AtomData objects, based on ncrystal_create_component_atomdata:
                raw_ad = _rawfct['ncrystal_create_component_atomdata'](self._rawobj,icomp)
                obj = AtomData(raw_ad)
                assert not obj.isTopLevel()#does not appear directly on Info object
                l += [(fraction,obj)]
        self.__comp = l
        return self.__comp

    def stateOfMatter(self):
        """State of matter, i.e. Solid, Liquid, Gas, ... as per the options in the
        StateOfMatter class. Note that the .isCrystalline() method can be used
        to additionally distinguish between amorphous and crystalline
        solids. Return value is an enum object, whose .name() method can be used
        in case a string value is desired.
        """
        if self._som is None:
            self._som = StateOfMatter(_rawfct['ncrystal_info_getstateofmatter'](self._rawobj))
        return self._som

    def isCrystalline(self):
        """Whether or not object is crystalline (i.e. has unit cell structure or list of
        reflection planes)."""
        return self.hasStructureInfo() or self.hasAtomInfo() or self.hasHKLInfo()

    def hasComposition(self):
        """OBSOLETE FUNCTION (always available now)."""
        return True

    def getComposition(self):
        """Get basic composition as list of (fraction,AtomData). For a single-phase
        object, the list is always consistent with AtomInfo/DynInfo (if
        present).
        """
        return self._initComp() if self.__comp is None else self.__comp
    composition=property(getComposition)

    def getFlattenedComposition( self,
                                 preferNaturalElements = True,
                                 naturalAbundProvider = None,
                                 asJSONStr=False ):
        """Break down the basic composition of the material into elements and
        isotopes. If an element only occurs as a natural element and has no
        specific isotopes, that element will be returned as a natural isotope
        unless preferNaturalElements=False. Generally, it is best if a
        naturalAbundProvider is given, since if a given Z value has a mix of
        isotopes and the natural elements, the natural element must always be
        broken up.

        If a naturalAbundProvider is given, it must be a function taking a Z
        value and returning the breakdown into isotopes,
        [(A1,frac1),...,(An,fracn)]. It can return None or an empty list to
        indicate missing information.

        Returns a list of (Z,<breakdown>) tuples, where <breakdown> is again a
        list of tuples of A-values and associated abundances (A=0 indicates
        natural elements). If asJSONStr=true, the data structure will be
        returned as a JSON-encoded string, instead of a Python dictionary.
        """
        _js = _rawfct['nc_info_getflatcompos'](self._rawobj,naturalAbundProvider,preferNaturalElements)
        return _js if asJSONStr else json.loads(_js)

    def dump(self,verbose=0):
        """Dump contained information to standard output. Use verbose argument to set
        verbosity level to 0 (minimal), 1 (middle), 2 (most verbose)."""
        sys.stdout.flush()
        sys.stderr.flush()
        _rawfct['ncrystal_dump_verbose'](self._rawobj,min(999,max(0,int(verbose))))

    def hasTemperature(self):
        """Whether or not material has a temperature available"""
        return _rawfct['ncrystal_info_gettemperature'](self._rawobj)>-1

    def getTemperature(self):
        """Material temperature (in kelvin)"""
        t=_rawfct['ncrystal_info_gettemperature'](self._rawobj)
        nc_assert(t>-1)
        return t

    def hasGlobalDebyeTemperature(self):
        """OBSOLETE FUNCTION: The concept of global versus per-element Debye
           temperatures has been removed. Please iterate over AtomInfo objects
           instead (see getAtomInfos() function) and get the Debye Temperature
           from those. This function will be removed in a future release.
        """
        return False

    def getGlobalDebyeTemperature(self):
        """OBSOLETE FUNCTION: The concept of global versus per-element Debye
           temperatures has been removed. Please iterate over AtomInfo objects
           instead (see getAtomInfos() function) and get the Debye Temperature
           from those. Calling this function will always result in an exception
           thrown for now, and the function will be removed in a future release..
        """
        raise NCLogicError('The concept of global Debye temperatures has been removed. Iterate over'
                           +' AtomInfo objects instead and get the Debye temperature values from those.')
        return None

    def hasAtomDebyeTemp(self):
        """Whether AtomInfo objects are present and have Debye temperatures available
        (they will either all have them available, or none of them will have
        them available).
        """
        if self.__atominfo is None:
            self.__initAtomInfo()
        return self.__atominfo[3]

    def hasDebyeTemperature(self):
        """Alias for hasAtomDebyeTemp()."""
        return self.hasAtomDebyeTemp()

    def hasAnyDebyeTemperature(self):
        """OBSOLETE FUNCTION which will be removed in a future release. Please
           call hasDebyeTemperature() instead.
        """
        return self.hasAtomDebyeTemp()

    def getDebyeTemperatureByElement(self,atomdata):
        """OBSOLETE FUNCTION which will be removed in a future release. Please access
           the AtomInfo objects instead and query the Debye temperature there.
        """
        if atomdata.isTopLevel():
            for ai in self.atominfos:
                if atomdata is ai.atomData:
                    return ai.debyeTemperature
        raise NCBadInput('Invalid atomdata object passed to Info.getDebyeTemperatureByElement'
                         +' (must be top-level AtomData from the same Info object)')

    def hasDensity(self):
        """OBSOLETE FUNCTION (densities are now always available)."""
        return True
    def getDensity(self):
        """Get density in g/cm^3. See also getNumberDensity()."""
        t=_rawfct['ncrystal_info_getdensity'](self._rawobj)
        nc_assert(t>0.0)
        return t
    density = property(getDensity)

    def hasNumberDensity(self):
        """OBSOLETE FUNCTION (densities are now always available)."""
        return True
    def getNumberDensity(self):
        """Get number density in atoms/angstrom^3. See also getDensity()."""
        t=_rawfct['ncrystal_info_getnumberdensity'](self._rawobj)
        nc_assert(t>0.0)
        return t
    numberdensity = property(getNumberDensity)

    def hasXSectAbsorption(self):
        """OBSOLETE FUNCTION"""
        return True
    def getXSectAbsorption(self):
        """Absorption cross section in barn (at 2200m/s)"""
        t=_rawfct['ncrystal_info_getxsectabsorption'](self._rawobj)
        nc_assert(t>-1)
        return t

    def hasXSectFree(self):
        """OBSOLETE FUNCTION"""
        return True
    def getXSectFree(self):
        """Saturated (free) scattering cross section in barn in the high-E limit"""
        t=_rawfct['ncrystal_info_getxsectfree'](self._rawobj)
        nc_assert(t>-1)
        return t

    def getSLD(self):
        """Get scattering length density in 1e-6/Aa^2"""
        return _rawfct['ncrystal_info_getsld'](self._rawobj)
    sld = property(getSLD)

    def hasStructureInfo(self):
        """Whether or not material has crystal structure information available."""
        return bool(_rawfct['ncrystal_info_getstructure'](self._rawobj))
    def getStructureInfo(self):
        """Information about crystal structure."""
        d=_rawfct['ncrystal_info_getstructure'](self._rawobj)
        nc_assert(d)
        return d

    def _provideAtomData(self,atomindex):
        if atomindex >= len(self.__atomdatas):
            if atomindex >= 65535:
                raise NCLogicError(f'Invalid atomindex ({atomindex}) provided to Info._provideAtomData')
            self.__atomdatas.extend([None,]*(atomindex+1-len(self.__atomdatas)))
        obj = self.__atomdatas[atomindex]
        if obj:
            return obj
        raw_ad = _rawfct['ncrystal_create_atomdata'](self._rawobj,atomindex)
        obj = AtomData(raw_ad)
        assert obj.isTopLevel()
        self.__atomdatas[atomindex] = obj
        return obj

    class AtomInfo:
        """Class with information about a particular atom in a unit cell, including the
        composition of atoms, positions, Debye temperature, and mean-squared-displacements.
        """

        def __init__(self,theinfoobj,atomidx,n,dt,msd,pos):
            """For internal usage only."""
            assert dt is None or ( isinstance(dt,float) and dt > 0.0 )
            assert msd is None or ( isinstance(msd,float) and msd > 0.0 )
            self._info_wr = weakref.ref(theinfoobj)
            self._atomidx,self.__n,self.__dt,self.__msd,=atomidx,n,dt,msd
            self.__pos = tuple(pos)#tuple, since it is immutable
            self.__atomdata = None
            self.__correspDI_wp = None

        def correspondingDynamicInfo(self):
            """
            Get corresponding DynamicInfo object from the same Info
            object. Returns None if Info object does not have dynamic info
            available
            """
            if self.__correspDI_wp is not None:
                if self.__correspDI_wp == False:
                    return None
                di = self.__correspDI_wp()
                nc_assert(di is not None,"AtomInfo.correspondingDynamicInfo can not"
                          +" be used after associated Info object is deleted")
                return di
            _info = self._info_wr()
            nc_assert(_info is not None,"AtomInfo.correspondingDynamicInfo can not"
                      +" be used after associated Info object is deleted")
            if not _info.hasDynamicInfo():
                self.__correspDI_wp = False
                return None
            for di in _info.dyninfos:
                if di._atomidx == self._atomidx:
                    self.__correspDI_wp = weakref.ref(di)
                    return di
            nc_assert(False,"AtomInfo.correspondingDynamicInfo: inconsistent internal state (bug?)")
        dyninfo = property(correspondingDynamicInfo)

        @property
        def atomData(self):
            """Return AtomData object with details about composition and relevant physics constants"""
            if self.__atomdata is None:
                _info = self._info_wr()
                nc_assert(_info is not None,"AtomInfo.atomData can not be used after associated Info object is deleted")
                self.__atomdata = _info._provideAtomData(self._atomidx)
                assert self.__atomdata.isTopLevel()
            return self.__atomdata

        @property
        def count(self):
            """Number of atoms of this type per unit cell"""
            return self.__n

        @property
        def debyeTemperature(self):
            """The Debye Temperature of the atom (kelvin). Returns None if not available."""
            return self.__dt

        @property
        def meanSquaredDisplacement(self):
            """The mean-squared-displacement of the atom (angstrom^2). Returns None if not
               available.
            """
            return self.__msd
        msd=meanSquaredDisplacement#alias

        @property
        def positions(self):
            """List (tuple actually) of positions of this atom in the unit cell. Each
            entry is given as a tuple of three values, (x,y,z)"""
            return self.__pos

        @property
        def atomIndex(self):
            """Index of atom on this material"""
            return self._atomidx

        def __str__(self):
            l=[str(self.atomData.displayLabel()),str(self.__n)]
            if self.__dt>0.0:
                l.append('DebyeT=%gK'%self.__dt if self.__dt else 'DebyeT=n/a')
            if self.__msd>0.0:
                l.append('MSD=%gAa^2'%self.__msd if self.__msd else 'MSD=n/a')
            l.append('hasPositions=%s'%('yes' if self.__pos else 'no'))
            return 'AtomInfo(%s)'%(', '.join(l))

    def hasAtomInfo(self):
        """Whether or no getAtomInfo()/atominfos are available"""
        if self.__atominfo is None:
            self.__initAtomInfo()
        return self.__atominfo[0]

    def hasAtomMSD(self):
        """Whether AtomInfo objects have mean-square-displacements available"""
        if self.__atominfo is None:
            self.__initAtomInfo()
        return self.__atominfo[1]

    def hasAtomPositions(self):
        """OBSOLETE FUNCTION: AtomInfo objects now always have positions
           available. Returns same as hasAtomInfo(). Will be removed in a future
           release.
        """
        return self.hasAtomInfo()

    def hasPerElementDebyeTemperature(self):
        """OBSOLETE FUNCTION which will be removed in a future
           release. Please use hasAtomDebyeTemp() instead.
        """
        return self.hasAtomDebyeTemp()

    def getAtomInfo(self):
        """Get list of AtomInfo objects, one for each atom. Returns empty list if unavailable."""
        if self.__atominfo is None:
            self.__initAtomInfo()
        return self.__atominfo[2]
    atominfos = property(getAtomInfo)

    def __initAtomInfo(self):
        assert self.__atominfo is None
        natoms = _rawfct['ncrystal_info_natominfo'](self._rawobj)
        hasmsd = bool(_rawfct['ncrystal_info_hasatommsd'](self._rawobj))
        hasperelemdt=False
        l=[]
        for iatom in range(natoms):
            atomidx,n,dt,msd = _rawfct['ncrystal_info_getatominfo'](self._rawobj,iatom)
            if dt:
                hasperelemdt=True
            assert hasmsd == (msd>0.0)
            pos=[]
            for ipos in range(n):
                pos.append( _rawfct['ncrystal_info_getatompos'](self._rawobj,iatom,ipos) )
            l.append( Info.AtomInfo(self,atomidx, n,
                                    ( dt if ( dt and  dt>0.0) else None),
                                    (msd if (msd and msd>0.0) else None),
                                    pos) )
        self.__atominfo = ( natoms>0, hasmsd, l, hasperelemdt )

    def hasHKLInfo(self):
        """Whether or not material has lists of HKL-plane info available"""
        return bool(_rawfct['ncrystal_info_nhkl'](self._rawobj)>-1)
    def nHKL(self):
        """Number of HKL planes available (grouped into families with similar d-spacing and f-squared)"""
        return int(_rawfct['ncrystal_info_nhkl'](self._rawobj))
    def hklDLower(self):
        """Lower d-spacing cutoff (angstrom)."""
        return float(_rawfct['ncrystal_info_hkl_dlower'](self._rawobj))
    def hklDUpper(self):
        """Upper d-spacing cutoff (angstrom)."""
        return float(_rawfct['ncrystal_info_hkl_dupper'](self._rawobj))
    def hklList(self,all_indices=False):
        """Iterator over HKL info, yielding tuples in the format
        (h,k,l,multiplicity,dspacing,fsquared). Running with all_indices=True to
        get the full list of hkl points in each group - in that case, h, k, and
        l will be numpy arrays of length multiplicity/2 (including just one of
        (h,k,l) and (-h,-k,-l) in the list).
        """
        nc_assert(self.hasHKLInfo())
        return _rawfct['iter_hkllist']( self._rawobj,
                                        all_indices = all_indices )
    def getBraggThreshold(self):
        """Get Bragg threshold in Aa (returns None if non-crystalline). This
        method is meant as a fast way to access the Bragg threshold without
        necessarily triggering a full initialisation of all HKL planes.
        """
        bt = float(_rawfct['ncrystal_info_braggthreshold'](self._rawobj))
        return bt if bt > 0.0 else None
    braggthreshold = property(getBraggThreshold)

    def hklInfoType(self):
        """What kind of information about plane normals and Miller indices are
        available in the hklList(). It is guaranteed to be the same for all
        HKLInfo entries, and will return "Minimal" when hklList() is present but
        empty. Like getBraggThreshold(), calling this method will not
        necessarily trigger a full initialisation of the hklList()."""
        return HKLInfoType(int(_rawfct['ncrystal_info_hklinfotype'](self._rawobj)))

    def dspacingFromHKL(self, h, k, l):
        """Convenience method, calculating the d-spacing of a given Miller
        index. Calling this incurs the overhead of creating a reciprocal lattice
        matrix from the structure info."""
        return float(_rawfct['ncrystal_info_dspacing_from_hkl'](self._rawobj,h,k,l))

    class DynamicInfo:
        """Class representing dynamic information (related to inelastic scattering)
           about a given atom"""

        def __init__(self,theinfoobj,fr,atomidx,tt,key):
            """internal usage only"""
            self._info_wr,self.__atomdata = weakref.ref(theinfoobj), None
            self.__fraction, self._atomidx, self._key, self.__tt = fr,atomidx,key,tt
            self.__correspAtomInfo_wp = None

        def correspondingAtomInfo(self):
            """Get corresponding AtomInfo object from the same Info object. Returns None if Info object does not have AtomInfo available"""
            if self.__correspAtomInfo_wp is not None:
                if self.__correspAtomInfo_wp == False:
                    return None
                ai = self.__correspAtomInfo_wp()
                nc_assert(ai is not None,"DynamicInfo.correspondingAtomInfo can not be used after associated Info object is deleted")
                return ai
            _info = self._info_wr()
            nc_assert(_info is not None,"DynamicInfo.correspondingAtomInfo can not be used after associated Info object is deleted")
            if not _info.hasAtomInfo():
                self.__correspAtomInfo_wp = False
                return None
            for ai in _info.atominfos:
                if ai._atomidx == self._atomidx:
                    self.__correspAtomInfo_wp = weakref.ref(ai)
                    return ai
            nc_assert(False,"DynamicInfo.correspondingAtomInfo: inconsistent internal state (bug?)")
        atominfo = property(correspondingAtomInfo)

        @property
        def atomIndex(self):
            """Index of atom on this material"""
            return self._atomidx

        @property
        def fraction(self):
            """Atom fraction in material (all fractions must add up to unity)"""
            return self.__fraction
        @property
        def temperature(self):
            """Material temperature (same value as on associated Info object)"""
            return self.__tt
        @property
        def atomData(self):
            """Return AtomData object with details about composition and relevant physics constants"""
            if self.__atomdata is None:
                _info = self._info_wr()
                nc_assert(_info is not None,"DynamicInfo.atomData can not be used after associated Info object is deleted")
                self.__atomdata = _info._provideAtomData(self._atomidx)
                assert self.__atomdata.isTopLevel()
            return self.__atomdata
        def _np(self):
            _ensure_numpy()
            return _np
        def _copy_cptr_2_nparray(self,cptr,n):
            np = self._np()
            return np.copy(np.ctypeslib.as_array(cptr, shape=(n,)))

        def __str__(self):
            n=self.__class__.__name__
            if n.startswith('DI_'):
                n=n[3:]
            s=', %s'%self._extradescr() if hasattr(self,'_extradescr') else ''
            return ('DynamicInfo(%s, fraction=%.4g%%, type=%s%s)'%(self.atomData.displayLabel(),
                                                                 self.__fraction*100.0,
                                                                 n,s))
    class DI_Sterile(DynamicInfo):
        """Class indicating atoms for which inelastic neutron scattering is absent
           or disabled."""
        pass

    class DI_FreeGas(DynamicInfo):
        """Class indicating atoms for which inelastic neutron scattering should be
           modelled as scattering on a free gas."""
        pass

    class DI_ScatKnl(DynamicInfo):
        """Base class indicating atoms for which inelastic neutron scattering will
           be, directly or indirectly, described by a scattering kernel,
           S(alpha,beta). This is an abstract class, and derived classes provide
           actual access to the kernels.
        """

        def __init__(self,theinfoobj,fr,atomidx,tt,key):
            """internal usage only"""
            super(Info.DI_ScatKnl, self).__init__(theinfoobj,fr,atomidx,tt,key)
            self.__lastknl,self.__lastvdoslux = None,None

        def _loadKernel( self, vdoslux = 3 ):
            assert isinstance(vdoslux,numbers.Integral) and 0<=vdoslux<=5
            vdoslux=int(vdoslux)
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
        def __init__(self,theinfoobj,fr,atomidx,tt,key):
            """internal usage only"""
            super(Info.DI_ScatKnlDirect, self).__init__(theinfoobj,fr,atomidx,tt,key)
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
        def __init__(self,theinfoobj,fr,atomidx,tt,key):
            """internal usage only"""
            super(Info.DI_VDOS, self).__init__(theinfoobj,fr,atomidx,tt,key)
            self.__vdosdata = None
            self.__vdosegrid_expanded = None
            self.__vdosorig = None

        def _extradescr(self):
            return 'npts=%i'%len(self.vdosOrigDensity())

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
                self.__vdosegrid_expanded = _np_linspace(_[0][0],_[0][1],len(_[1]))
            return self.__vdosegrid_expanded

        @property
        def vdos_density(self):
            """Access the VDOS density array"""
            return self.vdosData()[1]

        def loadKernel( self, vdoslux = 3 ):
            """Converts VDOS to S(alpha,beta) kernel with a luxury level given by the vdoslux parameter."""
            return self._loadKernel(vdoslux=vdoslux)

        def analyseVDOS(self):
            """Same as running the global analyseVDOS function on the contained VDOS."""
            return analyseVDOS(*self.vdos_egrid,self.vdos_density,
                               self.temperature,self.atomData.averageMassAMU())

    class DI_VDOSDebye(DI_ScatKnl):
        """Similarly to DI_VDOS, but instead of using a phonon VDOS spectrum provided
           externally, an idealised spectrum is used for lack of better
           options. This spectrum is based on the Debye Model, in which the
           spectrum rises quadratically with phonon energy below a cutoff value,
           kT, where T is the Debye temperature
        """

        def __init__(self,theinfoobj,fr,atomidx,tt,key):
            """internal usage only"""
            super(Info.DI_VDOSDebye, self).__init__(theinfoobj,fr,atomidx,tt,key)
            self.__vdosdata = None
            self.__debyetemp = None
            self.__vdosegrid_expanded = None

        def vdosData(self):
            """Access the idealised VDOS as ([egrid_min,egrid_max],vdos_density)"""
            if self.__vdosdata is None:
                self.__vdosdata = createVDOSDebye(self.debyeTemperature())
            return self.__vdosdata

        def debyeTemperature(self):
            """The Debye temperature of the atom"""
            if self.__debyetemp is None:
                self.__debyetemp = _rawfct['ncrystal_dyninfo_extract_vdosdebye'](self._key)
            return self.__debyetemp

        def _extradescr(self):
            return 'TDebye=%gK'%self.debyeTemperature()

        @property
        def vdos_egrid(self):
            """Access the VDOS energy grid as [egrid_min,egrid_max]"""
            return self.vdosData()[0]

        @property
        def vdos_egrid_expanded(self):
            """Access the egrid expanded into all actual egrid points"""
            if self.__vdosegrid_expanded is None:
                _ = self.vdosData()
                self.__vdosegrid_expanded = _np_linspace(_[0][0],_[0][1],len(_[1]))
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
        """Whether or not dynamic information for each atom is present"""
        return int(_rawfct['ncrystal_info_ndyninfo'](self._rawobj))>0 if self.__dyninfo is None else bool(self.__dyninfo)

    def getDynamicInfoList(self):
        """Get list of DynamicInfo objects (if available). One for each atom."""
        if self.__dyninfo is None:
            self.__dyninfo = []
            for idx in range(int(_rawfct['ncrystal_info_ndyninfo'](self._rawobj))):
                key = (self._rawobj,idx)
                fr,tt,atomidx,ditype = _rawfct['ncrystal_dyninfo_base'](key)
                args=(self,fr,atomidx,tt,key)
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
    dyninfos = property(getDynamicInfoList)

    def getAllCustomSections(self):
        """Custom information for which the core NCrystal code does not have any
        specific treatment. This is usually intended for usage by developers adding new
        experimental physics models."""

        if self.__custom is None:
            self.__custom = _rawfct['ncrystal_info_getcustomsections'](self._rawobj)
        return self.__custom
    customsections = property(getAllCustomSections)

class Process(RCBase):
    """Base class for calculations of processes in materials.

    Note that kinetic energies are in electronvolt and direction vectors are
    tuples of 3 numbers.

    """
    def getCalcName(self):
        """Obsolete alias for getName"""
        return self.getName()
    def getName(self):
        """Process name"""
        return _cstr2str(_rawfct['ncrystal_name'](self._rawobj))
    name = property(getName)

    def getUniqueID(self):
        """UID of underlying ProcImpl::Process object."""
        return _rawfct['procuid'](self._rawobj)
    uid = property(getUniqueID)

    def domain(self):
        """Domain where process has non-vanishing cross section.

        Returns the domain as (ekin_low,ekin_high). Outside this range of
        neutron kinetic energy, the process can be assumed to have vanishing
        cross sections. Thus, processes present at all energies will return
        (0.0,infinity).

        """
        return _rawfct['ncrystal_domain'](self._rawobj)

    def isNull(self):
        """Domain might indicate that this is a null-process, vanishing everywhere."""
        elow,ehigh = self.domain()
        #checking for inf like the following to avoid depending on numpy or math
        #modules just for this:
        return ( elow >= ehigh or ( elow>1e99 and elow==float('inf') ) )

    def isNonOriented(self):
        """opposite of isOriented()"""
        return bool(_rawfct['ncrystal_isnonoriented'](self._rawobj))
    def isOriented(self):
        """Check if process is oriented and results depend on the incident direction of the neutron"""
        return not self.isNonOriented()
    def crossSection( self, ekin, direction ):
        """Access cross sections."""
        return _rawfct['ncrystal_crosssection'](self._rawobj,ekin, direction)
    def crossSectionIsotropic( self, ekin, repeat = None ):
        """Access cross sections (should not be called for oriented processes).

        For efficiency it is possible to provide the ekin parameter as a numpy
        array of numbers and get a corresponding array of cross sections
        back. Likewise, the repeat parameter can be set to a positive number,
        causing the ekin value(s) to be reused that many times and a numpy array
        with results returned.

        """
        return _rawfct['ncrystal_crosssection_nonoriented'](self._rawobj,ekin,repeat)

    #Backwards compatible alias:
    crossSectionNonOriented = crossSectionIsotropic

    def xsect(self,ekin=None,direction=None,wl=None,repeat=None):
        """Convenience function which redirects calls to either crossSectionIsotropic
        or crossSection depending on whether or not a direction is given. It can
        also accept wavelengths instead of kinetic energies via the wl
        parameter. The repeat parameter is currently only supported when
        direction is not provided.
        """
        ekin = Process._parseekin( ekin, wl )
        if direction is None:
            return self.crossSectionIsotropic( ekin, repeat )
        else:
            if repeat is None:
                return self.crossSection( ekin, direction )
            else:
                raise NCBadInput('The repeat parameter is not currently supported when the direction parameter is also provided.')

    @staticmethod
    def _parseekin(ekin,wl):
        if wl is None:
            if ekin is None:
                raise NCBadInput('Please provide either one of the "ekin" or "wl" parameters.')
            return ekin
        else:
            if ekin is not None:
                raise NCBadInput('Do not provide both "ekin" and "wl" parameters')
            return wl2ekin(wl)

    def getSummary(self,short = False ):
        """By default access a high-level summary of the process in the form of
           a dictionary holding various information which is also available on
           the underlying C++ process object. If instead short==True, what is
           instead returned is simply a short process label along with a
           (recursive) list of sub-components and their scales (if
           appropriate). Finally, if short=='printable', the returned object
           will instead be a list of strings suitable for a quick printout (each
           string is one line of printout).

        """
        #Not caching, method is likely to be called sparringly.

        #short printable:
        if short == 'printable':
            toplbl,comps=self.getSummary(short=True)
            l=[ toplbl]
            def add_lines( comps, indentlvl = 1 ):
                ncomps = len(comps)
                prefix = '   '*indentlvl
                for i,(scale,(lbl,subcomps)) in enumerate(comps):
                    smb = '\--' if i+1==ncomps else '|--'
                    scale_str = '' if scale==1.0 else f'{scale:g} * '
                    l.append(f'{prefix}{smb} {scale_str}{lbl}')
                    if subcomps:
                        add_lines( subcomps, indentlvl + 1 )
            if comps:
                add_lines( comps )
            return l

        d=_rawfct['nc_dbg_proc'](self._rawobj)
        #full:
        if not short:
            return d
        #short:
        def fmt_lbl(proc):
            name = proc['name']
            summarystr = proc['specific'].get('summarystr','')
            return f'{name}({summarystr})' if summarystr else name
        def extract_subcomponents(proc):
            subprocs = proc.get('specific',{}).get('components',[])
            return (fmt_lbl(proc),list( (scl, extract_subcomponents(sp)) for scl,sp in subprocs ))
        return extract_subcomponents(d)

    def dump(self,prefix=''):
        """Print a quick high level summary of the process to stdout. What is printed is
        in fact the lines resulting from a call to self.getSummary(short='printable'),
        with an optional prefix prepended to each line.
        """
        print(prefix+f'\n{prefix}'.join(self.getSummary(short='printable')))

class Absorption(Process):
    """Base class for calculations of absorption in materials"""

    def __init__(self, cfgstr):
        """create Absorption object based on cfg-string (same as using createAbsorption(cfgstr))"""
        if isinstance(cfgstr,tuple) and len(cfgstr)==2 and cfgstr[0]=='_rawobj_':
            #Cloning:
            rawobj_abs = cfgstr[1]
        else:
            rawobj_abs = _rawfct['ncrystal_create_absorption'](_str2cstr(cfgstr))
        self._rawobj_abs = rawobj_abs
        rawobj_proc = _rawfct['ncrystal_cast_abs2proc'](rawobj_abs)
        super(Absorption, self).__init__(rawobj_proc)

    def clone(self):
        """Clone object. The clone will be using the same physics models and sharing any
         read-only data with the original, but will be using its own private copy of any
         mutable caches. All in all, this means that the objects are safe to use
         concurrently in multi-threaded programming, as long as each thread gets
         its own clone. Return value is the new Absorption object.
        """
        newrawobj = _rawfct['ncrystal_clone_absorption'](self._rawobj_abs)
        return Absorption( ('_rawobj_',newrawobj) )

class Scatter(Process):

    """Base class for calculations of scattering in materials.

    Note that kinetic energies are in electronvolt and direction vectors are
    tuples of 3 numbers.

    """

    def __init__(self, cfgstr):
        """create Scatter object based on cfg-string (same as using createScatter(cfgstr))"""
        if isinstance(cfgstr,tuple) and len(cfgstr)==2 and cfgstr[0]=='_rawobj_':
            #Already got an ncrystal_scatter_t object:
            self._rawobj_scat = cfgstr[1]
        else:
            self._rawobj_scat = _rawfct['ncrystal_create_scatter'](_str2cstr(cfgstr))
        rawobj_proc = _rawfct['ncrystal_cast_scat2proc'](self._rawobj_scat)
        super(Scatter, self).__init__(rawobj_proc)


    def clone(self,rng_stream_index=None,for_current_thread=False):
        """Clone object. The clone will be using the same physics models and sharing any
        read-only data with the original, but will be using its own private copy
        of any mutable caches and will get an independent RNG stream. All in
        all, this means that the objects are safe to use concurrently in
        multi-threaded programming, as long as each thread gets its own
        clone. Return value is the new Scatter object.

        If greater control over RNG streams are needed, it is optionally allowed
        to either set rng_stream_index to a non-negative integral value, or set
        for_current_thread=True.

        If rng_stream_index is set, the resulting object will use a specific
        rngstream index. All objects with the same indeed will share the same
        RNG state, so a sensible strategy is to use the same index for all
        scatter objects which are to be used in the same thread:

        If setting for_current_thread=True, the resulting object will use a
        specific rngstream which has been set aside for the current thread. Thus
        this function can be called from a given work-thread, in order to get
        thread-safe scatter handle, with all objects cloned within the same
        thread sharing RNG state.

        """
        if rng_stream_index is not None:
            if for_current_thread:
                raise NCBadInput('Scatter.clone(..): do not set both rng_stream_index and for_current_thread parameters')
            if not isinstance(rng_stream_index, numbers.Integral) or not 0 <= rng_stream_index <= 4294967295:
                raise NCBadInput('Scatter.clone(..): rng_stream_index must be integral and in range [0,4294967295]')
            newrawobj = _rawfct['ncrystal_clone_scatter_rngbyidx'](self._rawobj_scat,int(rng_stream_index))
        elif for_current_thread:
            newrawobj = _rawfct['ncrystal_clone_scatter_rngforcurrentthread'](self._rawobj_scat)
        else:
            newrawobj = _rawfct['ncrystal_clone_scatter'](self._rawobj_scat)
        return Scatter( ('_rawobj_',newrawobj) )

    def sampleScatter( self, ekin, direction, repeat = None ):
        """Randomly generate scatterings.

        Assuming a scattering took place, generate final state of neutron based
        on current kinetic energy and direction. Returns
        tuple(ekin_final,direction_final) where direct_final is itself a tuple
        (ux,uy,uz). The repeat parameter can be set to a positive number,
        causing the scattering to be sampled that many times and numpy arrays
        with results returned.

        """
        return _rawfct['ncrystal_samplesct'](self._rawobj_scat,ekin,direction,repeat)


    def sampleScatterIsotropic( self, ekin, repeat = None ):
        """Randomly generate scatterings (should not be called for oriented processes).

        Assuming a scattering took place, generate final state of
        neutron. Returns tuple(ekin_final,mu) where mu is the cosine of the
        scattering angle. For efficiency it is possible to provide the ekin
        parameter as a numpy array of numbers and get corresponding arrays of
        energies and mu back. Likewise, the repeat parameter can be
        set to a positive number, causing the ekin value(s) to be reused that
        many times and numpy arrays with results returned.

        """
        return _rawfct['ncrystal_samplesct_iso'](self._rawobj_scat,ekin,repeat)

    def generateScattering( self, ekin, direction, repeat = None ):
        """WARNING: Deprecated method. Please use the sampleScatter method instead.

        Randomly generate scatterings.

        Assuming a scattering took place, generate energy transfer (delta_ekin)
        and new neutron direction based on current kinetic energy and direction
        and return tuple(new_direction,delta_ekin). The repeat parameter can be
        set to a positive number, causing the scattering to be sampled that many
        times and numpy arrays with results returned.

        """
        return _rawfct['ncrystal_genscatter'](self._rawobj_scat,ekin,direction,repeat)

    def generateScatteringNonOriented( self, ekin, repeat = None ):
        """WARNING: Deprecated method. Please use the sampleScatterIsotropic method instead.

        Randomly generate scatterings (should not be called for oriented processes).

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

    def scatter(self,ekin=None,direction=None,wl=None,repeat=None):
        """Convenience function which redirects calls to either
        sampleScatterIsotropic or sampleScatter depending on whether
        or not a direction is given. It can also accept wavelengths instead of
        kinetic energies via the wl parameter.
        """
        ekin = Process._parseekin( ekin, wl )
        return self.sampleScatterIsotropic( ekin, repeat ) if direction is None else self.sampleScatter( ekin, direction, repeat )

    def genscat(self,ekin=None,direction=None,wl=None,repeat=None):
        """WARNING: Deprecated method. Please use the "scatter" method instead.

        Convenience function which redirects calls to either
        generateScatteringNonOriented or generateScattering depending on whether
        or not a direction is given. It can also accept wavelengths instead of
        kinetic energies via the wl parameter.
        """
        ekin = Process._parseekin( ekin, wl )
        return self.generateScatteringNonOriented( ekin, repeat ) if direction is None else self.generateScattering( ekin, direction, repeat )

    def rngSupportsStateManipulation(self):
        """Query whether associated RNG stream supports state manipulation"""
        return bool(_rawfct['ncrystal_rngsupportsstatemanip_ofscatter'](self._rawobj_scat))

    def getRNGState(self):
        """Get current RNG state (as printable hex-string with RNG type info
           embedded). This function returns None if RNG stream does not support
           state manipulation
        """
        return _rawfct['nc_getrngstate_scat'](self._rawobj_scat)

    def setRNGState(self,state):
        """Set current RNG state.

           Note that setting the rng state will affect all objects sharing the
           RNG stream with the given scatter object (and those subsequently cloned
           from any of those).

           Note that if the provided state originates in (the current version
           of) NCrystal's builtin RNG algorithm, it can always be used here,
           even if the current RNG uses a different algorithm (it will simply be
           replaced). Otherwise, a mismatch of RNG stream algorithms will result
           in an error.
        """
        _rawfct['ncrystal_setrngstate_ofscatter']( self._rawobj_scat,
                                                   _str2cstr(state) )


def createInfo(cfgstr):
    """Construct Info object based on provided configuration (using available factories)"""
    return Info(cfgstr)

def createScatter(cfgstr):
    """Construct Scatter object based on provided configuration (using available factories)"""
    return Scatter(cfgstr)

def createScatterIndependentRNG(cfgstr,seed = 0):
    """Construct Scatter object based on provided configuration (using available
    factories) and with its own independent RNG stream (using the builtin RNG
    generator and the provided seed)"""
    rawobj = _rawfct['ncrystal_create_scatter_builtinrng'](_str2cstr(cfgstr),seed)
    return Scatter(('_rawobj_',rawobj))

def createAbsorption(cfgstr):
    """Construct Absorption object based on provided configuration (using available factories)"""
    return Absorption(cfgstr)

def directMultiCreate( data, cfg_params='', *, dtype='',
                       doInfo = True, doScatter = True, doAbsorption = True ):
    """Convenience function which creates Info, Scatter, and Absorption objects
       directly from a data string rather than an on-disk or in-memory
       file. Such usage obviously precludes proper caching behind the scenes,
       and is intended for scenarios where the same data should not be used
       repeatedly.
    """
    if isinstance(data,TextData):
        _localkeepalive = data
        data = data.rawData
    if not dtype and not data.startswith('NCMAT') and 'NCMAT' in data:
        if data.strip().startswith('NCMAT'):
            raise NCBadInput('NCMAT data must have "NCMAT" as the first 5 characters (must not be preceded by whitespace)')

    rawi,raws,rawa = _rawfct['multicreate_direct'](data,dtype,cfg_params,doInfo,doScatter,doAbsorption)
    info = Info( ('_rawobj_',rawi) ) if rawi else None
    scatter = Scatter( ('_rawobj_',raws) ) if raws else None
    absorption = Absorption( ('_rawobj_',rawa) ) if rawa else None
    class MultiCreated:
        def __init__(self,i,s,a):
            self.__i,self.__s,self.__a = i,s,a
        @property
        def info(self):
            """Info object (None if not present)."""
            return self.__i
        @property
        def scatter(self):
            """Scatter object (None if not present)."""
            return self.__s
        @property
        def absorption(self):
            """Absorption object (None if not present)."""
            return self.__a
        def __str__(self):
            fmt = lambda x : str(x) if x else 'n/a'
            return 'MultiCreated(Info=%s, Scatter=%s, Absorption=%s)'%(fmt(self.__i),
                                                                       fmt(self.__s),
                                                                       fmt(self.__a))
    return MultiCreated(info,scatter,absorption)

def registerInMemoryFileData(virtual_filename,data):
    """Register in-memory file data. This needs a "filename" and the content of this
       virtual file. After registering such in-memory "files", they can be used
       as file names in cfg strings or MatCfg objects. Registering the same
       filename more than once, will simply override the content.

       As a special case data can specified as "ondisk://<path>",
       which will instead create a virtual alias for an on-disk file.
    """
    if ( isinstance(data,str) and data.startswith('ondisk://')):
        data = 'ondisk://'+str(pathlib.Path(data[9:]).resolve())
    _rawfct['ncrystal_register_in_mem_file_data'](virtual_filename,data)


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

_const_ekin2ksq_factor = _k4PiSq / _c_wl2ekin
def ekin2ksq(ekin):
    """Convert neutron kinetic energy in electronvolt to squared wavenumber (k^2) in 1/Angstrom^2"""
    return _const_ekin2ksq_factor*ekin

def wl2k(wl):
    """Convert neutron wavelength in Angstrom to wavenumber (k) in 1/Angstrom. This is simply k=2pi/wl"""
    if _np and hasattr(wl,'__len__'):
        #reciprocals without zero division:
        wlnonzero = ekin != 0.0
        ksafe = _k2Pi / _np.where( wlnonzero, wl, 1.0)#fallback 1.0 wont be used
        return  _np.where( wlnonzero, ksafe, _np.inf)
    else:
        return _k2Pi / wl if wl!=0.0 else float('inf')

def wl2ksq(wl):
    """Convert neutron wavelength in Angstrom to wavenumber (k^2) in 1/Angstrom^2. This is simply k^2=(2pi/wl)^2"""
    return wl2k(wl) ** 2

def clearCaches():
    """Clear various caches"""
    _rawfct['ncrystal_clear_caches']()
def clearInfoCaches():
    """Deprecated. Does the same as clearCaches()"""
    clearCaches()

def disableCaching():
    """Obsolete function. Instead call clearCaches() as needed."""
    raise RuntimeError('The disableCaching function has been removed. Users can'
                       +' instead call the clearCaches function if really needed to clear the caches.')

def enableCaching():
    """Obsolete function. Instead call clearCaches() as needed."""
    raise RuntimeError('The enableCaching function has been removed. Users can'
                       +' instead call the clearCaches function if really needed to clear the caches.')

def hasFactory(name):
    """Check if a factory of a given name exists"""
    return bool(_rawfct['ncrystal_has_factory'](_str2cstr(name)))

#Helper function, for scripts creating ncmat files:
def formatVectorForNCMAT(name,values):
    """Utility function for help in python scripts composing .ncmat files,
       transforming an array of of values into a properly formatted text string,
       with word-wrapping, usage of <val>r<n> syntax, etc.
    """
    def provideFormattedEntries():
        def _fmtnum(num):
            _ = '%g'%num if num else '0'#avoid 0.0, -0, etc.
            if _.startswith('0.'):
                _=_[1:]
            return _
        i=0
        v=values.flatten()
        nv=len(v)
        while i<nv:
            fmt_vi=_fmtnum(v[i])
            #check if is repeated:
            irepeat=i
            while irepeat+1<nv:
                if _fmtnum(v[irepeat+1])==fmt_vi:
                    irepeat+=1
                else:
                    break
            yield '%sr%i'%(fmt_vi,1+irepeat-i) if irepeat>i else '%s'%fmt_vi
            i=irepeat+1#advance
    out=''
    line='  %s'%name
    collim=80
    for e in provideFormattedEntries():
        snext=' %s'%e
        line_next=line+snext
        if len(line_next)>collim:
            out += line
            out += '\n'
            line = '   '+snext
        else:
            line = line_next
    if line:
        out += line
        out += '\n'
    return out

#Accept custom random generator:
def setDefaultRandomGenerator(rg, keepalive=True):
    """Set the default random generator.

    Note that this can only changes the random generator for those processes not
    already created.

    To ensure Python does not clean up the passed function object prematurely,
    the NCrystal python module will keep a reference to it eternally. To avoid
    this, call with keepalive=False. But in that case the caller is responsible
    for keeping a reference to the object for as long as NCrystal might use it
    to generate random numbers.

    """
    _rawfct['ncrystal_setrandgen'](rg)

__atomdb={}
def atomDB(Z,A=None,throwOnErrors=True):
    """Access internal database with data for isotopes and natural elements.

    If A is provided, both A and Z must be integers, thus defining a specific isotope.

    If Z is an integer and A is 0 or None, the corresponding natural element is provided.

    Finally, the function can be called with a string identifying either natural
    elements or isotopes: atomDB("Al"), atomDB("He3"), ...

    In all cases, in case of errors or missing entries in the database, either
    an NCBadInput exception is thrown (throwOnErrors==True) or None is
    returned (when throwOnErrors==False).
    """
    global __atomdb
    if isinstance(Z,numbers.Integral):
        Z=int(Z)
        key=(Z,int(A or 0))
        strkey=False
    else:
        assert A is None,"Do not supply two arguments unless the first argument is an integer"
        assert isinstance(Z,str),"The first argument to the function must either be of int or str type"
        key=Z
        strkey=True
    obj=__atomdb.get(key,None)
    if obj:
        return obj
    if strkey:
        rawatomdata=_rawfct['ncrystal_create_atomdata_fromdbstr'](_str2cstr(key))
    else:
        rawatomdata=_rawfct['ncrystal_create_atomdata_fromdb'](*key)
    if not _rawfct['ncrystal_valid'](rawatomdata):
        if not throwOnErrors:
            return None
        if strkey:
            s='key="%s"'%key
        else:
            if key[1]==0:
                s='Z=%i'%key[0]
            else:
                s='Z=%i,A=%i'%key
        raise NCBadInput('atomDB: Could not find entry for key (%s)'%s)
    ad = AtomData(rawatomdata)
    assert ad.isElement()
    Z,A = ad.Z(), (ad.A() if ad.isSingleIsotope() else 0)
    keys=[ (Z,A)]
    if Z==1 and A==2:
        keys+=['H2','D']
    elif Z==1 and A==3:
        keys+=['H3','T']
    else:
        assert ad.isNaturalElement() or ad.isSingleIsotope()
        keys += [ ad.description(False) ]#guaranteed to give just symbol for natelem/singleisotope!
    assert key in keys#Should always be true unless we forgot some keys above
    assert ad.description(False) in keys#Should also be true, given guarantees for AtomData::description(false)
    for k in keys:
        __atomdb[k] = ad
    return ad

def iterateAtomDB(objects=True):
    """Iterate over all entries in the internal database with data for isotopes and
       natural elements. If objects=True, AtomData objects are returned. If
       objects=False, (Z,A) values are returned (A=0 indicates a natural
       element)."""
    for z,a in _rawfct['atomdb_getall_za']():
        yield atomDB(z,a) if objects else (int(z),int(a))


class FileListEntry:
    """Entry in list returned by browseFiles."""
    def __init__(self,*,name,source,factName,priority):
        self.__n = name or None
        self.__f = factName or None
        self.__p = int(priority) if priority.isdigit() else priority
        self.__s = source or None

    @property
    def name(self):
        """The (possibly virtual) filename needed to select this entry"""
        return self.__n

    @property
    def source(self):
        """Description (such as the parent directory in case of on-disk files)"""
        return self.__s

    @property
    def factName(self):
        """Name of the factory delivering entry."""
        return self.__f

    @property
    def priority(self):
        """The priority value of the entry (important in case multiple factories
        delivers content with the the same name). Can be 'Unable',
        'OnlyOnExplicitRequest' or an integer priority value (entries with
        higher values will be preferred).
        """
        return self.__p

    @property
    def fullKey(self):
        """The string '%s::%s'%(self.factName,self.name), which can be used to
           explicitly request this entry without interference from similarly
           named entries in other factories.
        """
        return '%s::%s'%(self.__f,self.__n)

    def __str__(self):
        l=[]
        if self.__n:
            l+=['name=%s'%self.__n]
        if self.__s:
            l+=['source=%s'%self.__s]
        if self.__f:
            l+=['factory=%s'%self.__f]
        l+=['priority=%s'%self.__p]
        return 'FileListEntry(%s)'%(', '.join(l))

    def __lt__(self, other):
        if not isinstance(other, FileListEntry):
            return False
        return ( (self.__f,self.__n,self.__s,self.__p)
                 < (other.__f,other.__n,other.__s,other.__p) )

def browseFiles(dump=False,factory=None):
    """Browse list of available input files (virtual or on-disk). The list is not
       guaranteed to be exhaustive, but will usually include all files in
       supported files in the most obvious locations (the NCrystal data
       directory and other directories of the standard search path, the current
       working directory, virtual files embedded in the NCrystal library or
       registered dynamically.

       Returns a list of FileListEntry objects. If the dump flag is set to True,
       the list will also be printed to stdout in a human readable form.

       Setting factory parameter will only return / print entries from the
       factory of that name.

    """
    res=[]
    def sortkey(e):
        praw = e.priority
        if praw=='Unable':
            p=-2
        elif isinstance(praw,int):
            p=praw
        else:
            assert praw=='OnlyOnExplicitRequest'
            p=-1
        return (-p, e.factName,e.source,e.name)
    for n,s,f,p in _rawfct['ncrystal_get_filelist']():
        res.append( FileListEntry(name=n,source=s,factName=f,priority=p) )
    res.sort(key=sortkey)
    if dump:
        seen_names=set()
        groupfct = lambda e : (e.factName,e.source,e.priority)
        lastgroup = None
        pending=[]
        def print_pending():
            if not pending:
                return
            if factory is not None and lastgroup[0]!=factory:
                pending.clear()
                return
            n=len(pending) - 1
            pending[0] = pending[0]%('%s files'%n if n!=1 else '%s file'%n )
            for line in pending:
                print (line)
            pending.clear()
        for e in res:
            group = groupfct(e)
            if lastgroup != group:
                print_pending()
                lastgroup = group
                pending.append('==> %%s from "%s" (%s, priority=%s):'%group)
            hidden = e.name in seen_names
            seen_names.add(e.name)
            extra=''
            prname=e.name
            if e.priority=='OnlyOnExplicitRequest':
                prname='%s::%s'%(e.factName,e.name)
            elif hidden:
                extra=' <--- Hidden by higher priority entries (select as "%s::%s")'%(e.factName,e.name)
            pending.append(    '    %s%s'%(prname,extra))
        print_pending()
    if factory is None:
        return res
    return [e for e in res if e.factName==factory]

class TextData:
    """Text data accessible line by line, with associated meta-data. This always
       include a UID (useful for comparison and downstream caching purposes) and
       the data type (e.g. "ncmat"). Optionally available is the last known
       on-disk path to a file with the same content, which might be useful in
       case the data needs to be passed to 3rd party software which can only
       work with physical files.

       Text data objects are easily line-iterable, easily providing lines
       (without newline characters): for( auto& line : mytextdata ) {...}.  Of
       course, the raw underlying data buffer can also be accessed if needed.

       The raw data must be ASCII or UTF-8 text, with line endings \\n=CR=0x0A
       (Unix) or \\r\\n=LF+CR=0x0D0A (Windows/DOS). Other encodings might work
       only if 0x00, 0x0A, 0x0D bytes do not occur in them outside of line
       endings.

       Notice that ancient pre-OSX Mac line-endings \\r=LF=0x0D are not
       supported, and iterators will actually throw an error upon encountering
       them. This is done on purpose, since files with \\r on unix might hide
       content when inspected in a terminal can be either confusing, a potential
       security issue, or both.
    """

    def __init__(self,name):
        """create TextData object based on string (same as using createTextData(name))"""
        l=_rawfct['nc_gettextdata'](name)
        assert len(l)==5
        self.__rd = l[0]
        self.__uid = int(l[1])
        self.__dsn = l[2]
        self.__datatype= l[3]
        self.__rp = pathlib.Path(l[4]) if l[4] else None

    @property
    def uid(self):
        """Unique identifier. Objects only have identical UID if all contents and
           metadata are identical."""
        return self.__uid

    @property
    def dataType(self):
        """Data type ("ncmat", "lau", ...)."""
        return self.__datatype

    @property
    def dataSourceName(self):
        """Data source name. This might for instance be a filename."""
        return self.__dsn

    @property
    def rawData(self):
        """Raw access to underlying data."""
        return self.__rd

    @property
    def lastKnownOnDiskLocation(self):
        """Last known on-disk location (returns None if unavailable). Note that there
           is no guarantee against the file having been removed or modified since the
           TextData object was created.
        """
        return self.__rp

    def __str__(self):
        return 'TextData(%s, uid=%i, %i chars)'%(self.__dsn,self.__uid,len(self.__rd))

    def __iter__(self):
        """Line-iteration, yielding lines without terminating newline characters"""
        from io import StringIO
        def chomp(x):
            return x[:-2] if x.endswith('\r\n') else (x[:-1] if x.endswith('\n') else x)
        for l in StringIO(self.__rd):
            yield chomp(l)

def createTextData(name):
    """creates TextData objects based on requested name"""
    return TextData(name)

def getFileContents(name):
    """OBSOLETE FUNCTION: Use createTextData(..).rawData instead."""
    return createTextData(name).rawData

def addCustomSearchDirectory(dirpath):
    """Register custom directories to be monitored for data files."""
    _rawfct['ncrystal_add_custom_search_dir'](_str2cstr(str(pathlib.Path(dirpath).resolve())))

def removeCustomSearchDirectories():
    """Remove all search directories added with addCustomSearchDirectory."""
    _rawfct['ncrystal_remove_custom_search_dirs']()

def removeAllDataSources():
    """Disable all standard data sources, remove all TextData factories as well,
       clear all registered virtual files and custom search directories. Finish
       by calling global clearCaches function ("Ripley: I say we take off and
       nuke the entire site from orbit. It's the only way to be sure.").
    """
    _rawfct['ncrystal_remove_all_data_sources']()

def enableAbsolutePaths( enable = True ):
    """Whether or not absolute file paths are allowed."""
    _rawfct['ncrystal_enable_abspaths'](1 if enable else 0)

def enableRelativePaths( enable = True ):
    """Whether or not paths relative to current working directory are allowed."""
    _rawfct['ncrystal_enable_relpaths'](1 if enable else 0)

def enableStandardSearchPath( enable = True ):
    """Whether or not the standard search path should be searched. This standard
      search path is is by default searched *after* the standard data library,
      and is built by concatenating entries in the NCRYSTAL_DATA_PATH
      environment variables with entries in the compile time definition of the
      same name (in that order). Note that by default the standard search path
      is searched *after* the standard data library.
    """
    _rawfct['ncrystal_enable_stdsearchpath'](1 if enable else 0)

def enableStandardDataLibrary( enable = True, dirpath_override = None ):
    """Whether or not the standard data library shipped with NCrystal should be
       searched.

       Unless NCrystal is configured to have the standard data library embedded
       into the binary at compilation time, the location (directory path) of the
       standard data library is taken from the NCRYSTAL_DATADIR environment
       variable. If the environment variable is not set, the location is taken
       from the compile time definition of the same name. If neither is set, and
       data was not embedded at compilation time, the standard data library will
       be disabled by default and the location must be provided before it can be
       enabled. In all cases, the location can be overridden if explicitly
       provided by the user as the second parameter to this function.
    """
    d = _str2cstr(str(pathlib.Path(dirpath_override).resolve())) if dirpath_override else ctypes.cast(None, ctypes.c_char_p)
    _rawfct['ncrystal_enable_stddatalib'](1 if enable else 0, d)

def browsePlugins(dump=False):

    """Return list of plugins [(pluginname,filename,plugintype),...].

    If the dump flag is set to True, the list will not be returned. Instead it
    will be printed to stdout.
    """
    l=_rawfct['ncrystal_get_pluginlist']()
    if not dump:
        return l
    print('NCrystal has %i plugins loaded.'%len(l))
    for i in range(len(l)):
        pluginname, filename, plugintype = l[i]
        print('==> %s (%s%s)'%(pluginname,plugintype,
                             ' from %s'%filename if filename else ''))

def debyeIsotropicMSD( *, debye_temperature, temperature, mass ):
    """Estimate (isotropic, harmonic) atomic mean-squared-displacement using the
       Debye Model (eq. 11+12 in R.J. Glauber, Phys. Rev. Vol98 num 6,
       1955). Unit of returned MSD value is Aa^2. Input temperatures should be
       in Kelvin, and input atomic mass should be in amu.
    """
    return float(_rawfct['ncrystal_debyetemp2msd'](debye_temperature, temperature, mass))

def debyeTempFromIsotropicMSD( *, msd, temperature, mass ):
    """The inverse of debyeIsotropicMSD (implemented via root-finding), allowing to
       get the Debye temperature which will give rise to a given
       mean-squared-displacement.
    """
    return float(_rawfct['ncrystal_msd2debyetemp'](msd, temperature, mass))

def analyseVDOS(emin,emax,density,temperature,atom_mass_amu):
    """Analyse VDOS curve to extract mean-squared-displacements, Debye temperature,
    effective temperature, gamma0 and integral. Input VDOS must be defined via
    an array of density values, over an equidistant energy grid over [emin,emax]
    (in eV). Additionally, it is required that emin>0, and a parabolic trend
    towards (0,0) will be assumed for energies in [0,emin]. Units are kelvin and
    eV where appropriate.
    """
    return _rawfct['nc_vdoseval'](emin,emax,density,temperature,atom_mass_amu)

def normaliseCfg(cfgstr):
    """Returns normalised version of cfgstr. This is done behind the scenes by
       loading the specified cfg-string into a C++ MatCfg object and then
       re-encoding it as a string.
    """
    return _rawfct['nc_normcfgstr'](cfgstr)

def decodeCfg(cfgstr,*,asJSONStr=False):
    """Decodes cfg-string and returns as Python data structure (a dictionary). The
       format of this data structure should be mostly self-evident by
       inspection, and is not guaranteed to stay the same across NCrystal
       versions. If asJSONStr=true, the data structure will be returned as a
       JSON-encoded string, instead of a Python dictionary.
    """
    _js = _rawfct['nc_cfgstr2json'](cfgstr)
    if asJSONStr:
        return _js
    return json.loads(_js)

def _rawParseNCMAT(text_data_name,*,asJSONStr=False):
    """Parses NCMAT content and returns as Python data structure (a dictionary). The
       format of this data structure should be mostly self-evident by
       inspection, and is not guaranteed to stay the same across NCrystal
       versions. If asJSONStr=true, the data structure will be returned as a
       JSON-encoded string, instead of a Python dictionary.

       WARNING: This function is considered experimental and is currently NOT
       feature complete. It only returns data from a few select NCMAT sections."""
    _js = _rawfct['nc_ncmat2json'](text_data_name)
    if asJSONStr:
        return _js
    return json.loads(_js)

def generateCfgStrDoc( mode = "print" ):
    """Generates documentation about the available cfg-str variables. Mode can
    either be 'print' (print detailed explanation to stdout), 'txt_full' (return
    detailed explanations as string), 'txt_short' (return short explanations as
    string), 'json' (return json-encoded string), or 'python' (return python
    data structures).
    """
    modemap={'print':0,'txt_full':0,'txt_short':1,'json':2,'python':2}
    modeint = modemap.get(mode,None)
    if modeint is None:
        raise NCBadInput('mode must be one of %s'%list(sorted(modemap.keys())))
    _=_rawfct['nc_gencfgdoc'](modeint)
    if mode == 'print':
        print(_)
    else:
        return json.loads(_) if mode=='python' else _

def test():
    """Quick test that NCrystal works as expected in the current installation."""
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

    al = createInfo('stdlib::Al_sg225.ncmat;dcutoff=1.4')
    require(hasFactory('stdncmat'))
    require(al.hasTemperature() and require_flteq(al.getTemperature(),293.15))
    require_flteq(al.getXSectFree(),1.39667)
    require_flteq(al.getXSectAbsorption(),0.231)
    require_flteq(al.getDensity(),2.69864547673)
    require_flteq(al.getNumberDensity(),0.06023238256131625)
    require(al.hasDebyeTemperature())

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
    expected_hkl = { 0  : (1, 1, 1, 8, 2.3380261031049243, 1.7731590373262052),
                     1  : (2, 0, 0, 6, 2.02479, 1.7317882793764163),
                     2  : (2, 2, 0, 12, 1.4317427394787092, 1.5757351418107877) }
    for idx,hkl in enumerate(al.hklList()):
        h,k,l,mult,dsp,fsq = hkl
        require(idx<len(expected_hkl))
        e = expected_hkl[idx]
        require( list(e)[0:4] == [h,k,l,mult] )
        require_flteq(dsp, e[4])
        require_flteq(fsq, e[5])

    #We do all createScatter... here with independent RNG, for reproducibility
    #and to avoid consuming random numbers from other streams.
    alpc = createScatterIndependentRNG('stdlib::Al_sg225.ncmat;dcutoff=1.4;incoh_elas=0;inelas=0')
    require( alpc.name == 'PCBragg' )
    require( isinstance(alpc.name,str) )
    require( alpc.refCount() in (1,2) and type(alpc.refCount()) == int )
    require( alpc.isNonOriented() )
    #print(alpc.xsect(wl=4.0))
    require_flteq(1.632435821586171,alpc.crossSectionIsotropic(wl2ekin(4.0)) )
    require_flteq(1.632435821586171,alpc.crossSection(wl2ekin(4.0),(1,0,0)))
    require( alpc.crossSectionIsotropic(wl2ekin(5.0)) == 0.0 )

    require( alpc.rngSupportsStateManipulation() )
    require(alpc.getRNGState()=='a79fd777407ba03b3d9d242b2b2a2e58b067bd44')

    alpc.setRNGState('deadbeefdeadbeefdeadbeefdeadbeefb067bd44')
    require(alpc.getRNGState()=='deadbeefdeadbeefdeadbeefdeadbeefb067bd44')
    alpc_clone = alpc.clone()
    require(alpc.getRNGState()=='deadbeefdeadbeefdeadbeefdeadbeefb067bd44')
    require(alpc_clone.getRNGState()=='e0fd16d42a2aced7706cffa08536d869b067bd44')
    alpc_clone2 = alpc_clone.clone(for_current_thread=True)
    require(alpc_clone2.getRNGState()=='cc762bb1160a0be514300da860f6d160b067bd44')
    alpc_clone3 = alpc_clone.clone(rng_stream_index = 12345 )
    require(alpc_clone3.getRNGState()=='3a20660a10fd581bd7cddef8fc3f32a2b067bd44')

    #Pick Nickel at 1.2 angstrom, to also both vdos + incoherent-elastic + coherent-elastic:
    nipc = createScatterIndependentRNG('stdlib::Ni_sg225.ncmat;dcutoff=0.6;vdoslux=2',2543577)
    nipc_testwl = 1.2
    #print(nipc.xsect(wl=nipc_testwl),nipc.xsect(wl=5.0))
    require_flteq(16.76322537767633,nipc.xsect(wl=nipc_testwl))
    require_flteq(16.76322537767633,nipc.xsect(wl=nipc_testwl,direction=(1,0,0)))
    require_flteq(5.958094744249944,nipc.xsect(wl=5.0))

    require( nipc.name == 'ProcComposition' )

    expected = [ ( 0.056808478892590906, 0.5361444826572666 ),
                 ( 0.056808478892590906, 0.5361444826572666 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.056808478892590906, 0.8391056916029316 ),
                 ( 0.032001313096074194, -0.3726211530494784 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.056808478892590906, -0.15963879335683306 ),
                 ( 0.056808478892590906, 0.8260541809964751 ),
                 ( 0.07799891342244511, -0.5293689067509328 ),
                 ( 0.05348597239804991, -0.09542895891111686 ),
                 ( 0.056808478892590906, 0.8260541809964751 ),
                 ( 0.04125596596989546, -0.2114086959439411 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.056808478892590906, 0.5361444826572666 ),
                 ( 0.056808478892590906, -0.3915665520281999 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.05750930296126236, -0.5221485562238027 ),
                 ( 0.056808478892590906, 0.36219866365374176 ),
                 ( 0.0812282226184265, -0.9893430983107131 ),
                 ( 0.056808478892590906, -0.5655123710317247 ),
                 ( 0.058097184485747494, -0.951408724433637 ),
                 ( 0.056808478892590906, 0.3042167239859003 ),
                 ( 0.056808478892590906, 0.7378808571510718 ),
                 ( 0.056808478892590906, -0.10165685368899147 ),
                 ( 0.08045270890565197, -0.8062090812584963 ),
                 ( 0.056808478892590906, -0.5655123710317247 ),
                 ( 0.06930621305459778, 0.0790013056899895 ),
                 ( 0.04019454300337846, -0.9619857378392679 ),
                 ( 0.08983663449895618, -0.5822245509723732 ) ]

    if _np is None:
        ekin,mu=[],[]
        for i in range(30):
            _ekin,_mu=nipc.sampleScatterIsotropic(wl2ekin(nipc_testwl))
            mu += [_mu]
            ekin += [_ekin]
    else:
        ekin,mu = nipc.sampleScatterIsotropic(wl2ekin(nipc_testwl),repeat=30)

    for i in range(len(ekin)):
        #print ( f'    ( {ekin[i]}, {mu[i]} ),');continue
        require_flteq(ekin[i],expected[i][0])
        require_flteq(mu[i],expected[i][1])

    expected = [ ( 0.056808478892590906, (0.07228896531453344, -0.5190173207165885, 0.8517014302500192) ),
                 ( 0.056808478892590906, (-0.9249112255344181, -0.32220112076758217, -0.20180600252850442) ),
                 ( 0.056808478892590906, (-0.15963879335683306, -0.8486615569734178, 0.5042707778277745) ),
                 ( 0.0492224669270449, (-0.9779916301402904, 0.140975569549056, 0.15381241876342955) ),
                 ( 0.056808478892590906, (0.07228896531453344, 0.7905105193171594, -0.6081672667471253) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.8869759070713323, -0.4504882066969593) ),
                 ( 0.056808478892590906, (0.07228896531453344, -0.39741541395284924, -0.914787021249449) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.9768880366798581, -0.1880309758785167) ),
                 ( 0.025610418737037184, (-0.884783168996826, -0.4657283985847863, 0.015993830422524287) ),
                 ( 0.056808478892590906, (0.8260541809964751, 0.539797243436807, 0.16202909009269678) ),
                 ( 0.07443181306716587, (-0.6036941700256581, -0.06282145095069493, -0.7947369466543514) ),
                 ( 0.056808478892590906, (0.8260541809964751, 0.10854661864786977, 0.5530389874487663) ),
                 ( 0.056808478892590906, (0.5361444826572666, 0.7795115518549294, 0.3238994199452849) ),
                 ( 0.056808478892590906, (0.07228896531453344, 0.746175597107444, 0.6618128767069312) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.4247181868490453, 0.8996001033001911) ),
                 ( 0.056808478892590906, (0.5361444826572666, 0.5555760611065321, -0.6355189486093415) ),
                 ( 0.05736918247226906, (-0.17265143322057086, -0.684984059616355, 0.7078052844380157) ),
                 ( 0.056808478892590906, (0.3042167239859003, -0.8706122815482211, -0.3866347631352975) ),
                 ( 0.056808478892590906, (-0.7384733804796917, 0.6322144258925643, -0.23443972789660028) ),
                 ( 0.056808478892590906, (-0.15963879335683306, 0.21525619037302965, -0.9634211063505222) ),
                 ( 0.056808478892590906, (0.41359447569500096, 0.4927058865194684, 0.7656242675514158) ),
                 ( 0.056808478892590906, (0.25796367721315083, 0.48520231047621615, 0.8354839670198411) ),
                 ( 0.056808478892590906, (0.5785005938702705, 0.8104481067271115, -0.09225469740985966) ),
                 ( 0.04320288783696474, (-0.030385860971878235, -0.49547867364837517, 0.8680884651996275) ),
                 ( 0.05428746831317616, (-0.3602629693021255, -0.9063804616575692, 0.22064689364463652) ),
                 ( 0.056808478892590906, (0.36219866365374176, -0.8822186430862216, 0.3008361577978114) ),
                 ( 0.056808478892590906, (0.7680722413286334, 0.5975216576265994, -0.23028873347945303) ),
                 ( 0.056808478892590906, (0.32922859149927786, -0.9426419619170849, 0.0550878042084668) ),
                 ( 0.056808478892590906, (-0.10165685368899147, -0.2489220191768986, -0.9631737706493833) ),
                 ( 0.0670456981700729, (-0.8979842133391931, 0.34668021323231085, 0.2709929562678522) ) ]

    for i in range(30):
        out_ekin,outdir = nipc.sampleScatter(wl2ekin(nipc_testwl),(1.0,0.0,0.0))
        #print ( f'    ( {out_ekin}, {outdir} ),');continue
        require_flteq(out_ekin,expected[i][0])
        require_flteq(outdir[0],expected[i][1][0])
        require_flteq(outdir[1],expected[i][1][1])
        require_flteq(outdir[2],expected[i][1][2])
    gesc = createScatterIndependentRNG("""stdlib::Ge_sg227.ncmat;dcutoff=0.5;mos=40.0arcsec
                            ;dir1=@crys_hkl:5,1,1@lab:0,0,1
                            ;dir2=@crys_hkl:0,-1,1@lab:0,1,0""",3453455)
    require_flteq(591.025731949681,gesc.crossSection(wl2ekin(1.540),( 0., 1., 1. )))
    require_flteq(1.666984885615526,gesc.crossSection(wl2ekin(1.540),( 1., 1., 0. )))

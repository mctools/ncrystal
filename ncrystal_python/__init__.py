#!/usr/bin/env python
"""Python module for using the NCrystal library for thermal neutron transport in crystals

Please find more information about NCrystal at the website:

   https://mctools.github.io/ncrystal/

In particular, a small example using the NCrystal python module can be found at:

   https://github.com/mctools/ncrystal/blob/master/examples/ncrystal_example_py

A substantial effort went into developing NCrystal. If you use it for your work,
we would appreciate it if you would use the following reference in your work
(note, this reference will change once a proper peer-reviewed publication for
NCrystal exists):

  X. X. Cai and T. Kittelmann, NCrystal, https://doi.org/10.5281/zenodo.853186

For detailed usage conditions and licensing of this open source project, see:

   https://github.com/mctools/ncrystal/blob/master/NOTICE
   https://github.com/mctools/ncrystal/blob/master/LICENSE
   https://github.com/mctools/ncrystal/blob/master/ncrystal_extra/LICENSE

"""

################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2017 NCrystal developers                                   ##
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
__copyright__ = "Copyright 2017"
__version__ = '0.9.9'
__status__ = "Production"
__author__ = "NCrystal developers (Thomas Kittelmann, Xiao Xiao Cai)"
__copyright__ = "Copyright 2015-2017 %s"%__author__
__maintainer__ = __author__
__email__ = "ncrystal-developers@cern.ch"
__all__ = ['NCException','NCFileNotFound','NCDataLoadError','NCCalcError','NCLogicError','NCBadInput',
           'RCBase','Info','CalcBase','Process','Absorption','Scatter',
           'createInfo','createScatter','createAbsorption',
           'test','setDefaultRandomGenerator','RandomCtxMgr','wl2ekin','ekin2wl',
           'clearInfoCaches','disableCaching','enableCaching','clearFactoryRegistry','hasFactory',
           'version_num']

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
    #converts any string (str,bytes,unicode) to bytes
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
                     '/usr/lib' ]

    #Look inside searched directories for libNCrystal:
    for d in searchpaths:
        for f in sorted(glob.glob(os.path.join(d,'libNCrystal.*'))):
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

_globalstates = {}

def _load(nclib_filename):
    import ctypes

    _nclib = ctypes.CDLL(nclib_filename)
    _int,_intp,_uint,_uintp,_dbl,_dblp,_cstr,_voidp = (ctypes.c_int, ctypes.POINTER(ctypes.c_int),
                                                       ctypes.c_uint,ctypes.POINTER(ctypes.c_uint), ctypes.c_double,
                                                       ctypes.POINTER(ctypes.c_double), ctypes.c_char_p, ctypes.c_void_p)
    ndarray_to_dblp = lambda a : a.ctypes.data_as(_dblp)

    _npempty = [None]
    def _create_numpy_array(n):
        if not _np:
            raise NCException("numpy not available - array versions of crossSection/genScatter methods are unavailable")
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
        assert _ncerror()
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
        print("WARNING: Version mismatch detected between NCrystal python code (v%s) and loaded binary library (v%s). Control which NCrystal library to load with the NCRYSTAL_LIB env var."%(__version__,lib_version))

    _wrap('ncrystal_sethaltonerror',_int,(_int,),hide=True,error_check=False)(False)
    _wrap('ncrystal_setquietonerror',_int,(_int,),hide=True,error_check=False)(True)
    _ncerror       = _wrap('ncrystal_error',_int,tuple(),hide=True,error_check=False)
    _ncerror_msg   = _wrap('ncrystal_lasterror',_cstr,tuple(),hide=True,error_check=False)
    _ncerror_type  = _wrap('ncrystal_lasterrortype',_cstr,tuple(),hide=True,error_check=False)
    _ncerror_clear = _wrap('ncrystal_clearerror',None,tuple(),hide=True,error_check=False)

    _wrap('ncrystal_refcount',_int,(_voidp,),take_ref=True)
    #_wrap('ncrystal_ref',None,(_voidp,),take_ref=True)
    _wrap('ncrystal_unref',None,(_voidp,),take_ref=True)
    #_wrap('ncrystal_unrefnodelete',None,(_voidp,),take_ref=True)
    #_wrap('ncrystal_valid',_int,(_voidp,),take_ref=True)
    #_wrap('ncrystal_invalidate',None,(_voidp,),take_ref=True)

    _wrap('ncrystal_cast_scat2proc',ncrystal_process_t,(ncrystal_scatter_t,))
    _wrap('ncrystal_cast_abs2proc',ncrystal_process_t,(ncrystal_absorption_t,))


    _wrap('ncrystal_dump',None,(ncrystal_info_t,))
    _wrap('ncrystal_ekin2wl',_dbl,(_dbl,))
    _wrap('ncrystal_wl2ekin',_dbl,(_dbl,))
    _wrap('ncrystal_isnonoriented',_int,(ncrystal_process_t,))
    _wrap('ncrystal_name',_cstr,(ncrystal_process_t,))

    for s in ('temperature','xsectabsorption','xsectfree','debyetemp','density'):
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
    functions['ncrystal_info_gethkl_setuppars'] = lambda : (_int(),_int(),_int(),_int(),_dbl(),_dbl())

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
    _wrap('ncrystal_setsimplerandgen',None,tuple())

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
    _wrap('ncrystal_clear_info_caches',None,tuple())
    _wrap('ncrystal_disable_caching',None,tuple())
    _wrap('ncrystal_enable_caching',None,tuple())
    _wrap('ncrystal_clear_factory_registry',None,tuple())
    _wrap('ncrystal_has_factory',_int,(_cstr,))

    return functions

_rawfct = _load(_find_nclib())

def decodecfg_packfact(cfgstr):
    return _rawfct['ncrystal_decodecfg_packfact'](_str2cstr(cfgstr))

class RCBase:
    """Base class for all NCrystal objects"""
    def __init__(self, rawobj):
        """internal usage only"""
        self._rawobj = rawobj
        #do not ref here, since ncrystal_create_xxx functions in C-interface already did so.
        self._rawunref = _rawfct['ncrystal_unref']#keep fct reference
    def __del__(self):
        if hasattr(self,'_rawunref'):
            self._rawunref(self._rawobj)
    def refCount(self):
        """Access reference count of wrapped C++ object"""
        return _rawfct['ncrystal_refcount'](self._rawobj)

def nc_assert(b,msg=""):
    """Assertion which throws NCLogicError on failure"""
    if not bool(b):
        raise NCLogicError(msg if msg else 'assertion failed')

class Info(RCBase):
    """Class representing information about a given crystal"""
    def __init__(self, cfgstr):
        """create Info object based on cfg-string (same as using createInfo(cfgstr))"""
        rawobj = _rawfct['ncrystal_create_info'](_str2cstr(cfgstr))
        super(Info, self).__init__(rawobj)
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

    def hasDebyeTemperature(self):
        return _rawfct['ncrystal_info_getdebyetemp'](self._rawobj)>-1
    def getDebyeTemperature(self):
        t=_rawfct['ncrystal_info_getdebyetemp'](self._rawobj)
        nc_assert(t>-1)
        return t

    def hasDensity(self):
        return _rawfct['ncrystal_info_getdensity'](self._rawobj)>-1
    def getDensity(self):
        t=_rawfct['ncrystal_info_getdensity'](self._rawobj)
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
        """Domain of has non-vanishing cross-section.

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
    """Clear all cached Info objects held by factory infrastructure"""
    _rawfct['ncrystal_clear_info_caches']()
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

    al = createInfo('Al_sg225.ncmat;dcutoff=1.4')
    require(hasFactory('stdncmat'))
    require(al.hasTemperature() and flteq(al.getTemperature(),293.15))
    require(al.hasXSectFree() and flteq(al.getXSectFree(),1.39667))
    require(al.hasXSectAbsorption() and flteq(al.getXSectAbsorption(),0.231))
    require(al.hasDensity() and flteq(al.getDensity(),2.69864547673))
    require(not al.hasDebyeTemperature())
    #TODO for NC2: expose and check hasPerElementDebyeTemperature()

    require(al.hasStructureInfo())
    si=al.getStructureInfo()
    require( si['spacegroup'] == 225 )
    require( flteq(si['a'],4.04958) )
    require( flteq(si['b'],4.04958) )
    require( flteq(si['c'],4.04958) )
    require( si['alpha'] == 90.0 )
    require( si['beta'] == 90.0 )
    require( si['gamma'] == 90.0 )
    require( si['n_atoms'] == 4 )
    require( flteq(si['volume'],66.4094599932) )
    require( al.hasHKLInfo() )
    require( al.nHKL() == 3 )
    require( flteq(al.hklDLower(),1.4) )
    require( al.hklDUpper() > 1e36 )
    expected_hkl = { 0  : (1, -1, -1, 8, 2.3380261031049243, 1.7720759166647517),
                     1  : (0, 0, 2, 6, 2.02479, 1.730377956791613),
                     2  : (0, 2, -2, 12, 1.4317427394787094, 1.5731697127736115) }
    for idx,hkl in enumerate(al.hklList()):
        h,k,l,mult,dsp,fsq = hkl
        require(idx<len(expected_hkl))
        e = expected_hkl[idx]
        require( list(e)[0:4] == [h,k,l,mult] )
        require( flteq(dsp, e[4]) )
        require( flteq(fsq, e[5]) )
    #TODO for NC2: Missing is atominfo, bkgd cross-sections, ...

    alpc = createScatter('Al_sg225.ncmat;dcutoff=1.4;bkgd=none')
    require( alpc.name == 'PCBragg' )
    require( isinstance(alpc.name,str) )
    require( alpc.refCount() in (1,2) and type(alpc.refCount()) == int )
    require( alpc.isNonOriented() )
    require( flteq(1.63130945209,alpc.crossSectionNonOriented(wl2ekin(4.0)) ) )
    require( flteq(1.63130945209,alpc.crossSection(wl2ekin(4.0),(1,0,0))))
    require( alpc.crossSectionNonOriented(wl2ekin(5.0)) == 0.0 )

    #use simple fall-back rng for guaranteed reproducibility (set explicitly to avoid warning):
    _rawfct['ncrystal_setsimplerandgen']()
    alpc.generateScatteringNonOriented(wl2ekin(4.0))

    alpc = createScatter('Al_sg225.ncmat;dcutoff=1.4')
    require( alpc.name == 'ScatterComp' )
    expected  = [(2.8283092712311082, 0.0), (2.8283092712311082, 0.0), (2.8283092712311082, 0.0),
                 (2.0527318521221001, 0.0), (2.0527318521221001, 0.0), (2.0527318521221001, 0.0),
                 (2.8283092712311082, 0.0), (0.72997630752329012, 0.038570589541834406),
                 (2.0527318521221001, 0.0), (2.0527318521221001, 0.0)]
    if _np is None:
        ang,de=[],[]
        for i in range(10):
            _ang,_de=alpc.generateScatteringNonOriented(wl2ekin(4.0))
            ang += [_ang]
            de += [_de]
    else:
        ang,de = alpc.generateScatteringNonOriented(wl2ekin(4.0),repeat=10)
    for i in range(10):
        require(flteq(ang[i],expected[i][0]))
        require(flteq(de[i],expected[i][1]))

    expected = [((0.6341262223410173, 0.13834601785989772, -0.7607524653143226), 0.0),
                ((0.6341262223410173, -0.5926699932081182, 0.49661475339562755), 0.0),
                ((0.2990135239748528, -0.13526863370357844, -0.9446127826872275), 0.01260216452284341),
                ((0.5121682964546899, 0.25499685307902603, -0.8201586682017661), 0.0),
                ((-0.6102623490574222, -0.4411233110504387, 0.6580198247551626), -0.00931496091586174),
                ((0.6341262223410173, -0.6188438264965689, 0.46359060877739444), 0.0),
                ((-0.3161005332440677, -0.8400696255604475, -0.440866733938438), -0.008399889907757894),
                ((0.6341262223410173, 0.1298845525938583, -0.7622427022523759), 0.0),
                ((0.6341262223410173, 0.5449549585211818, 0.5485508429696264), 0.0),
                ((0.6341262223410173, 0.5391607731500441, -0.554246871741968), 0.0)]

    for i in range(10):
        ang,de = alpc.generateScattering(wl2ekin(2.0),(1.0,0.0,0.0))
        require(flteq(de,expected[i][1]))
        require(flteq(ang[0],expected[i][0][0]))
        require(flteq(ang[1],expected[i][0][1]))
        require(flteq(ang[2],expected[i][0][2]))

    gesc = createScatter("""Ge_sg227.ncmat;dcutoff=0.5;mos=40.0arcsec
                            ;dir1=@crys_hkl:5,1,1@lab:0,0,1
                            ;dir2=@crys_hkl:0,-1,1@lab:0,1,0""")
    require(flteq(587.78062659,gesc.crossSection(wl2ekin(1.540),( 0., 1., 1. ))))
    require(flteq(1.76682279301,gesc.crossSection(wl2ekin(1.540),( 1., 1., 0. ))))

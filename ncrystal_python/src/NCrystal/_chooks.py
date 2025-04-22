
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


"""Internal module providing ctypes-based hooks into the compiled NCrystal
shared library"""

__all__ = ['_get_raw_cfcts','_str2cstr','_cstr2str']

_rawfcts = [None]
_namespace = [None]
def _get_raw_cfcts():
    if _rawfcts[0] is None:
        from ._locatelib import get_libpath_and_namespace
        _rawfcts[0] = False
        thelib, namespace = get_libpath_and_namespace()
        _rawfcts[0] = _load(thelib, namespace)
        assert _rawfcts[0] is not None
        _namespace[0] = namespace or ''
    return _rawfcts[0]

def _get_build_namespace():
    if _namespace[0] is None:
        _get_raw_cfcts()
    return _namespace[0]

def _str2cstr(s):
    #converts any string (str,bytes,unicode,path) to bytes
    if hasattr(s,'__fspath__'):
        s = str(s)
    try:
        return s if isinstance(s,bytes) else s.encode('utf8')
    except UnicodeEncodeError as e:
        from .exceptions import NCBadInput
        raise NCBadInput("Only unicode strings are supported") from e

def _cstr2str(s):
    #converts bytes object to str
    try:
        return s if isinstance(s,str) else s.decode('utf8')
    except UnicodeDecodeError as e:
        from .exceptions import NCBadInput
        raise NCBadInput("Only UTF8-encoded C-strings are supported") from e

_keepalive = []#for python based callback functions which we need to keep alive

def _load(nclib_filename, ncrystal_namespace_protection ):

    import ctypes
    from ._numpy import _np, _ensure_numpy

    try:
        _nclib = ctypes.CDLL(nclib_filename)
    except TypeError:
        _nclib = None

    if _nclib is None:
        #For some reason, on windows we get a TypeError and must pass a string
        #rather than a pathlib object:
        _nclib = ctypes.CDLL(str(nclib_filename))

    _int,_intp,_uint,_uintp,_dbl,_dblp,_cstr,_voidp = (ctypes.c_int, ctypes.POINTER(ctypes.c_int),
                                                       ctypes.c_uint,ctypes.POINTER(ctypes.c_uint), ctypes.c_double,
                                                       ctypes.POINTER(ctypes.c_double), ctypes.c_char_p, ctypes.c_void_p)
    _ulong = ctypes.c_ulong
    _charptr = ctypes.POINTER(ctypes.c_char)
    _charptrptr = ctypes.POINTER(_charptr)

    _cstrp = ctypes.POINTER(_cstr)
    _cstrpp = ctypes.POINTER(_cstrp)
    _dblpp = ctypes.POINTER(_dblp)
    def ndarray_to_dblp(a):
        return a.ctypes.data_as(_dblp)
    def ndarray_to_uintp(a):
        return a.ctypes.data_as(_uintp)
    def ndarray_to_intp(a):
        return a.ctypes.data_as(_intp)

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
    from .exceptions import NCException, NCFileNotFound, NCDataLoadError, NCMissingInfo, NCCalcError, NCLogicError, NCBadInput, nc_assert
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
    ncrystal_namespace_prefix = 'ncrystal_' if not ncrystal_namespace_protection else 'ncrystal%s_'%ncrystal_namespace_protection
    def _wrap(fct_name,restype,argtypes,take_ref = False, hide=False, error_check=True):
        assert isinstance(argtypes,tuple)
        assert fct_name.startswith('ncrystal_')
        if ncrystal_namespace_protection:
            fct_name_actual = ncrystal_namespace_prefix + fct_name[9:]
        else:
            fct_name_actual = fct_name

        raw=getattr(_nclib,fct_name_actual)
        raw.argtypes=argtypes
        raw.restype=restype

        if take_ref:
            assert len(argtypes)==1
            def fct(arg):
                return raw(ctypes.byref(arg))
        else:
            def fct(*args):
                return raw(*args)
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
    from . import __version__ as _nc_version
    if lib_version != _nc_version:
        raise NCException("ERROR: Version mismatch detected between NCrystal python code (v%s)"
                          " and loaded binary"" library (v%s). Control which NCrystal library"
                          " to load with the NCRYSTAL_LIB env var."%(_nc_version,lib_version))

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
        h,k,ll,mult,dsp,fsq = _int(),_int(),_int(),_int(),_dbl(),_dbl()
        nhkl = int(functions['ncrystal_info_nhkl'](nfo))
        for idx in range(nhkl):
            functions['ncrystal_info_gethkl'](nfo,idx,h,k,ll,mult,dsp,fsq)
            if not all_indices:
                yield h.value,k.value,ll.value,mult.value,dsp.value,fsq.value
            else:
                nc_assert( mult.value % 2 == 0 )
                n = mult.value // 2
                hvals, hvalsptr = _create_numpy_int_array( n )
                kvals, kvalsptr = _create_numpy_int_array( n )
                lvals, lvalsptr = _create_numpy_int_array( n )
                _raw_gethkl_allindices(nfo,idx,hvalsptr,kvalsptr,lvalsptr)
                yield hvals,kvals,lvals,mult.value,dsp.value,fsq.value
    functions['iter_hkllist']=iter_hkllist

    _raw_dealloc_dblptr = _wrap('ncrystal_dealloc_doubleptr',None,(_dblp,),hide=True)
    def _cptr_to_nparray( cptr, n, free_after_with_ncrystal_dealloc = True ):
        import numpy.ctypeslib as _np_ctypeslib
        np_arr = _np_ctypeslib.as_array(cptr,shape=(int(n.value if hasattr(n,'value') else n),)).copy()
        if free_after_with_ncrystal_dealloc:
            _raw_dealloc_dblptr(cptr)
        return np_arr

    _wrap('ncrystal_info_ndyninfo',_uint,(ncrystal_info_t,))
    _raw_di_base = _wrap('ncrystal_dyninfo_base',None,(ncrystal_info_t,_uint,_dblp,_uintp,_dblp,_uintp),hide=True)
    _raw_di_scatknl = _wrap('ncrystal_dyninfo_extract_scatknl',None,(ncrystal_info_t,_uint,_uint,_dblp,_uintp,_uintp,_uintp,
                                                                     _dblpp,_dblpp,_dblpp,_dblpp),hide=True)
    _raw_di_vdos = _wrap('ncrystal_dyninfo_extract_vdos',None,(ncrystal_info_t,_uint,_dblp,_dblp,_uintp,_dblpp),hide=True)
    _raw_di_vdosdebye = _wrap('ncrystal_dyninfo_extract_vdosdebye',None,(ncrystal_info_t,_uint,_dblp),hide=True)
    _raw_di_vdos_input = _wrap('ncrystal_dyninfo_extract_vdos_input',None,(ncrystal_info_t,_uint,_uintp,_dblpp,_uintp,_dblpp),hide=True)
    _raw_vdos2gn = _wrap('ncrystal_raw_vdos2gn',None,(_dblp,_dblp,_uint,_uint,_dbl,_dbl,_dbl,_uint,_dblp,_dblp,_uintp,_dblpp),hide=True)

    def raw_vdos2gn( egrid, density, scatxs, mass_amu, temperature, nvalue ):
        _ensure_numpy()
        _egrid = _np.asarray(egrid,dtype=float)
        _density = _np.asarray(density,dtype=float)
        _s = _dbl(float(scatxs))
        _m = _dbl(float(mass_amu))
        _t = _dbl(float(temperature))
        _n = _uint(int(nvalue))
        xmin, xmax, ny = _dbl(), _dbl(), _uint()
        yraw = _dblp()
        _raw_vdos2gn( ndarray_to_dblp(_egrid),
                      ndarray_to_dblp(_density),
                      _uint(len(_egrid)),
                      _uint(len(_density)),
                      _s,_m,_t,_n, xmin, xmax, ny, yraw )
        return float(xmin.value), float(xmax.value), _cptr_to_nparray( yraw, ny )
    functions['raw_vdos2gn'] = raw_vdos2gn

    _ORDERWEIGHTFCTTYPE = ctypes.CFUNCTYPE( _dbl, _uint )
    _raw_vdos2knl = _wrap('ncrystal_raw_vdos2kernel',None,
                          (_dblp,_dblp,_uint,_uint,_dbl,_dbl,_dbl,_uint,
                           _ORDERWEIGHTFCTTYPE,_uintp,_uintp,_dblpp,_dblpp,
                           _dblpp,_dbl,_dblp),
                          hide=True)
    def raw_vdos2knl( egrid, density, scatxs, mass_amu, temperature,
                      vdoslux, order_weight_fct, target_emax ):
        _ensure_numpy()
        _egrid = _np.asarray(egrid,dtype=float)
        _density = _np.asarray(density,dtype=float)
        _s = _dbl(float(scatxs))
        _m = _dbl(float(mass_amu))
        _t = _dbl(float(temperature))
        _vdl = _uint(int(vdoslux))
        _tgtemax = target_emax or 0.0
        _tgtemax = _dbl( float(_tgtemax if _tgtemax>0.0 else 0.0 ) )
        nalpha, nbeta, suggest_emax = _uint(), _uint(), _dbl()
        agrid, bgrid, sab = _dblp(), _dblp(), _dblp()
        if order_weight_fct:
            def owf(order):
                return float(order_weight_fct( int(order.value
                                                   if hasattr(order,'value')
                                                   else order) ))
            #NB: No need to keep owf alive, after our function exists, since it is only used during evaluation.
            _owf = _ORDERWEIGHTFCTTYPE(owf)
        else:
            _owf = ctypes.cast(None, _ORDERWEIGHTFCTTYPE)
        _raw_vdos2knl( ndarray_to_dblp(_egrid), ndarray_to_dblp(_density),
                       _uint(len(_egrid)),_uint(len(_density)),
                       _s,_m,_t,_vdl,_owf, nalpha, nbeta,
                       ctypes.byref(agrid), ctypes.byref(bgrid), ctypes.byref(sab),
                       _tgtemax, suggest_emax )
        return ( _cptr_to_nparray( agrid, nalpha ),
                 _cptr_to_nparray( bgrid, nbeta ),
                 _cptr_to_nparray( sab, nalpha.value * nbeta.value ),
                 float( suggest_emax.value ) or None )
    functions['raw_vdos2knl'] = raw_vdos2knl

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
        if _np is None and repeat is not None:
            raise NCBadInput('Can not use "repeat" parameter when Numpy is absent on the system')
        if repeat is None and not hasattr(ekin,'__len__'):
            return None#scalar case, array interface not triggered
        repeat = 1 if repeat is None else repeat
        ekin = (ekin if hasattr(ekin,'ctypes') else _np.asarray(ekin,dtype=float) ) if hasattr(ekin,'__len__') else _np.ones(1)*ekin
        #NB: returning the ekin object itself is important in order to keep a reference to it after the call:
        return ndarray_to_dblp(ekin),len(ekin),repeat,ekin

    _raw_xs_no = _wrap('ncrystal_crosssection_nonoriented',None,(ncrystal_process_t,_dbl,_dblp),hide=True)
    _raw_xs_no_many = _wrap('ncrystal_crosssection_nonoriented_many',None,(ncrystal_process_t,_dblp,_ulong,
                                                                           _ulong,_dblp),hide=True)
    _empty_arrayf = _np.empty(shape=(0,),dtype=float) if _np else tuple()
    _empty_arrayf_2tuple = ( ( _np.empty(shape=(0,),dtype=float),
                               _np.empty(shape=(0,),dtype=float) )
                             if _np else (tuple(),tuple()) )

    def ncrystal_crosssection_nonoriented(scat,ekin,repeat=None):
        many = _prepare_many(ekin,repeat)
        if many is None:
            res = _dbl()
            _raw_xs_no(scat,ekin,res)
            return res.value
        else:
            ekin_ct,n_ekin,repeat,ekin_nparr = many
            if repeat == 0:
                return _empty_arrayf
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
            if repeat == 0:
                #special case which happens often in our nctool usage:
                return _empty_arrayf_2tuple
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

    _raw_dump_tostr = _wrap('ncrystal_dump_tostr',_charptr,(ncrystal_info_t,_uint),hide=True)
    def nc_dump_tostr( raw_infoobj, verbosity ):
        return _decode_and_dealloc_raw_str( _raw_dump_tostr( raw_infoobj, verbosity ) )
    functions['nc_dump_tostr'] = nc_dump_tostr

    _raw_dbg_process = _wrap('ncrystal_dbg_process',_charptr,(ncrystal_process_t,),hide=True)
    def nc_dbg_proc(rawprocobj):
        import json
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
    _raw_deallocstrlist = _wrap('ncrystal_dealloc_stringlist',None,(_uint,_cstrp),hide=True)

    def nc_gettextdata(name):
        ll = _raw_gettextdata(_str2cstr(str(name)))
        assert ll is not None
        n = 5
        def _decode( s ):
            s=s.decode('utf-8')
            return ( s.replace('\r\n','\n').replace('\r','\n')
                     if '\r' in s else s )
        res = [_decode(ll[i]) for i in range(n)]
        assert isinstance(res[0],str)
        _raw_deallocstrlist(n,ll)
        return res
    functions['nc_gettextdata'] = nc_gettextdata

    _raw_getfilelist = _wrap('ncrystal_get_file_list',None,(_uintp,_cstrpp),hide=True)
    def ncrystal_get_filelist():
        n,ll = _uint(),_cstrp()
        _raw_getfilelist(n,ctypes.byref(ll))
        assert n.value%4==0
        res=[]
        for i in range(n.value//4):
            res += [ (ll[i*4].decode(),ll[i*4+1].decode(),
                      ll[i*4+2].decode(),ll[i*4+3].decode()) ]
        _raw_deallocstrlist(n,ll)
        return res
    functions['ncrystal_get_filelist'] = ncrystal_get_filelist

    _raw_getpluginlist = _wrap('ncrystal_get_plugin_list',None,(_uintp,_cstrpp),hide=True)
    def ncrystal_get_pluginlist():
        n,ll = _uint(),_cstrp()
        _raw_getpluginlist(n,ctypes.byref(ll))
        assert n.value%3==0
        res=[]
        for i in range(n.value//3):
            pluginname,filename,plugintype=ll[i*3].decode(),ll[i*3+1].decode(),ll[i*3+2].decode()
            res+=[(pluginname,filename,plugintype)]
        _raw_deallocstrlist(n,ll)
        return res
    functions['ncrystal_get_pluginlist'] = ncrystal_get_pluginlist

    _wrap('ncrystal_add_custom_search_dir',None,(_cstr,))
    _wrap('ncrystal_remove_custom_search_dirs',None,tuple())
    _wrap('ncrystal_enable_abspaths',None,(_int,))
    _wrap('ncrystal_enable_relpaths',None,(_int,))
    _wrap('ncrystal_enable_stddatalib',None,(_int,_cstr))
    _wrap('ncrystal_enable_stdsearchpath',None,(_int,))
    _wrap('ncrystal_remove_all_data_sources',None,tuple())
    _wrap('ncrystal_enable_factory_threadpool',None,(_uint,))

    _raw_benchloadcfg = _wrap('ncrystal_benchloadcfg',_dbl,(_cstr,_int,_int),hide=True)
    def ncrystal_benchloadcfg( cfgstr, do_scatter, nrepeat ):
        arg_doscatter = _int( 1 if do_scatter else 0)
        arg_nrepeat = _int( int(nrepeat) )
        res=_raw_benchloadcfg(_str2cstr(cfgstr),arg_doscatter,arg_nrepeat)
        return float(res)
    functions['benchloadcfg'] = ncrystal_benchloadcfg

    _MSGHANDLERFCTTYPE = ctypes.CFUNCTYPE( None, _cstr, _uint )
    _raw_setmsghandler    = _wrap('ncrystal_setmsghandler',None,(_MSGHANDLERFCTTYPE,),hide=True)
    def ncrystal_setmsghandler(pyhandler):
        #Set msg handler function, keeping references as needed (otherwise fct ptrs
        #kept on C++ side will suddenly stop working!) and casting None to a null-ptr.
        def handler( msg, msgtype ):
            #NB: We could instead consider converting to bytes rather than str
            #in case msgtype==2 (raw output):
            pyhandler( _cstr2str(msg), int(msgtype) )
        if not handler:
            keepalive=(None,None,ctypes.cast(None, _MSGHANDLERFCTTYPE))
        else:
            keepalive=(pyhandler,handler,_MSGHANDLERFCTTYPE(handler))#keep refs!
        _keepalive.append(keepalive)
        _raw_setmsghandler(keepalive[-1])
    functions['setmsghandler'] = ncrystal_setmsghandler

    _raw_runmmcsim_stdengine    = _wrap('ncrystal_runmmcsim_stdengine',None,
                                        (_uint,_uint,
                                         _cstr,_cstr,_cstr,
                                         _charptrptr,_uintp,_dblpp,_dblpp),
                                        hide=True)
    def nc_runmmcsim_stdengine( nthreads,
                                tally_detail_lvl,
                                mat_cfgstr,
                                mmc_geomcfg,
                                mmc_srccfg ):
        nthreads_ = _uint(int(nthreads))
        assert 0<=tally_detail_lvl<=2
        tally_detail_lvl_ = _uint(int(tally_detail_lvl))
        mat_cfgstr_ = _str2cstr(mat_cfgstr)
        mmc_geomcfg_ = _str2cstr(mmc_geomcfg)
        mmc_srccfg_ = _str2cstr(mmc_srccfg)
        t_json = _charptr()
        t_exitangle_nbins = _uint()
        t_exitangle_ct = _dblp()
        t_exitangle_errsq = _dblp()
        _raw_runmmcsim_stdengine( nthreads_,
                                  tally_detail_lvl_,
                                  mat_cfgstr_,
                                  mmc_geomcfg_,
                                  mmc_srccfg_,
                                  t_json,
                                  t_exitangle_nbins,
                                  ctypes.byref(t_exitangle_ct),
                                  ctypes.byref(t_exitangle_errsq) )
        if t_json is not None:
            t_json_cstr = ctypes.cast(t_json,_cstr).value
            if t_json_cstr is not None:
                _ = _cstr2str(t_json_cstr)
                _raw_deallocstr(t_json)
                t_json = _
            else:
                t_json = None

        tally_exitangle_contents = _cptr_to_nparray( t_exitangle_ct, t_exitangle_nbins )
        tally_exitangle_errsq = _cptr_to_nparray( t_exitangle_errsq, t_exitangle_nbins )
        return tally_exitangle_contents,tally_exitangle_errsq,t_json

    functions['runmmcsim_stdengine']=nc_runmmcsim_stdengine


    return functions

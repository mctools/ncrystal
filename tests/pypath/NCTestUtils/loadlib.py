
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2024 NCrystal developers                                   ##
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

"""

   A utility Lib class which can be used to access compiled test libraries, and
   the extern "C" functions defined within. If the library defines a
   nctest_ctypes_dictionary() function, it can even be loaded directly and the
   returned dictionary string will be inspected to set up functions
   automatically.

   Example of such a dictionary function:

   extern "C" const char * nctest_ctypes_dictionary()
   {
     return
       "int nctest_file_exists( const char * );"
       "const char * nctest_ncgetcwd();"
     ;
   }

"""

import ctypes

class Lib:
    def __init__( self, test_shlib_name ):
        self.__name = _normalise_testlib_name( test_shlib_name )
        self.__lib = _ctypes_load_testlib( test_shlib_name )
        if not hasattr(self.__lib,'nctest_ctypes_dictionary'):
            print("Warning: No nctest_ctypes_dictionary symbol"
                  " in testlib %s"%self.__name)
        else:
            dictfct = _ctypes_create_fct( self.__lib,
                                          'nctest_ctypes_dictionary',
                                          ctypes.c_char_p )
            dictstr = dictfct()
            for e in dictfct().split(';'):
                e=e.strip()
                if e:
                    self.add_signature(e)
            del dictfct

    def add_signature( self, fct_signature ):
        return self.add(fct_signature,
                        '__include_fct_name__')

    def add( self, fctname, restype, *argtypes ):
        if restype=='__include_fct_name__':
            ( fctname,
              restype,
              argtypes) = _decode_signature_str(fctname,
                                                include_fct_name = True)
        elif len(argtypes)==0 and isinstance(restype,str) and '(' in restype:
            ( restype,
              argtypes ) = _decode_signature_str(restype,
                                                 include_fct_name = False)
        fct = _ctypes_create_fct( self.__lib,
                                  fctname,
                                  restype,
                                  *argtypes )
        #Fixme test: fct.__name__ = fctname
        assert not hasattr(self,fctname),f'Fct {repr(fctname)} already added!'
        setattr(self,fctname,fct)

def _decode_type_str( s ):
    s=' '.join(s.replace('*',' * ').strip().split())
    m = { 'const char *':ctypes.c_char_p,
          'void' : 'void',
          'int':ctypes.c_int,
          'uint':ctypes.c_uint,
          'double':ctypes.c_double }
    return m.get(s,s)

def _decode_signature_str( signature, include_fct_name ):
    signature=signature.strip()
    assert signature.count('(')==1
    assert signature.count(')')==1
    assert signature.index(')')+1==len(signature)
    r,args = signature[:-1].split('(',2)
    args = list( a for a in args.split(',') ) if args.strip() else []

    if include_fct_name:
        r = r.split()
        assert len(r) >= 2
        fctname = r[-1]
        r = ' '.join(r[:-1])
    else:
        fctname = None
    r = _decode_type_str(r)
    args = tuple( _decode_type_str(a) for a in args )
    if fctname is None:
        return r,args
    else:
        return fctname,r,args

def _normalise_testlib_name(name):
    return name if name.startswith('TestLib_') else f'TestLib_{name}'
def _find_testlib(name):
    tln = _normalise_testlib_name(name)
    import pathlib
    import os
    locdir = pathlib.Path(os.environ.get('MCTOOLS_TESTMODULES_LOCDIR'))
    assert locdir.is_dir()
    f = locdir / f'module_loc_{tln}.txt'
    lib = pathlib.Path(f.read_text().strip())
    assert lib.is_file()
    return lib

def _load_lib_with_ctypes( path ):
    try:
        lib = ctypes.CDLL(path)
    except TypeError:
        lib = None

    if lib is None:
        #For some reason, on windows we get a TypeError and must pass a string
        #rather than a pathlib object:
        lib = ctypes.CDLL(str(path))
    return lib

def _ctypes_load_testlib( test_shlib_name ):
    libpath = _find_testlib(test_shlib_name)
    return _load_lib_with_ctypes(libpath)

def _ctypes_create_fct( lib, fctname, restype, *argtypes ):
    from NCrystal._chooks import _str2cstr, _cstr2str

    def resolve_type( tpe ):
        if tpe is None or tpe=='void':
            return None
        return ( getattr(ctypes,tpe)
                 if isinstance(tpe,str) and hasattr(ctypes,tpe)
                 else tpe )
    assert hasattr(lib,fctname), f"Missing symbol: {fctname}"
    rawfct = getattr(lib,fctname)
    rawfct.restype = resolve_type(restype)
    argtypes = tuple( resolve_type(a) for a in argtypes )
    rawfct.argtypes = argtypes
    def fct( *args ):
        if len(args) != len(argtypes):
            raise RuntimeError(f"{fctname}(..) takes {len(argtypes)} "
                               f"args ({len(args)} provided)")
        al = []
        for a,at in zip(args,argtypes):
            if at == ctypes.c_char_p:
                al.append( _str2cstr(a) )
            else:
                al.append( a )
        rv = rawfct( *al )
        if restype == ctypes.c_char_p:
            rv = _cstr2str( rv )
        return rv
    #fct.__name__== fixme
    return fct

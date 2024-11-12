#ifndef NCrystal_CFileUtils_hh
#define NCrystal_CFileUtils_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
//                                                                            //
//  Licensed under the Apache License, Version 2.0 (the "License");           //
//  you may not use this file except in compliance with the License.          //
//  You may obtain a copy of the License at                                   //
//                                                                            //
//      http://www.apache.org/licenses/LICENSE-2.0                            //
//                                                                            //
//  Unless required by applicable law or agreed to in writing, software       //
//  distributed under the License is distributed on an "AS IS" BASIS,         //
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  //
//  See the License for the specific language governing permissions and       //
//  limitations under the License.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// This file provides a C/C++ compatible backend for cross-platform filesystem
// utilities. It is NOT meant to be included directly (other NCrystal code
// should use the NCFileUtils.hh header instead!), but is made available here
// for easier testing, etc.
//
// This file (and its corresponding implementation file) are deliberately
// written as to be possible to use in both C and C++ projects (the latter
// should set MCFILEUTILS_CPPNAMESPACE to an appropriate namespace (handled
// automatically for NCrystal below). It is also used to implement the
// ncrystal-config command (written in C for maximal robustness).

#if defined(NCrystal_EXPORTS) || defined(NCRYSTAL_PRETEND_EXPORTS)
#  include "NCrystal/ncapi.h"
#  define MCFILEUTILS_CPPNAMESPACE NCRYSTAL_NAMESPACE
#endif

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////                                                                       ////
////                     C-compatible string utilities                     ////
////                                                                       ////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#ifdef MCFILEUTILS_CPPNAMESPACE
  //C++ mode, with namespace.
#  include <cstddef>
#  include <cstdio>
#  include <string>
#  include <fstream>
namespace MCFILEUTILS_CPPNAMESPACE {
#else
  //C mode
#  include <stddef.h>
#  include <stdio.h>
#endif

  //A string struct "mcu8str", for somewhat more convenient and safe string
  //handling in C. Note that client code can access the fields directly, but
  //should normally restrict themselves to accessing the three first fields. If
  //editing the string in place, one should be aware that "buflen" is the size
  //of the current buffer, and that the "size" field should ALWAYS indicate the
  //size which would be computed by strlen(c_str) (i.e. it is the index of the
  //first null-char). If c_str content has been modified (while not overstepping
  //the buffer size given in buflen of course!), mcu8str_update_size(..) can be
  //used to set size to strlen(c_str) (avoid doing this manually!).

  typedef struct {
    char * c_str; //Buffer, which is always a null-terminated string.
    unsigned int size;//strlen(c_str)
    unsigned int buflen;//Actual allocated size of buffer in c_str (including
                        //terminating null char).
    int owns_memory;//(for internal usage only).
  } mcu8str;

  void mcu8str_dealloc( mcu8str* );
  void mcu8str_update_size( mcu8str* );//sets size to strlen(c_str) after
                                       //editing c_str manually
  mcu8str mcu8str_create_empty();
  mcu8str mcu8str_create( size_t prealloc_size );//fixme std::size_t if c++
  mcu8str mcu8str_create_from_staticbuffer( char * buf, size_t buflen );
  mcu8str mcu8str_create_from_cstr( const char * );
  void mcu8str_append( mcu8str* str, const mcu8str * otherstr );
  void mcu8str_swap( mcu8str* str1, mcu8str* str2 );
  void mcu8str_assign( mcu8str* dest, const mcu8str* src );
  mcu8str mcu8str_copy( const mcu8str* );
  void mcu8str_reserve( mcu8str* str, size_t nsize );
  void mcu8str_clear( mcu8str* );
  void mcu8str_append_cstr( mcu8str* str, const char * );
  int mcu8str_is_ascii( const mcu8str* str );//ASCII (no chars have high bit set)
  void mcu8str_ensure_dynamic_buffer( mcu8str* );//stop referencing static buffer
  int mcu8str_contains( const mcu8str* str, char needle );
  void mcu8str_replace( mcu8str* str, char from, char to );

#ifdef MCFILEUTILS_CPPNAMESPACE
  std::string mcu8str_to_stdstring( const mcu8str* );
  mcu8str mcu8str_create_from_stdstring( const std::string& );
#endif

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////                                                                       ////
////                         File system utilities                         ////
////                                                                       ////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

  //Actual file utilities, meant to be portable:
  mcu8str mctools_determine_exe_self_path( int argc, char** argv );
  int mctools_file_exists_and_readable( const mcu8str* );
  mcu8str mctools_get_current_working_dir();
  mcu8str mctools_make_path_absolute( const mcu8str* );

  //Portable version of_fopen. As usual, this returns null ptr on failure and
  //otherwise a FILE handle:
#ifdef MCFILEUTILS_CPPNAMESPACE
  typedef std::FILE MCTOOLS_FILE_t;
#else
  typedef FILE MCTOOLS_FILE_t;
#endif
  MCTOOLS_FILE_t * mctools_fopen( const mcu8str* path, const char * mode );

#ifdef MCFILEUTILS_CPPNAMESPACE
}
#endif

#endif

#ifndef NCrystal_CFileUtils_hh
#define NCrystal_CFileUtils_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
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
// utilities. It is meant as a private backend and NOT to be used directly in
// user code. It is envisioned to be shared between NCrystal and MCPL projects.
//
// This file (and its corresponding implementation file) are deliberately
// written as to be possible to use in both C and C++ projects (the latter
// should set MCFILEUTILS_CPPNAMESPACE to an appropriate namespace.

//Special lines supporting usage in NCrystal:
#if !defined(MCFILEUTILS_CPPNAMESPACE) && defined(__cplusplus)
#  if defined(NCrystal_EXPORTS) || defined(NCRYSTAL_PRETEND_EXPORTS)
#    include "NCrystal/ncapi.h"
#    define MCFILEUTILS_CPPNAMESPACE NCRYSTAL_NAMESPACE
#  endif
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
  using mcu8str_size_t = std::size_t;
#else
  //C mode
#  include <stddef.h>
#  include <stdio.h>
  typedef size_t mcu8str_size_t;
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
  mcu8str mcu8str_create_empty(void);//empty path (safe to not deallocate)
  mcu8str mcu8str_create( mcu8str_size_t prealloc_size );
  mcu8str mcu8str_create_from_staticbuffer( char * buf, mcu8str_size_t buflen );
  mcu8str mcu8str_create_from_cstr( const char * );
  mcu8str mcu8str_copy( const mcu8str* );

  //Views (input strings should outlive resulting objects). Note they are
  //technically const only in C++, but should always be treated as such:
#ifdef MCFILEUTILS_CPPNAMESPACE
  const
#endif
  mcu8str mcu8str_view_cstr( const char * );
#ifdef MCFILEUTILS_CPPNAMESPACE
  const
#endif
  mcu8str mcu8str_view_str( const mcu8str * );

  void mcu8str_ensure_dynamic_buffer( mcu8str* );//replace any static with
                                                 //dynamic buffers
  void mcu8str_append( mcu8str* str, const mcu8str * otherstr );
  void mcu8str_swap( mcu8str* str1, mcu8str* str2 );
  void mcu8str_assign( mcu8str* dest, const mcu8str* src );
  void mcu8str_reserve( mcu8str* str, mcu8str_size_t nsize );
  void mcu8str_clear( mcu8str* );
  void mcu8str_append_cstr( mcu8str*, const char * );
  void mcu8str_replace( mcu8str*, char from, char to );
  int mcu8str_contains( const mcu8str*, char needle );
  int mcu8str_is_ascii( const mcu8str* );//ASCII (no chars have high bit set)
  int mcu8str_equal( const mcu8str*, const mcu8str* );

  //Convenience macro to create an mcu8str object viewing a string literal:
#define MCU8STR_FROMLITERAL( varame, str ) \
  mcu8str varname; { \
    varname.c_str = str; \
    varname.buflen = sizeof(str); \
    varname.size = sizeof(str)-1; \
    varname.owns_memory = 0; }

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

  //TODO: Make clear (also in naming) which functions are purely string based
  //and which are operating on the filesystem (e.g. a path-normalisation could
  //be purely string based, or it could actually resolve symlinks etc.)

  //Actual file utilities, meant to be portable:
  int mctools_exists( const mcu8str* );//path exists
  int mctools_is_dir( const mcu8str* );//path exists and is directory
  int mctools_is_file( const mcu8str* );//path exists and is not a directory
  mcu8str mctools_get_current_working_dir(void);//like unix getcwd

  //Check if two files exist and are actually the same file (i.e. same inode on
  //same device) in the system. Returns false for directories:
  int mctools_is_same_file( const mcu8str*, const mcu8str* );

  //For usage in main fcts to get path to executable self:
  mcu8str mctools_determine_exe_self_path( int argc, char** argv );


  //Try to normalise path representation by turning it into an absolute path and
  //bringing it into a normalised form. If it is a symlink, it might be resolved
  //as well if possible (depending on the platform). Depending on the platform,
  //this might also result in an error (i.e. return an empty path), for instance
  //if the path is non-existing or too long. The function will also replace any
  //path separators with platform-specific separators (if generic forward
  //slashes are desired, use a subsequent call to mctools_pathseps_generic):
  mcu8str mctools_real_path( const mcu8str* );

  //Bring to portable form, i.e. any forward or backwards slashes become forward
  //slashes. Also upper-cases any Windows drive letter:
  void mctools_pathseps_generic( mcu8str* );

  //Bring to platform specific form, i.e. any forward or backwards slashes
  //become backwards slashes if current platform is Windows, otherwise they
  //become forward slashes. Also upper-cases any Windows drive letter:
  void mctools_pathseps_platform( mcu8str* );

  //Check if path is absolute or relative:
  int mctools_path_is_absolute( const mcu8str* );
  int mctools_path_is_relative( const mcu8str* );

  //Get (uppercased) Windows drive letter from path (e.g. 'C' if path is
  //"c:\bla.txt"). Returns 0 if absent.
  char mctools_drive_letter( const mcu8str* path );

  //Various parts:
  mcu8str mctools_basename( const mcu8str* path );
  mcu8str mctools_dirname( const mcu8str* path );
  mcu8str mctools_fileextension( const mcu8str* path );

  //View versions working without any malloc, but results can only be used while
  //arguments are alive and unmodified
  const char * mctools_basename_view( const mcu8str* );
  const char * mctools_fileextension_view( const mcu8str* );

  //Join paths (with native path separator).
  mcu8str mctools_path_join( const mcu8str*, const mcu8str* );

  //Turn path into absolute path (does not resolve symlinks):
  mcu8str mctools_absolute_path( const mcu8str* );

  //On unix this expands a leading '~/' to the users home directory (only if the
  //HOME env var is set), and on Windows this expands paths with short 8.3-form
  //filenames to their long form. It also affects slashes and drive letters like
  //mctools_pathseps_platform:
  mcu8str mctools_expand_path( const mcu8str* );

  //Portable version of_fopen. As usual, this returns null ptr on failure and
  //otherwise a FILE handle:
#ifdef MCFILEUTILS_CPPNAMESPACE
  typedef std::FILE MCTOOLS_FILE_t;
#else
  typedef FILE MCTOOLS_FILE_t;
#endif
  MCTOOLS_FILE_t * mctools_fopen( const mcu8str* path, const char * mode );

#ifdef _WIN32
  //Get resolved wpath. Returns null pointer in case of problems. Caller must
  //free(..) the returned value.
  wchar_t* mctools_path2wpath( const mcu8str* );
  //Generic utf-16 wchar_t* to utf-8 char* string conversion:
  mcu8str mctool_wcharstr_to_u8str( const wchar_t* );
#endif

#ifdef MCFILEUTILS_CPPNAMESPACE
}
#endif

#endif

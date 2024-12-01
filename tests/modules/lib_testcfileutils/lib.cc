
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

#ifndef NCRYSTAL_PRETEND_EXPORTS
// Let NCCFileUtils.hh know it is compiled into NCrystal, not some other
// project.
#  define NCRYSTAL_PRETEND_EXPORTS
#endif
#include "NCrystal/internal/utils/NCCFileUtils.hh"

#include "NCTestUtils/NCTestModUtils.hh"
#include <cstring>

NCTEST_CTYPE_DICTIONARY
{
  return
    "const char * nctest_get_current_working_dir();"
    "const char * nctest_real_path( const char * );"
    "const char * nctest_absolute_path( const char * );"
    "const char * nctest_fopen_and_read_text( const char * );"
    "const char * nctest_pathseps_platform( const char* );"
    "const char * nctest_pathseps_generic( const char* );"
    "int nctest_path_is_absolute( const char* );"
    "int nctest_path_is_relative( const char* );"
    "const char * nctest_drive_letter( const char* );"
    "const char * nctest_basename( const char* );"
    "const char * nctest_basename_view( const char* );"
    "const char * nctest_dirname( const char* );"
    "const char * nctest_fileextension( const char* );"
    "const char * nctest_fileextension_view( const char* );"
    "const char * nctest_path_join( const char*, const char* );"
    "int nctest_is_same_file( const char *, const char * );"
    "int nctest_is_dir( const char* );"
    "int nctest_is_file( const char* );"
    "int nctest_exists( const char* );"
    "const char * nctest_expand_path( const char* );"
    ;
}

namespace {
  template<int BUFSIZE = 65536>
  class NCTestCharBuf {
    char m_buf[BUFSIZE+1];
  public:
    const char * consume( NC::mcu8str* str )
    {
      std::size_t size = str->size;
      std::size_t calc_size = (std::size_t)strlen(str->c_str);
      if ( size != calc_size ) {
        printf("Inconsistent size of mcu8str (size=%i but strlen=%i): \"%s\"\n",
               (int)size, int(calc_size), str->c_str );
        return "ERROR: inconsistent string size in mc8ustr!!";
      }
      assert(size < BUFSIZE);
      m_buf[0] = '\0';
      std::strncat(m_buf,str->c_str,BUFSIZE);
      NC::mcu8str_dealloc(str);
      return &m_buf[0];
    }
  };
}

NCTEST_CTYPES const char * nctest_expand_path( const char* path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  auto res =  NC::mctools_expand_path( &mpath );
  static NCTestCharBuf<> buf;
  return buf.consume(&res);
}

NCTEST_CTYPES const char * nctest_pathseps_platform( const char* path )
{
  auto res = NC::mcu8str_create_from_cstr(path);
  NC::mctools_pathseps_platform( &res );
  static NCTestCharBuf<> buf;
  return buf.consume(&res);
}

NCTEST_CTYPES const char * nctest_pathseps_generic( const char* path )
{
  auto res = NC::mcu8str_create_from_cstr(path);
  NC::mctools_pathseps_generic( &res );
  static NCTestCharBuf<> buf;
  return buf.consume(&res);
}

NCTEST_CTYPES int nctest_path_is_absolute( const char* path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  return mctools_path_is_absolute( &mpath );
}

NCTEST_CTYPES int nctest_path_is_relative( const char* path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  return mctools_path_is_relative( &mpath );
}

NCTEST_CTYPES const char * nctest_drive_letter( const char* path )
{
  auto mpath = NC::mcu8str_view_cstr( path );
  static char buf[2];
  buf[0] = NC::mctools_drive_letter( &mpath );
  buf[1] = '\0';
  return buf;
}

NCTEST_CTYPES const char * nctest_basename( const char* path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  auto res =  NC::mctools_basename( &mpath );
  static NCTestCharBuf<> buf;
  return buf.consume(&res);
}

NCTEST_CTYPES const char * nctest_basename_view( const char* path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  const char * view =  NC::mctools_basename_view( &mpath );
  auto res = NC::mcu8str_create_from_cstr(view);
  static NCTestCharBuf<> buf;
  return buf.consume(&res);
}

NCTEST_CTYPES const char * nctest_fileextension( const char* path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  auto res =  NC::mctools_fileextension( &mpath );
  static NCTestCharBuf<> buf;
  return buf.consume(&res);
}

NCTEST_CTYPES const char * nctest_fileextension_view( const char* path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  const char * view =  NC::mctools_fileextension_view( &mpath );
  auto res = NC::mcu8str_create_from_cstr(view);
  static NCTestCharBuf<> buf;
  return buf.consume(&res);
}

NCTEST_CTYPES const char * nctest_dirname( const char* path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  auto res =  NC::mctools_dirname( &mpath );
  static NCTestCharBuf<> buf;
  return buf.consume(&res);
}

NCTEST_CTYPES const char * nctest_path_join( const char* p1, const char* p2 )
{
  auto mp1 = NC::mcu8str_view_cstr(p1);
  auto mp2 = NC::mcu8str_view_cstr(p2);
  auto res = NC::mctools_path_join( &mp1, &mp2 );
  static NCTestCharBuf<> buf;
  return buf.consume(&res);
}

NCTEST_CTYPES const char * nctest_fopen_and_read_text( const char * path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  auto fh = NC::mctools_fopen( &mpath, "r" );
  if (!fh)
    return "<<READ_FAILURE>>";
  static char buf[65536];//ok for these testing purposes
  auto nread = std::fread( &buf[0], 1, sizeof(buf), fh );
  if ( std::feof(fh) ) {
    //reached end of file
    std::fclose(fh);
    if ( nread < sizeof(buf) )
      buf[nread] = '\0';
    if ( std::strlen(buf) != nread )
      return "<<READ_FAILURE>>";
    return buf;
  }
  fclose(fh);
  return "<<READ_FAILURE>>";
}

NCTEST_CTYPES const char * nctest_get_current_working_dir()
{
  auto res = NC::mctools_get_current_working_dir();
  static NCTestCharBuf<> buf;
  return buf.consume(&res);
}

NCTEST_CTYPES const char * nctest_real_path( const char * path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  auto res = NC::mctools_real_path( &mpath );
  static NCTestCharBuf<> buf;
  return buf.consume(&res);
}

NCTEST_CTYPES const char* nctest_absolute_path( const char* path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  auto res = NC::mctools_absolute_path( &mpath );
  static NCTestCharBuf<> buf;
  return buf.consume(&res);
}

NCTEST_CTYPES int nctest_is_same_file( const char * p1, const char * p2 )
{
  auto mp1 = NC::mcu8str_view_cstr(p1);
  auto mp2 = NC::mcu8str_view_cstr(p2);
  return NC::mctools_is_same_file( &mp1, &mp2 );
}

NCTEST_CTYPES int nctest_is_dir( const char* path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  return NC::mctools_is_dir( &mpath );
}

NCTEST_CTYPES int nctest_is_file( const char* path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  return NC::mctools_is_file( &mpath );
}

NCTEST_CTYPES int nctest_exists( const char* path )
{
  auto mpath = NC::mcu8str_view_cstr(path);
  return NC::mctools_exists( &mpath );
}

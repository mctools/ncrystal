
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

#include "NCTestUtils/NCTestModUtils.hh"
#include "NCrystal/internal/utils/NCFileUtils.hh"
#include <cstring>

NCTEST_CTYPE_DICTIONARY
{
  return
    "int nctest_file_exists( const char * );"
    "const char * nctest_ncgetcwd();"
    "const char * nctest_readEntireFileToString( const char * );"
    "const char * nctest_ncglob( const char * );"
    "const char * nctest_dirname( const char *);"
    "const char * nctest_basename( const char *);"
    "const char * nctest_getfileext( const char *);"
    "const char * nctest_normalise( const char *);"
    "int nctest_path_is_absolute( const char* );"
    ;
}

//All I/O from these functions are in UTF-8.





NCTEST_CTYPES const char * nctest_readEntireFileToString( const char * path )
{
  //Returns "<<NULL>>" in case readEntireFileToString did not return a value.

  //For testing, we can get away with just using a large static buffer:

  static char buf[10485760];//10MB, hopefully enough for even the largest file
                            //in our tests
  try {
    nc_assert_always( path );
    auto opt_content = NC::readEntireFileToString( path );
    const char * content = ( opt_content.has_value()
                             ? opt_content.value().c_str()
                             : "<<NULL>>" );
    buf[0] = '\0';
    std::strncat(buf,content,sizeof(buf)-1);
  } NCCATCH;
  return &buf[0];
}


NCTEST_CTYPES int nctest_file_exists( const char * f )
{
  int res = -1;
  try {
    nc_assert_always( f );
    res = NC::file_exists( f ) ? 1 : 0;
  } NCCATCH;
  return res;
}

NCTEST_CTYPES const char * nctest_ncgetcwd()
{
  //For testing, we can get away with just using a large static buffer:
  static char buf[16384];
  try {
    std::string thecwd = NC::ncgetcwd();
    buf[0] = '\0';
    std::strncat(buf,thecwd.c_str(),sizeof(buf)-1);
  } NCCATCH;
  return &buf[0];
}


NCTEST_CTYPES const char * nctest_ncglob(const char * pattern)
{
  //For testing, we can get away with just using a large static buffer:
  static char buf[10485760];//10MB, plenty!
  try {
    std::ostringstream ss;
    bool first = true;
    for (auto&e : NC::ncglob(pattern) ) {
      if ( !first )
        ss << "<<@>>";//separator
      first = false;
      ss << e;
    }
    buf[0] = '\0';
    std::strncat(buf,ss.str().c_str(),sizeof(buf)-1);
  } NCCATCH;
  return &buf[0];
}


NCTEST_CTYPES const char * nctest_dirname( const char * path )
{
  static char buf[16384];
  try {
    auto res = NC::dirname(path);
    buf[0] = '\0';
    std::strncat(buf,res.c_str(),sizeof(buf)-1);
  } NCCATCH;
  return &buf[0];
}

NCTEST_CTYPES const char * nctest_basename( const char * path )
{
  static char buf[16384];
  try {
    auto res = NC::basename(path);
    buf[0] = '\0';
    std::strncat(buf,res.c_str(),sizeof(buf)-1);
  } NCCATCH;
  return &buf[0];
}

NCTEST_CTYPES const char * nctest_getfileext( const char * path )
{
  static char buf[16384];
  try {
    auto res = NC::getfileext(path);
    buf[0] = '\0';
    std::strncat(buf,res.c_str(),sizeof(buf)-1);
  } NCCATCH;
  return &buf[0];
}

NCTEST_CTYPES const char * nctest_normalise( const char * path )
{
  static char buf[16384];
  try {
    auto res = NC::normalise(path);
    buf[0] = '\0';
    std::strncat(buf,res.c_str(),sizeof(buf)-1);
  } NCCATCH;
  return &buf[0];
}

NCTEST_CTYPES int nctest_path_is_absolute( const char* path )
{
  bool res( false );
  try {
    res = NC::path_is_absolute(path);
  } NCCATCH;
  return ( res ? 1 : 0 );
}

#include <cstring>
#include "NCrystal/NCDefs.hh"
#include "NCrystal/internal/NCFileUtils.hh"
#include <iostream>
namespace NC = NCrystal;

extern "C" const char * nctest_ctypes_dictionary()
{
  return
    "int nctest_file_exists( const char * );"
    "const char * nctest_ncgetcwd();"
    ;
}

extern "C" int nctest_file_exists( const char * f )
{
  nc_assert_always(f);
  std::cout<<"TKTEST "<<f<<" -> "<<NC::file_exists( f )<<std::endl;
  return NC::file_exists( f ) ? 1 : 0;
}

extern "C" const char * nctest_ncgetcwd()
{
  //For testing, we can get away with just using a large static buffer:
  std::string thecwd = NC::ncgetcwd();
  static char buf[16384];
  buf[0] = '\0';
  std::strncat(buf,thecwd.c_str(),sizeof(buf)-1);
  return &buf[0];
}


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
#include "NCrystal/internal/NCCFileUtils.hh"

#include <stdio.h>

namespace NC=NCrystal;

#define ncrystal_xstr(a) ncrystal_str(a)
#define ncrystal_str(a) #a

int main(int argc, char** argv) {

  NC::mcu8str selfpath = NC::mctools_determine_exe_self_path( argc, argv );
  printf("Self path:          \"%s\"\n",selfpath.c_str);

  if ( !NC::mctools_file_exists_and_readable( &selfpath ) ) {
    printf("ERROR: self path is not a (readable) file!\n");
    return 1;
  }

  //Test versus the expected path. We get the expected path from CMake, which
  //has set the MCTOOLS_TESTAPP_FILE define to a string with the path:

  NC::mcu8str expected_selfpath =
    NC::mcu8str_create_from_cstr( MCTOOLS_TESTAPP_FILE );

  //Note that CMake might have used forward instead of backwards slashes on
  //windows (e.g. "D:/a/some/where/selfpathc.exe" instead of
  //"D:\a\some\where\selfpathc.exe") so we should be sure to not get a spurious
  //failure here, by normalising:



  printf("Expected self path: \"%s\"\n",expected_selfpath.c_str);
  if ( !NC::mctools_file_exists_and_readable( &expected_selfpath ) ) {
    printf("ERROR: expected self path is not a (readable) file!\n");
    return 1;
  }

  mctools_pathseps_generic(&selfpath);
  mctools_pathseps_generic(&expected_selfpath);

  if ( !NC::mcu8str_equal(&selfpath,&expected_selfpath) ) {
    printf("ERROR: Mismatch!!\n");
    return 1;
  }
  NC::mcu8str_dealloc( &selfpath );
  NC::mcu8str_dealloc( &expected_selfpath );
  return 0;
}

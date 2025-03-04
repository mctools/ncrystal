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

#ifndef NCRYSTAL_PRETEND_EXPORTS
#  define NCRYSTAL_PRETEND_EXPORTS
#endif
#include "NCrystal/internal/utils/NCCFileUtils.hh"

#include <cstdio>
#ifdef NCRYSTAL_SIMPLEBUILD_DEVEL_MODE
#  include <cstdlib>//for getenv
#endif

namespace NC = NCrystal;

int main(int argc, char** argv) {

  NC::mcu8str selfpath = NC::mctools_determine_exe_self_path( argc, argv );
  printf("Self path:          \"%s\"\n",selfpath.c_str);

  if ( !NC::mctools_is_file( &selfpath ) ) {
    printf("ERROR: self path is not a file!\n");
    return 1;
  }

#ifndef NCRYSTAL_SIMPLEBUILD_DEVEL_MODE
  //CMake/CTest has provided MCTOOLS_TESTAPP_FILE with out location
  NC::mcu8str expected_selfpath
    = NC::mcu8str_create_from_cstr( MCTOOLS_TESTAPP_FILE );
#else
  //Expect in $SBLD_INSTALL_PREFIX/bin
  const char * env_sbldinstdir = getenv("SBLD_INSTALL_PREFIX");
  if (!env_sbldinstdir) {
    printf("ERROR: SBLD_INSTALL_PREFIX not set!\n");
    return 1;
  }
  NC::mcu8str expected_selfpath = NC::mcu8str_create_from_cstr( env_sbldinstdir );
  NC::mcu8str_append_cstr( &expected_selfpath, "/bin/sb_nctestapps_testselfpath" );
#endif

  //Note that CMake might have used forward instead of backwards slashes on
  //windows (e.g. "D:/a/some/where/selfpath.exe" instead of
  //"D:\a\some\where\selfpath.exe") so we should be sure to not get a spurious
  //failure here, by normalising:

  printf("Expected self path: \"%s\"\n",expected_selfpath.c_str);
  if ( !NC::mctools_is_file( &expected_selfpath ) ) {
    printf("ERROR: expected self path is not a file!\n");
    return 1;
  }

  NC::mctools_pathseps_generic(&selfpath);
  NC::mctools_pathseps_generic(&expected_selfpath);

  if ( !NC::mcu8str_equal(&selfpath,&expected_selfpath) ) {
    printf("ERROR: Mismatch!!\n");
    return 1;
  }
  NC::mcu8str_dealloc( &selfpath );
  NC::mcu8str_dealloc( &expected_selfpath );
  return 0;
}

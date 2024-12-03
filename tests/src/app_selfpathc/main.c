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

#include "NCrystal/ncapi.h"//For NCRYSTAL_SIMPLEBUILD_DEVEL_MODE define
#include "NCCFileUtils.h"

#include <stdio.h>
#ifdef NCRYSTAL_SIMPLEBUILD_DEVEL_MODE
#  include <stdlib.h>//for getenv
#endif

int main(int argc, char** argv) {

  mcu8str selfpath = mctools_determine_exe_self_path( argc, argv );
  printf("Self path:          \"%s\"\n",selfpath.c_str);

  if ( !mctools_is_file( &selfpath ) ) {
    printf("ERROR: self path is not a file!\n");
    return 1;
  }

#ifndef NCRYSTAL_SIMPLEBUILD_DEVEL_MODE
  //CMake/CTest has provided MCTOOLS_TESTAPP_FILE with out location
  mcu8str expected_selfpath
    = mcu8str_create_from_cstr( MCTOOLS_TESTAPP_FILE );
#else
  //Expect in $SBLD_INSTALL_PREFIX/bin
  const char * env_sbldinstdir = getenv("SBLD_INSTALL_PREFIX");
  if (!env_sbldinstdir) {
    printf("ERROR: SBLD_INSTALL_PREFIX not set!\n");
    return 1;
  }
  mcu8str expected_selfpath = mcu8str_create_from_cstr( env_sbldinstdir );
  mcu8str_append_cstr( &expected_selfpath, "/bin/sb_nctestapps_testselfpathc" );
#endif

  //Note that CMake might have used forward instead of backwards slashes on
  //windows (e.g. "D:/a/some/where/selfpathc.exe" instead of
  //"D:\a\some\where\selfpathc.exe") so we should be sure to not get a spurious
  //failure here, by normalising:

  printf("Expected self path: \"%s\"\n",expected_selfpath.c_str);
  if ( !mctools_is_file( &expected_selfpath ) ) {
    printf("ERROR: expected self path is not a file!\n");
    return 1;
  }

  mctools_pathseps_generic(&selfpath);
  mctools_pathseps_generic(&expected_selfpath);

  if ( !mcu8str_equal(&selfpath,&expected_selfpath) ) {
    printf("ERROR: Mismatch!!\n");
    return 1;
  }
  mcu8str_dealloc( &selfpath );
  mcu8str_dealloc( &expected_selfpath );
  return 0;
}

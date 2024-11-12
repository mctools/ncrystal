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
  printf("Self path:          %s\n",selfpath.c_str);

  NC::mcu8str expected_selfpath =
    NC::mcu8str_create_from_cstr( MCTOOLS_TESTAPP_FILE );
  printf("Expected self path: %s\n",expected_selfpath.c_str);
  if ( std::string(selfpath.c_str) != expected_selfpath.c_str ) {
    printf("ERROR: Mismatch!!\n");
    return 1;
  }
  NC::mcu8str_dealloc( &selfpath );
  NC::mcu8str_dealloc( &expected_selfpath );
  return 0;
}


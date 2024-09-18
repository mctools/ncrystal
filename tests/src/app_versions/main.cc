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

#include <string>
#include <sstream>
//Check that version macros in two headers are identical, and self-consistent.

#include "NCrystal/ncrystal.h"
namespace ncv {
  int ncrystal_version_major = NCRYSTAL_VERSION_MAJOR;
  int ncrystal_version_minor = NCRYSTAL_VERSION_MINOR;
  int ncrystal_version_patch = NCRYSTAL_VERSION_PATCH;
  int ncrystal_version = NCRYSTAL_VERSION;
  std::string ncrystal_version_str = NCRYSTAL_VERSION_STR;
}
#undef NCRYSTAL_VERSION_MAJOR
#undef NCRYSTAL_VERSION_MINOR
#undef NCRYSTAL_VERSION_PATCH
#undef NCRYSTAL_VERSION
#undef NCRYSTAL_VERSION_STR
#include "NCrystal/NCVersion.hh"

int main(int,char**) {
  if ( ncv::ncrystal_version_major != NCRYSTAL_VERSION_MAJOR) return 1;
  if ( ncv::ncrystal_version_minor != NCRYSTAL_VERSION_MINOR) return 1;
  if ( ncv::ncrystal_version_patch != NCRYSTAL_VERSION_PATCH) return 1;
  if ( ncv::ncrystal_version != NCRYSTAL_VERSION) return 1;
  if ( ncv::ncrystal_version_str != NCRYSTAL_VERSION_STR) return 1;
  if ( (NCRYSTAL_VERSION_MAJOR * 1000000
        + NCRYSTAL_VERSION_MINOR * 1000
        + NCRYSTAL_VERSION_PATCH ) != NCRYSTAL_VERSION)
    return 1;

  std::stringstream s;
  s << NCRYSTAL_VERSION_MAJOR << "." << NCRYSTAL_VERSION_MINOR << "." << NCRYSTAL_VERSION_PATCH;
  if (s.str() != NCRYSTAL_VERSION_STR)
    return 1;
  return 0;
}

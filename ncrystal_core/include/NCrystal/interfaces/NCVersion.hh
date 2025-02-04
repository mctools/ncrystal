#ifndef NCrystal_Version_hh
#define NCrystal_Version_hh

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

//NB Several VERSION defines are set in ncapi.h:
// #define NCRYSTAL_VERSION_MAJOR (integer)
// #define NCRYSTAL_VERSION_MINOR (integer)
// #define NCRYSTAL_VERSION_PATCH (integer)
// #define NCRYSTAL_VERSION  (integer = 1000000*MAJOR+1000*MINOR+PATCH)
// #define NCRYSTAL_VERSION_STR (string ="major.minor.patch")

#include "NCrystal/ncapi.h"
#include <stdexcept>


#ifndef NCRYSTAL_VERSION
#  error "inconsistency detected"
#endif

namespace NCRYSTAL_NAMESPACE {

  //Function which returns NCRYSTAL_VERSION. If it does not, it indicates symbol
  //clashes from multiple installations of NCrystal.
  NCRYSTAL_API int getVersion();

  //If compiled with NCRYSTAL_NAMESPACE_PROTECTION, return the namespace here
  //(will be an empty string in default installations):
  NCRYSTAL_API const char* getBuildNameSpace();


  //Call in client code to detect broken installations where there is a mismatch
  //in versions in NCrystal headers included and NCrystal library loaded. Raise
  //generic exception rather than NCrystal exception in this case, since we
  //can't trust custom exceptions to work in broken environments:
  inline void libClashDetect() {
    if ( getVersion() != NCRYSTAL_VERSION )
      throw std::runtime_error( "Broken NCrystal installation detected "
                                "(the NCrystal header files included when building your code "
                                "are incompatible with the linked NCrystal library). "
                                "Most likely you have multiple conflicting NCrystal installations." );
  }

}

#endif

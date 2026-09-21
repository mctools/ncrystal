////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

////////////////////////////////////////////////////////////////////////////////
// Verifies that std::fma() gives the exact, correctly-rounded-fused result   //
// on this platform.                                                          //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/core/NCDefs.hh"
#include <iostream>
#include <cmath>

int main()
{
  //With a=1+2^-27, the exact product a*a = 1 + 2^-26 + 2^-54 needs one bit
  //more than a double has, so with c=-(1+2^-26): fma(a,a,c) is 2^-54 if
  //(and only if) the product is not rounded before the addition. Taken from
  //a volatile so it cannot be constant-folded:
  volatile double va = 1.0 + std::ldexp(1.0,-27);
  const double a = va;
  const double c = -( 1.0 + std::ldexp(1.0,-26) );

  if ( std::fma( a, a, c ) != std::ldexp(1.0,-54) ) {
    std::cout << "ERROR: std::fma is not exact on this platform" << std::endl;
    return 1;
  }
  std::cout << "std::fma is exact" << std::endl;
  return 0;
}

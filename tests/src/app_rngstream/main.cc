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

#include "NCrystal/interfaces/NCRNG.hh"
#include <iostream>
namespace NC = NCrystal;

void testrng( NC::RNG& rng )
{
  for ( auto i : NC::ncrange( 100000000 ) ) {
    auto val = rng.generate64RndmBits();
    auto n = i + 1;
    if ( n <= 10
         || ( n <= 1000000 && n % 10000 == 0 )
         || ( n <= 100000000 && n % 1000000 == 0 )
         || ( n % 10000000 == 0 ) )
      {
        std::cout << "Gen 64bits #"<<n<<" : "<< val<<std::endl;
      }
  }
}

void testrngdbl( NC::RNG& rng )
{
  for ( auto i : NC::ncrange( 1000000 ) ) {
    auto val = rng.generate();
    auto n = i + 1;
    if ( n <= 10 || ( n <= 1000000 && n % 10000 == 0 ) ) {
      std::cout << "Gen Dbl #"<<n<<" : "<< NC::fmt(val)<<std::endl;
    }
  }
}

int main ()
{
  auto rng = NC::getRNG();
  testrng( rng );
  testrngdbl( rng );
  return 0;
}

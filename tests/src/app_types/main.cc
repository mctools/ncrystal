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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include <iostream>
namespace NC = NCrystal;

namespace {
  void fct1( NC::NeutronEnergy a, NC::NeutronWavelength b ) {
    std::cout<<"Fct1 got: "<<a<<" and "<<b<<std::endl;
  }
  void fct2( const NC::NeutronDirection& dir ) {
    std::cout<<"Fct2 got: "<<dir<<std::endl;
  }
}

int main() {

  static_assert(sizeof(NC::NeutronWavelength)==sizeof(double),"");
  static_assert(sizeof(NC::NeutronEnergy)==sizeof(double),"");
  static_assert(sizeof(NC::NeutronDirection)==3*sizeof(double),"");//not actually guaranteed in case of padding...

  std::cout << NC::NeutronEnergy(0.025)<<std::endl;
  std::cout << NC::NeutronWavelength(1.8)<<std::endl;
  std::cout << NC::NeutronWavelength(NC::NeutronEnergy(0.025))<<std::endl;
  std::cout << NC::NeutronEnergy(NC::NeutronWavelength(1.8))<<std::endl;
  std::cout << NC::NeutronEnergy(0.025).get()<<std::endl;
  //  fct1(0.1,0.3);//does not compile (and it should not!)
  fct1(NC::NeutronWavelength(1.8),NC::NeutronEnergy(0.025));
  fct1(NC::NeutronEnergy(0.025),NC::NeutronEnergy(0.025));
  fct1(NC::NeutronWavelength(1.8),NC::NeutronWavelength(1.8));
  //fct1(NC::NeutronEnergy(0.025),1.0);

  NC::NeutronEnergy ekin(0.025);
  std::cout<<ekin<<std::endl;
  std::cout<<NC::NeutronWavelength(ekin)<<std::endl;
  std::cout<<(NC::NeutronWavelength)ekin<<std::endl;
  std::cout<<(const NC::NeutronWavelength)ekin<<std::endl;

  NC::NeutronWavelength wl(1.8);
  std::cout<<wl<<std::endl;
  std::cout<<NC::NeutronEnergy(wl)<<std::endl;
  std::cout<<(NC::NeutronEnergy)wl<<std::endl;
  std::cout<<(const NC::NeutronEnergy)wl<<std::endl;

  std::cout << NC::CosineScatAngle(0.99) << std::endl;
  std::cout << NC::CrossSect(5.55) << std::endl;
  //does not compile (good!):  std::cout << NC::CrossSect(NC::CosineScatAngle(0.99)) << std::endl;

  NC::NeutronDirection dir(0,0,1);
  std::cout << dir << std::endl;
  std::array<double,3> arr = { 0,1,0 };
  NC::NeutronDirection dir2(arr);
  std::cout << dir2 << std::endl;
  fct2(NC::NeutronDirection{1,0,0});
  fct2(NC::NeutronDirection{arr});
  double carr[3] = { 1, 0, 0 };
  fct2(NC::NeutronDirection{carr});

  NC::Vector vv( 1, 0, 0 );
  NC::NeutronDirection dir3(vv.as<NC::NeutronDirection>());
  NC::Vector v2(dir.as<NC::Vector>());
  (void)v2;//for clang
  (void)dir3;//for clang

  fct2(NC::NeutronDirection(std::array<double,3>()));

}

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

#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/internal/cfgutils/NCExtinctionCfg.hh"
#include <iostream>
namespace NC=NCrystal;

int main() {
  std::cout<<"TESTa"<<std::endl;
  NC::Cfg::ExtinctionCfg ecfg( NC::StrView("10e-6m") );//fixme why StrView??
  std::cout<<"BLA:";
  ecfg.stream(std::cout);
  std::cout<<":"<<std::endl;




// cfgutils/NCExtinctionCfg.hh

  auto cfg = NC::MatCfg("Al_sg225.ncmat;extn=10e-6m");
//   auto ecfg =
//  std::cout<<"TESTb"<< cfg.get_extn()<<std::endl;
  std::cout<<cfg<<std::endl;
   std::cout<<"TESTc"<<std::endl;
  return 0;
}

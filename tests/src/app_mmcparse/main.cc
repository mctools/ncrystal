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

#include "NCrystal/NCrystal.hh"
#include "NCrystal/internal/minimc/NCMMC_EngineOpts.hh"
#include <iostream>

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

void testcfg( const char * eopt_str )
{
  std::cout << std::endl;
  std::cout << "Test engine opt string: \""<<eopt_str<<'"' << std::endl;
  auto eopt = NCMMC::parseEngineOpts(eopt_str);
  std::cout << "    -> obj: "<<eopt<<std::endl;
  std::cout << "    -> json: ";
  NCMMC::engineOptsToJSON(std::cout,eopt);
  std::cout<<std::endl;
}

int main()
{
  //Loading a few simple materials in this non-Python test is useful for various
  //debugging scenarios:
  testcfg("");
  testcfg("nthreads=4");
  testcfg("tallybreakdown  =  0");
  testcfg("nthreads=0 ;   ignoremiss =1");
  testcfg("std;nthreads=0;ignoremiss=1");
  return 0;
}

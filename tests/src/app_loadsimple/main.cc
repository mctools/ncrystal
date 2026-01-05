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

#include "NCrystal/NCrystal.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include <iostream>

namespace NC = NCrystal;

void testcfg( const char * cfgstr )
{
  std::cout<<">>> Calling createInfo(\""<<cfgstr<<"\")"<<std::endl;
  auto info = NC::FactImpl::createInfo(cfgstr);
  std::cout<<">>> Calling dump(info)"<<std::endl;
  NC::dump( std::cout, info );
  std::cout<<">>> Dump complete"<<std::endl;

  std::cout<<">>> Calling createScatter(\""<<cfgstr<<"\")"<<std::endl;
  auto scat = NC::FactImpl::createScatter(cfgstr);
  (void)scat;

  std::cout<<">>> Calling createAbsorption(\""<<cfgstr<<"\")"<<std::endl;
  auto absn = NC::FactImpl::createAbsorption(cfgstr);
  (void)absn;
}

int main()
{
  //Loading a few simple materials in this non-Python test is useful for various
  //debugging scenarios:
  testcfg("stdlib::Al_sg225.ncmat");
  testcfg("stdlib::Al_sg225.ncmat;dcutoff=0.6;vdoslux=2");
  testcfg("stdlib::Ni_sg225.ncmat");
  testcfg("stdlib::Ni_sg225.ncmat;dcutoff=0.6;vdoslux=2");
  testcfg("Polyethylene_CH2.ncmat;density=0.95gcm3");
  testcfg("gasmix::air/-10C/0.8atm/0.30relhumidity");
  testcfg("Ge_sg227.ncmat;mos=20.0arcmin;dir1=@crys_hkl:5,1,1"
          "@lab:0,0,1;dir2=@crys_hkl:0,-1,1@lab:0,1,0");
  testcfg("phases<0.1*PbS_sg225_LeadSulfide.ncmat"
          "&0.9*Epoxy_Araldite506_C18H20O3.ncmat>;temp=296K");
  return 0;
}

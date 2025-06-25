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
#include "NCrystal/internal/cfgutils/NCCfgExtn.hh"
#include "NCrystal/internal/cfgutils/NCCfgVars.hh"
#include <iostream>

namespace NC=NCrystal;
namespace NCCE = NCrystal::Cfg::Extn;

void test_extncfg( const char * cfgstrval )
{
  const auto varname = "extn";
  std::cout <<"\nTESTING \""<<varname<<"="<<cfgstrval<<"\""<<std::endl;
  auto varid_extn = NC::Cfg::constexpr_varIdFromName( varname );
  NC::Cfg::CfgKeyValMap cfgdata = NCCE::decode_cfgstr( varid_extn, cfgstrval );
  std::cout <<"  Normalised: \"";
  NCCE::stream_to_cfgstr(std::cout,cfgdata);
  std::cout<<"\""<<std::endl;

  auto extncfg_base = NCCE::ExtnCfg_Base::decode( cfgdata );
  std::cout <<"  Base: \""<<extncfg_base<<"\""<<std::endl;
  nc_assert_always( extncfg_base.model == NCCE::Model::Sabine );
  if ( extncfg_base.model == NCCE::Model::Sabine ) {
    auto extncfg_sabine = NCCE::ExtnCfg_Sabine::decode( cfgdata );
    std::cout <<"  Sabine: \""<<extncfg_sabine<<"\""<<std::endl;
  } else {
    nc_assert_always(false);//new model, need to add tests.
  }

}



int main() {
  std::cout<<NC::MatCfg("Al_sg225.ncmat;extn=10e-6m")<<std::endl;
  test_extncfg("10e-6m");
  test_extncfg("1e-5m");
  test_extncfg("1");
  test_extncfg("0.01mm/1mm/2deg");
  test_extncfg("   0.01   mm  /   1  mm   / 2.0   deg  ");
  test_extncfg("   0.01   mm  /   1  mm   / .123456789123456789123456789   deg  ");
  test_extncfg("   0.01   mm  /   1  mm   / .123456789123456789123456789   rad  ");

  //fixme: many more tests!

  return 0;
}

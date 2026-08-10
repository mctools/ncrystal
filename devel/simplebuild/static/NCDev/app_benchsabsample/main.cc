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

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include <iostream>


namespace NC = NCrystal;

int main( int argc, char** argv )
{
  //Parse args:
  if ( argc!=4 ) {
    std::cout << "Please provide args: [cfgstr] [ekin_eV] [nsample]"
              << std::endl;
    return 1;
  }
  const std::string cfgstr(argv[1]);

  auto opt_ekin = NC::StrView(argv[2]).toDbl();
  nc_assert_always(opt_ekin.has_value());
  const NC::NeutronEnergy ekin{ NC::DoValidate, opt_ekin.value() };

  auto opt_nsample = NC::StrView(argv[3]).toUInt64();
  nc_assert_always(opt_nsample.has_value());
  const std::uint64_t nsample = opt_nsample.value();

  auto scatter = NC::FactImpl::createScatter(cfgstr);
  auto rng = NC::getRNG();

  std::cout<<"Sampling "<<nsample<<" times at ekin = "<<ekin<<std::endl;

  NC::CachePtr cache;
  for ( std::uint64_t i = 0; i < nsample; ++i )
    (void)scatter->sampleScatterIsotropic(cache,rng,ekin);
  return 0;
}

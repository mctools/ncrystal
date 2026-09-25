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
  if ( argc!=4 && argc!=5 ) {
    std::cout << "Please provide args: [cfgstr] [ekin_eV] [nsample] [ekin_hi_eV]"
              << std::endl;
    std::cout << "(ekin_hi_eV is optional: if given, a fresh energy is drawn"
                 " log-uniformly in [ekin_eV,ekin_hi_eV] for every sample,"
                 " rather than reusing the single ekin_eV value throughout)."
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

  NC::Optional<double> opt_ekin_hi;
  if ( argc==5 ) {
    auto v = NC::StrView(argv[4]).toDbl();
    nc_assert_always(v.has_value());
    nc_assert_always(v.value()>=ekin.dbl());
    opt_ekin_hi = v.value();
  }

  auto scatter = NC::FactImpl::createScatter(cfgstr);
  auto rng = NC::getRNG();

  if ( !opt_ekin_hi.has_value() ) {
    std::cout<<"Sampling "<<nsample<<" times at ekin = "<<ekin<<std::endl;
    NC::CachePtr cache;
    for ( std::uint64_t i = 0; i < nsample; ++i )
      (void)scatter->sampleScatterIsotropic(cache,rng,ekin);
    return 0;
  }

  //Fresh (log-uniform) random energy per call rather than reusing one fixed
  //ekin throughout: a single repeated ekin always lands in the same
  //cell/grid region, which is unrepresentative of an actual simulation (and
  //can make some costs -- e.g. anything that only runs once per distinct
  //region and is then effectively free on every subsequent identical call --
  //look artificially cheap in a profile).
  const double ekin_lo = ekin.dbl();
  const double ekin_hi = opt_ekin_hi.value();
  std::cout<<"Sampling "<<nsample<<" times at ekin = "
           <<ekin_lo<<" .. "<<ekin_hi<<" eV"<<std::endl;
  const double logratio = std::log(ekin_hi/ekin_lo);
  NC::CachePtr cache;
  for ( std::uint64_t i = 0; i < nsample; ++i ) {
    const NC::NeutronEnergy ekin_i{ NC::DoValidate,
      ekin_lo*std::exp(logratio*rng->generate()) };
    (void)scatter->sampleScatterIsotropic(cache,rng,ekin_i);
  }
  return 0;
}

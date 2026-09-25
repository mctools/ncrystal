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

//Draws a fresh (log-uniform) random energy per call rather than reusing one
//fixed ekin throughout: a single repeated ekin always lands in the same
//cell/grid region, which is unrepresentative of an actual simulation (and
//can make some costs -- e.g. anything that only runs once per distinct
//region and is then effectively free on every subsequent identical call --
//look artificially cheap in a profile).

namespace NC = NCrystal;

int main( int argc, char** argv )
{
  //Parse args:
  if ( argc!=4 ) {
    std::cout << "Please provide args: [cfgstr] [ekin_eV|lo_eV:hi_eV] [nsample]"
              << std::endl;
    return 1;
  }
  const std::string cfgstr(argv[1]);

  const NC::StrView ekinarg(argv[2]);
  double ekin_lo, ekin_hi;
  auto isep = ekinarg.find(':');
  if ( isep == NC::StrView::npos ) {
    auto opt_ekin = ekinarg.toDbl();
    nc_assert_always(opt_ekin.has_value());
    ekin_lo = ekin_hi = opt_ekin.value();
  } else {
    auto opt_lo = ekinarg.substr(0,isep).toDbl();
    auto opt_hi = ekinarg.substr(isep+1).toDbl();
    nc_assert_always(opt_lo.has_value()&&opt_hi.has_value());
    ekin_lo = opt_lo.value();
    ekin_hi = opt_hi.value();
    nc_assert_always(ekin_lo>0.0&&ekin_hi>=ekin_lo);
  }

  auto opt_nsample = NC::StrView(argv[3]).toUInt64();
  nc_assert_always(opt_nsample.has_value());
  const std::uint64_t nsample = opt_nsample.value();

  auto scatter = NC::FactImpl::createScatter(cfgstr);
  auto rng = NC::getRNG();

  std::cout<<"Sampling "<<nsample<<" times at ekin = "
           <<ekin_lo<<" .. "<<ekin_hi<<" eV"<<std::endl;

  const bool varyekin = ( ekin_hi > ekin_lo );
  const double logratio = ( varyekin ? std::log(ekin_hi/ekin_lo) : 0.0 );

  NC::CachePtr cache;
  for ( std::uint64_t i = 0; i < nsample; ++i ) {
    const NC::NeutronEnergy ekin{ NC::DoValidate,
      varyekin ? ekin_lo*std::exp(logratio*rng->generate()) : ekin_lo };
    (void)scatter->sampleScatterIsotropic(cache,rng,ekin);
  }
  return 0;
}

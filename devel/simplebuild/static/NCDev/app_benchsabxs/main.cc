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

//Sibling of app_benchsabsample, benchmarking crossSectionIsotropic instead of
//sampleScatterIsotropic (i.e. cross-section evaluation rather than scattering
//sampling). If ekin_eV is a single value, every call hits the same energy; if
//given as "lo:hi", each call uses a different energy, log-uniformly cycled
//across [lo,hi] over nsample calls, to also exercise the energy-grid lookup
//(as opposed to always hitting the same cached bin/index).

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include <iostream>

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

  std::cout<<"Evaluating cross-sections "<<nsample<<" times at ekin = "
           <<ekin_lo<<" .. "<<ekin_hi<<" eV"<<std::endl;

  NC::CachePtr cache;
  double xssum = 0.0;//accumulate to prevent the call from being optimised away
  const bool varyekin = ( ekin_hi > ekin_lo );
  const double logratio = ( varyekin
                            ? std::log(ekin_hi/ekin_lo) / nsample
                            : 0.0 );
  for ( std::uint64_t i = 0; i < nsample; ++i ) {
    const NC::NeutronEnergy ekin{ NC::DoValidate,
      varyekin ? ekin_lo*std::exp(logratio*i) : ekin_lo };
    xssum += scatter->crossSectionIsotropic(cache,ekin).dbl();
  }
  std::cout<<"(sum of evaluated cross-sections: "<<xssum<<")"<<std::endl;
  return 0;
}

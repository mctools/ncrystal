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

// Tests related to the generation of lists of HKL planes.

#include "NCrystal/NCrystal.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>

namespace NC = NCrystal;

////////////////////////////////////////////////////////////////////////////////
// Bragg threshold and partial HKL lists
////////////////////////////////////////////////////////////////////////////////

// Tests that Info::getBraggThreshold() and Info::hklListPartialCalc(..) work
// for all dcutoffup values, in particular those equal to the lower cutoffs
// used for partial HKL calculations internally (5, 1.5 and 0.75 Aa, which once
// lead to invalid empty-range calculations and an assertion failure), and that
// the Bragg threshold obtained without full initialisation of the HKL lists is
// the same as after a full initialisation.

namespace test_partial {

namespace {
  void checkInvalidRequests()
  {
    auto info = NC::FactImpl::createInfo( "stdlib::Al_sg225.ncmat;"
                                          "dcutoff=0.5;dcutoffup=3" );
    const double nan = std::numeric_limits<double>::quiet_NaN();
    using OD = NC::Optional<double>;
    auto isBad = [&info]( OD dl, OD du )
    {
      try {
        info->hklListPartialCalc( dl, du );
      } catch ( NC::Error::BadInput& ) {
        return true;
      }
      return false;
    };
    nc_assert_always( isBad( nan, NC::NullOpt ) );
    nc_assert_always( isBad( NC::NullOpt, nan ) );
    nc_assert_always( isBad( nan, nan ) );
    nc_assert_always( isBad( 2.0, 1.0 ) );
    //Valid requests outside [0.5,3] give empty lists:
    for ( auto dd : { NC::PairDD(4.0,5.0), NC::PairDD(0.1,0.2),
                      NC::PairDD(3.0,5.0), NC::PairDD(0.1,0.5) } ) {
      auto l = info->hklListPartialCalc( dd.first, dd.second );
      nc_assert_always( l.has_value() && l.value().empty() );
    }
    //Partial lists must be consistent with the full list:
    auto lp = info->hklListPartialCalc( 1.0 );
    std::size_t nexpect = 0;
    for ( auto& e : info->hklList() )
      if ( e.dspacing >= 1.0 )
        ++nexpect;
    nc_assert_always( lp.has_value() && lp.value().size() == nexpect );
    std::cout << "Invalid/out-of-range requests OK ("
              << nexpect << " planes with d>=1Aa)" << std::endl;
  }
}

void run()
{
  for ( auto mat : { "stdlib::Al_sg225.ncmat",
                     "stdlib::AlN_sg186_AluminumNitride.ncmat",
                     "stdlib::PbS_sg225_LeadSulfide.ncmat" } ) {
    for ( auto dcutoffup : { "5", "1.5", "0.75", "1.4",
                             "5.0000001", "4.9999999" } ) {
      std::string cfgstr = ( std::string(mat) + ";dcutoff=0.3;dcutoffup="
                             + dcutoffup );
      //Bragg threshold on fresh object (without full hkl initialisation):
      NC::clearCaches();
      auto bt_partial = NC::FactImpl::createInfo( cfgstr )->getBraggThreshold();
      //Bragg threshold after full initialisation:
      NC::clearCaches();
      auto info = NC::FactImpl::createInfo( cfgstr );
      const auto nhkl = info->hklList().size();
      auto bt_full = info->getBraggThreshold();
      nc_assert_always( bt_partial.has_value() == bt_full.has_value() );
      nc_assert_always( !bt_full.has_value()
                        || bt_partial.value().dbl() == bt_full.value().dbl() );
      std::cout << cfgstr << " : " << nhkl << " HKL groups, Bragg threshold ";
      if ( bt_full.has_value() )
        std::cout << NC::fmt( bt_full.value().dbl(), "%.10g" ) << " Aa"
                  << std::endl;
      else
        std::cout << "none" << std::endl;
      //Degenerate [d,d] ranges must give exactly the planes with that d:
      for ( double d : { 0.3, 0.75, 1.0, 1.5, 10.0 } ) {
        auto l = info->hklListPartialCalc( d, d );
        nc_assert_always( l.has_value() && l.value().empty() );
      }
      for ( auto& e : info->hklList() ) {
        std::size_t nexpect = 0;
        for ( auto& e2 : info->hklList() )
          if ( e2.dspacing == e.dspacing )
            ++nexpect;
        auto l = info->hklListPartialCalc( e.dspacing, e.dspacing );
        nc_assert_always( l.has_value() && l.value().size() == nexpect );
        for ( auto& e2 : l.value() )
          nc_assert_always( e2.dspacing == e.dspacing );
      }
    }
  }
  checkInvalidRequests();
  std::cout << "All OK" << std::endl;
}

}

int main()
{
  std::cout << "==== Bragg threshold and partial HKL lists ====" << std::endl;
  test_partial::run();
  return 0;
}

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
#include "NCrystal/internal/extd_utils/NCFillHKL.hh"
#include "NCrystal/internal/phys_utils/NCEqRefl.hh"
#include <iostream>
#include <map>
#include <set>

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

////////////////////////////////////////////////////////////////////////////////
// HKL families without space group
////////////////////////////////////////////////////////////////////////////////

//Checks that grouping of HKL planes into families for crystals without space
//group does not depend on how many hkl points are buffered before being merged
//into families (FillHKLCfg::max_buffered_points), that it never splits the
//true symmetry families, and that families get the average values of their
//members.

namespace test_nosymbuffer {

namespace {
  struct Fam { unsigned mult; double d; double fsq; };
  using FamMap = std::map<std::vector<NC::HKL>,Fam>;//key: sorted hkl list

  FamMap calcNoSym( const NC::Info& info, double dcutoff, std::size_t nbuf,
                     double merge_tolerance = 1e-6 )
  {
    NC::StructureInfo si = info.getStructureInfo();
    si.spacegroup = 0;
    NC::FillHKLCfg cfg;
    cfg.dcutoff = dcutoff;
    cfg.max_buffered_points = nbuf;
    cfg.merge_tolerance = merge_tolerance;
    FamMap res;
    for ( auto& e : NC::calculateHKLPlanes( si, info.getAtomInfos(), cfg ) ) {
      nc_assert_always( e.explicitValues != nullptr );
      auto v = e.explicitValues->list.get<std::vector<NC::HKL>>();
      nc_assert_always( 2*v.size() == e.multiplicity );
      std::sort( v.begin(), v.end() );
      nc_assert_always( e.hkl == v.front() );
      Fam f{ e.multiplicity, e.dspacing, e.fsquared };
      nc_assert_always( res.emplace( std::move(v), f ).second );
    }
    return res;
  }

  void test( const char * mat, double dcutoff )
  {
    auto info = NC::FactImpl::createInfo( std::string("stdlib::") + mat );
    const FamMap ref = calcNoSym( *info, dcutoff, 0 );
    std::size_t npts = 0;
    for ( auto& e : ref )
      npts += e.first.size();

    //Same families and values for any buffer size:
    for ( std::size_t nbuf : { std::size_t(1000), std::size_t(101),
                               std::size_t(npts), std::size_t(npts-1) } ) {
      const FamMap fm = calcNoSym( *info, dcutoff, nbuf );
      nc_assert_always( fm.size() == ref.size() );
      for ( auto& e : fm ) {
        auto it = ref.find( e.first );
        nc_assert_always( it != ref.end() );
        nc_assert_always( it->second.mult == e.second.mult );
        nc_assert_always( it->second.d == e.second.d );
        nc_assert_always( it->second.fsq == e.second.fsq );
      }
    }

    //Values must be the averages of the members (from a calculation where
    //nothing is merged):
    std::map<NC::HKL,Fam> pt2vals;
    for ( auto& e : calcNoSym( *info, dcutoff, 0, 0.0 ) ) {
      nc_assert_always( e.first.size() == 1 );
      pt2vals[e.first.front()] = e.second;
    }
    nc_assert_always( pt2vals.size() == npts );
    std::size_t nfam_varying = 0;
    for ( auto& e : ref ) {
      NC::StableSum sd, sf;
      double fmin(NC::kInfinity), fmax(-NC::kInfinity);
      for ( auto& hkl : e.first ) {
        const Fam& f = pt2vals.at( hkl );
        sd.add( f.d );
        sf.add( f.fsq );
        fmin = NC::ncmin( fmin, f.fsq );
        fmax = NC::ncmax( fmax, f.fsq );
      }
      const double n = e.first.size();
      nc_assert_always( NC::floateq( e.second.d, sd.sum() / n, 1e-15, 0.0 ) );
      nc_assert_always( NC::floateq( e.second.fsq, sf.sum() / n, 1e-15, 0.0 ) );
      nc_assert_always( e.second.fsq >= fmin && e.second.fsq <= fmax );
      if ( fmin == fmax )
        nc_assert_always( e.second.fsq == fmin );//exact if all agree
      else
        ++nfam_varying;
    }

    //Symmetry families must never be split:
    std::map<NC::HKL,const std::vector<NC::HKL>*> hkl2fam;
    for ( auto& e : ref )
      for ( auto& hkl : e.first )
        hkl2fam[hkl] = &e.first;
    const auto& si = info->getStructureInfo();
    NC::EqRefl sym( static_cast<int>(si.spacegroup) );
    std::set<NC::HKL> symfams;//representatives
    for ( auto& e : ref ) {
      for ( auto& hkl : e.first ) {
        const NC::HKL mhkl( -hkl.h, -hkl.k, -hkl.l );
        auto eqv = sym.getEquivalentReflections( hkl );
        symfams.insert( std::min( eqv.front(),
                                  sym.getEquivalentReflections( mhkl ).front() ) );
        for ( auto& h2 : eqv ) {
          //Found either directly or via the Friedel partner (-h,-k,-l):
          auto it = hkl2fam.find( h2 );
          if ( it == hkl2fam.end() )
            it = hkl2fam.find( NC::HKL( -h2.h, -h2.k, -h2.l ) );
          nc_assert_always( it != hkl2fam.end() && it->second == &e.first );
        }
      }
    }
    const std::size_t nsymfam = symfams.size();
    std::cout << mat << " (dcutoff=" << dcutoff << "): " << npts
              << " hkl points in " << ref.size() << " families ("
              << nsymfam << " symmetry families, " << nfam_varying
              << " with varying F2), OK" << std::endl;
  }
}

void run()
{
  test( "Al_sg225.ncmat", 0.2 );
  test( "GaSe_sg194_GalliumSelenide.ncmat", 0.4 );
  test( "BaO_sg225_BariumOxide.ncmat", 0.35 );
  test( "Al2O3_sg167_Corundum.ncmat", 0.5 );
  test( "UF6_sg62_UraniumHexaflouride.ncmat", 0.8 );
  std::cout << "All OK" << std::endl;
}

}

int main()
{
  std::cout << "==== Bragg threshold and partial HKL lists ====" << std::endl;
  test_partial::run();
  std::cout << "==== HKL families without space group ====" << std::endl;
  test_nosymbuffer::run();
  return 0;
}

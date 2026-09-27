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
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"
#include "NCrystal/internal/utils/NCVector.hh"
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

////////////////////////////////////////////////////////////////////////////////
// Range of HKL indices
////////////////////////////////////////////////////////////////////////////////

//Tests estimateHKLRange: All planes with d>=dcutoff must have indices within
//the returned range (checked by brute force for many lattices, including
//cases where d is exactly dcutoff), and the range must be tight (it is based
//on |h|<=a/d, with equality when G is parallel to the lattice vector a).

namespace test_range {

namespace {

  //Simple deterministic generator of values in [0,1) (splitmix64):
  class SimpleRand {
  public:
    double generate()
    {
      std::uint64_t z = ( m_state += 0x9E3779B97F4A7C15ull );
      z = ( z ^ ( z >> 30 ) ) * 0xBF58476D1CE4E5B9ull;
      z = ( z ^ ( z >> 27 ) ) * 0x94D049BB133111EBull;
      z ^= ( z >> 31 );
      return static_cast<double>( z >> 11 ) * ( 1.0 / 9007199254740992.0 );
    }
  private:
    std::uint64_t m_state = 12345;
  };

  struct Lattice { double a, b, c, alpha, beta, gamma; };

  bool isValid( const Lattice& lt )
  {
    const double ca = std::cos( lt.alpha * NC::kDeg );
    const double cb = std::cos( lt.beta * NC::kDeg );
    const double cg = std::cos( lt.gamma * NC::kDeg );
    return 1.0 - ca*ca - cb*cb - cg*cg + 2.0*ca*cb*cg > 1e-3;
  }

  //Returns max(|h|,|k|,|l|) found by brute force:
  NC::MaxHKL check( const Lattice& lt, double dcutoff )
  {
    const auto mx = NC::estimateHKLRange( dcutoff, lt.a, lt.b, lt.c );
    //Tight (just the 0.1% safety margin on top of the exact bound):
    nc_assert_always( mx.h <= std::max( 1.0, 1.001 * lt.a / dcutoff ) + 1e-9 );
    nc_assert_always( mx.k <= std::max( 1.0, 1.001 * lt.b / dcutoff ) + 1e-9 );
    nc_assert_always( mx.l <= std::max( 1.0, 1.001 * lt.c / dcutoff ) + 1e-9 );
    nc_assert_always( mx.h + 1 > lt.a / dcutoff );
    nc_assert_always( mx.k + 1 > lt.b / dcutoff );
    nc_assert_always( mx.l + 1 > lt.c / dcutoff );
    const auto rec = NC::getReciprocalLatticeRot( lt.a, lt.b, lt.c,
                                                  lt.alpha * NC::kDeg,
                                                  lt.beta * NC::kDeg,
                                                  lt.gamma * NC::kDeg );
    const int bh = 2 * mx.h + 2;
    const int bk = 2 * mx.k + 2;
    const int bl = 2 * mx.l + 2;
    NC::MaxHKL found{ 0, 0, 0 };
    for ( int h = -bh; h <= bh; ++h ) {
      for ( int k = -bk; k <= bk; ++k ) {
        for ( int l = -bl; l <= bl; ++l ) {
          if ( !h && !k && !l )
            continue;
          //Allow for rounding errors in the d-spacing calculation:
          if ( NC::dspacingFromHKL( h, k, l, rec ) < dcutoff * ( 1.0 - 1e-12 ) )
            continue;
          nc_assert_always( std::abs(h) <= mx.h );
          nc_assert_always( std::abs(k) <= mx.k );
          nc_assert_always( std::abs(l) <= mx.l );
          found.h = std::max( found.h, std::abs(h) );
          found.k = std::max( found.k, std::abs(k) );
          found.l = std::max( found.l, std::abs(l) );
        }
      }
    }
    return found;
  }

  //The analytic basis: rows of the lattice matrix (mapping G/2pi to hkl) are
  //the lattice vectors, with lengths a, b and c:
  void checkAnalyticBasis( const Lattice& lt )
  {
    const auto lat = NC::getLatticeRot( lt.a, lt.b, lt.c, lt.alpha * NC::kDeg,
                                        lt.beta * NC::kDeg, lt.gamma * NC::kDeg );
    const auto rec = NC::getReciprocalLatticeRot( lt.a, lt.b, lt.c,
                                                  lt.alpha * NC::kDeg,
                                                  lt.beta * NC::kDeg,
                                                  lt.gamma * NC::kDeg );
    const double len[3] = { lt.a, lt.b, lt.c };
    for ( auto i : NC::ncrange( 3 ) ) {
      //Row i of lat is its product with the unit vector along axis i:
      NC::Vector e( i==0 ? 1.0 : 0.0, i==1 ? 1.0 : 0.0, i==2 ? 1.0 : 0.0 );
      NC::Vector row( 0.0, 0.0, 0.0 );
      for ( auto j : NC::ncrange( 3 ) ) {
        NC::Vector ej( j==0 ? 1.0 : 0.0, j==1 ? 1.0 : 0.0, j==2 ? 1.0 : 0.0 );
        row[j] = e.dot( lat * ej );
      }
      nc_assert_always( NC::floateq( row.mag(), len[i], 1e-12, 0.0 ) );
    }
    //And lat maps G/2pi back to hkl:
    for ( auto& hkl : { NC::Vector(1,0,0), NC::Vector(2,-3,5),
                        NC::Vector(-7,1,4) } ) {
      NC::Vector back = lat * ( rec * hkl );
      back *= ( 1.0 / NC::k2Pi );
      for ( auto j : NC::ncrange( 3 ) )
        nc_assert_always( NC::floateq( back[j], hkl[j], 1e-12, 1e-12 ) );
    }
  }
}

void run()
{
  //Special lattices, with dcutoff giving exact integer ratios (so planes with
  //d exactly equal to dcutoff exist):
  const Lattice special[] = {
    { 4.0, 4.0, 4.0, 90.0, 90.0, 90.0 },//cubic
    { 3.0, 3.0, 5.0, 90.0, 90.0, 120.0 },//hexagonal
    { 5.0, 5.0, 5.0, 70.0, 70.0, 70.0 },//rhombohedral axes
    { 4.0, 6.0, 5.0, 90.0, 105.0, 90.0 },//monoclinic
    { 4.0, 5.0, 6.0, 80.0, 95.0, 110.0 },//triclinic
    { 3.0, 3.0, 40.0, 90.0, 90.0, 90.0 },//long c-axis
  };
  for ( auto& lt : special ) {
    checkAnalyticBasis( lt );
    for ( double dcut : { 1.0, 0.5, 0.25, 0.8, 0.123 } ) {
      if ( lt.c / dcut > 60 )
        continue;//keep brute force fast
      auto mx = NC::estimateHKLRange( dcut, lt.a, lt.b, lt.c );
      auto found = check( lt, dcut );
      std::cout << "a,b,c=" << lt.a << "," << lt.b << "," << lt.c
                << " angles=" << lt.alpha << "," << lt.beta << ","
                << lt.gamma << " dcutoff=" << dcut << " : range=("
                << mx.h << "," << mx.k << "," << mx.l << "), max found=("
                << found.h << "," << found.k << "," << found.l << ")"
                << std::endl;
    }
  }

  //Many random lattices:
  SimpleRand rng;
  unsigned ntested = 0;
  while ( ntested < 200 ) {
    Lattice lt{ 2.0 + 8.0 * rng.generate(), 2.0 + 8.0 * rng.generate(),
                2.0 + 8.0 * rng.generate(), 50.0 + 80.0 * rng.generate(),
                50.0 + 80.0 * rng.generate(), 50.0 + 80.0 * rng.generate() };
    if ( !isValid( lt ) )
      continue;
    ++ntested;
    checkAnalyticBasis( lt );
    check( lt, 0.4 + 0.8 * rng.generate() );
  }
  std::cout << "Checked " << ntested << " random lattices" << std::endl;

  //Minimum value is 1, and huge values are capped:
  auto m1 = NC::estimateHKLRange( 100.0, 3.0, 4.0, 5.0 );
  nc_assert_always( m1.h == 1 && m1.k == 1 && m1.l == 1 );
  auto m2 = NC::estimateHKLRange( 1e-300, 3.0, 4.0, 5.0 );
  nc_assert_always( m2.h == std::numeric_limits<int>::max() );
  std::cout << "All OK" << std::endl;
}

}

int main()
{
  std::cout << "==== Bragg threshold and partial HKL lists ====" << std::endl;
  test_partial::run();
  std::cout << "==== HKL families without space group ====" << std::endl;
  test_nosymbuffer::run();
  std::cout << "==== Range of HKL indices ====" << std::endl;
  test_range::run();
  return 0;
}

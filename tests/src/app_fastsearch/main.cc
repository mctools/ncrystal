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

////////////////////////////////////////////////////////////////////////////////
// Correctness stress-test for NC::fastLowerBoundIdx/fastUpperBoundIdx
// (NCFastSearch.hh): checks, across a wide variety of grid shapes and query
// points, that they return *exactly* the same index as
// (std::lower_bound/upper_bound(first,first+n,val)-first). This is a pure
// correctness check, not a performance benchmark -- see
// nocheck/fastsearch/bench_search.cc (scratch, not part of the repo) for the
// disassembly-verified performance comparison that motivated this. Grid
// shapes deliberately include cases designed to stress the loop's boundary
// arithmetic: n=1, n=2, grids starting exactly at 0.0 or going negative,
// grids with values within a few ULP of each other, and deliberately
// irregular (non-uniformly spaced) grids.
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/utils/NCFastSearch.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/interfaces/NCRNG.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>
#include <algorithm>
#include <sstream>

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  using NC::VectD;
  using NC::ncrange;

  std::uint64_t g_ntrials = 0;

  void checkOne( const VectD& x, double t )
  {
    const auto fullLo = static_cast<std::size_t>(
      std::lower_bound( x.begin(), x.end(), t ) - x.begin() );
    const auto fullUp = static_cast<std::size_t>(
      std::upper_bound( x.begin(), x.end(), t ) - x.begin() );

    const auto fastLo = NC::fastLowerBoundIdx( x.data(), x.size(), t );
    const auto fastUp = NC::fastUpperBoundIdx( x.data(), x.size(), t );

    REQUIRE( fastLo == fullLo );
    REQUIRE( fastUp == fullUp );
    ++g_ntrials;
  }

  void stressGrid( const char* label, const VectD& x )
  {
    REQUIRE( !x.empty() );
    for ( auto v : x )
      checkOne( x, v );
    for ( auto i : ncrange(std::size_t(1),x.size()) )
      checkOne( x, 0.5*(x[i-1]+x[i]) );
    //Out-of-range queries (below/above the grid) are well-defined for
    //std::lower_bound/upper_bound too (0 or n respectively), so check them:
    checkOne( x, x.front() - 1.0 );
    checkOne( x, x.back() + 1.0 );

    auto rng = NC::getRNG();
    const double span = x.back() - x.front();
    for ( unsigned i = 0; i < 3000; ++i )
      checkOne( x, x.front() + rng->generate()*span );
    std::cout << "OK: " << label << " (n=" << x.size() << ")" << std::endl;
  }

  VectD logspaceGrid( double lo, double hi, std::size_t n )
  {
    REQUIRE( lo > 0.0 && hi > lo && n >= 2 );
    VectD x(n);
    const double logratio = std::log(hi/lo) / (n-1);
    for ( auto i : ncrange(n) )
      x[i] = ( i+1==n ? hi : lo*std::exp(i*logratio) );
    return x;
  }

  VectD linspaceGrid( double lo, double hi, std::size_t n )
  {
    REQUIRE( hi > lo && n >= 2 );
    VectD x(n);
    const double step = (hi-lo) / (n-1);
    for ( auto i : ncrange(n) )
      x[i] = ( i+1==n ? hi : lo+i*step );
    return x;
  }

  //A cumulative-sum array with some zero-weight entries, exactly the shape
  //of e.g. PowderBragg::m_fdm_commul (symmetry-forbidden reflections
  //contribute zero, so consecutive cumulative values tie) -- unlike a
  //grid, this is only non-decreasing, not strictly increasing:
  VectD cumulWithTiesGrid( NC::RNG& rng, std::size_t n, double zeroFrac )
  {
    VectD x(n);
    double s = 0.0;
    for ( auto i : ncrange(n) ) {
      if ( rng.generate() >= zeroFrac )
        s += rng.generate();
      x[i] = s;
    }
    return x;
  }

  //Deliberately irregular grid, bounded growth to avoid overflow (see
  //app_gridindex's identical helper -- caught a real inf-overflow bug there
  //with an earlier, unbounded version of this):
  VectD irregularGrid( NC::RNG& rng, std::size_t n )
  {
    VectD x(n);
    double v = 1e-3;
    for ( auto i : ncrange(n) ) {
      x[i] = v;
      v *= ( 1.0 + 0.15*rng.generate() );
    }
    return x;
  }

}

int main()
{
  stressGrid( "single point", VectD{ 1.0 } );
  stressGrid( "two points", VectD{ 1.0, 2.0 } );
  stressGrid( "two points, negative", VectD{ -5.0, -1.0 } );
  stressGrid( "log-uniform, small", logspaceGrid(1e-5,1e3,3) );
  stressGrid( "log-uniform, medium", logspaceGrid(1e-5,1e3,300) );
  stressGrid( "log-uniform, large", logspaceGrid(1e-8,1e8,2400) );
  stressGrid( "linear-uniform, medium", linspaceGrid(0.0,1.0,300) );
  stressGrid( "linear-uniform, from 0", linspaceGrid(0.0,50.0,140) );
  stressGrid( "linear-uniform, negative range", linspaceGrid(-10.0,10.0,140) );

  //Explicit small cases with duplicate/tied values (fastLowerBoundIdx/
  //fastUpperBoundIdx only require non-decreasing, unlike GridIndex's
  //stricter grid requirement -- verified here, not just asserted):
  stressGrid( "all-equal", VectD(50,3.0) );
  stressGrid( "leading duplicates", VectD{1,1,1,1,2,3,4,5} );
  stressGrid( "trailing duplicates", VectD{1,2,3,4,5,5,5,5} );
  stressGrid( "many interior ties", VectD{0,0,0,1,2,3,3,3,3,4,5} );

  auto rng = NC::getRNG();
  //Cumulative-sum-with-ties arrays (PowderBragg::m_fdm_commul's shape),
  //including at the ~1e5 scale PowderBragg tables can reach:
  for ( double zeroFrac : { 0.0, 0.3, 0.7, 0.95 } ) {
    std::ostringstream oss;
    oss << "cumul-with-ties, zeroFrac=" << zeroFrac << ", n=2000";
    stressGrid( oss.str().c_str(), cumulWithTiesGrid(*rng,2000,zeroFrac) );
  }
  stressGrid( "cumul-with-ties, n=150000", cumulWithTiesGrid(*rng,150000,0.3) );

  for ( unsigned trial = 0; trial < 8; ++trial ) {
    std::ostringstream oss;
    oss << "irregular, trial " << trial;
    stressGrid( oss.str().c_str(), irregularGrid(*rng,500) );
  }

  //Grid with values within a few ULP of each other:
  {
    VectD x;
    double v = 1.0;
    for ( unsigned i = 0; i < 50; ++i ) {
      x.push_back(v);
      v = std::nextafter(v,2.0);
    }
    stressGrid( "near-ULP-spaced", x );
  }

  std::cout << "All " << g_ntrials << " fastLowerBoundIdx/fastUpperBoundIdx"
               " checks passed." << std::endl;
  return 0;
}

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
// Tests of reducePtsByEquidistribution, the alternative to the greedy        //
// reducePtsInDistribution (which is tested in app_ptreduce), which selects   //
// points at equal steps of the cumulative integral of an importance density. //
// Apart from the basic contract (number of points, ordering, being a subset  //
// of the input), the tests concern the properties which motivate its         //
// existence: points are concentrated where needed, gaps are bounded, the     //
// result is exactly invariant under scaling of y, and above all is robust    //
// against tiny changes of the input (noise, or the removal of single input   //
// points) which the greedy algorithm in the original reducePtsInDistribution //
// amplifies to macroscopic changes of the result.                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCTestUtils/NCTestFindData.hh"
#include <iostream>
#include <fstream>
#include <functional>

namespace NC=NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  using NC::VectD;
  using Res = std::pair<VectD,VectD>;

  //The whole test battery below is run once per algorithm (see main()), by
  //rebinding this to either reducePtsByEquidistribution or
  //reducePtsByEquidistributionRobust before each run. This means every
  //contract/quality test that reducePtsByEquidistribution has always been
  //run under also exercises reducePtsByEquidistributionRobust:
  using ReduceFct = std::function
    <Res(const VectD&,const VectD&,std::size_t,const NC::PtReduceCfg&)>;
  ReduceFct g_reduce;

  //Deterministic pseudo-random numbers in [-1,1]:
  //Fixme: Do not repeat here, put in NCTestUtils or use NCrystal's own.
  struct Rng {
    std::uint64_t s;
    explicit Rng( std::uint64_t seed = 88172645463325252ULL ) : s(seed) {}
    double u()
    {
      s ^= s << 13;
      s ^= s >> 7;
      s ^= s << 17;
      return static_cast<double>( s % 2000001 ) / 1000000.0 - 1.0;
    }
    std::size_t idx( std::size_t n )
    {
      s ^= s << 13;
      s ^= s >> 7;
      s ^= s << 17;
      return static_cast<std::size_t>( s % n );
    }
  };

  VectD linspaceN( double a, double b, std::size_t n )
  {
    return NC::linspace( a, b, static_cast<unsigned>(n) );
  }

  template<class F>
  VectD applyF( const VectD& x, F f )
  {
    VectD y;
    y.reserve( x.size() );
    for ( double v : x )
      y.push_back( f(v) );
    return y;
  }

  double gauss( double x, double mu, double sigma )
  {
    const double t = ( x - mu ) / sigma;
    return std::exp( -0.5 * t * t );
  }

  //Check the basic contract of any point reduction:
  Res checkedReduce( const VectD& x, const VectD& y, std::size_t k,
                     const NC::PtReduceCfg& cfg = {} )
  {
    auto r = g_reduce( x, y, k, cfg );
    REQUIRE( r.first.size() == std::min( k, x.size() ) );
    REQUIRE( r.second.size() == r.first.size() );
    REQUIRE( r.first.front() == x.front() );
    REQUIRE( r.first.back() == x.back() );
    REQUIRE( r.second.front() == y.front() );
    REQUIRE( r.second.back() == y.back() );
    //increasing, and each point (with its y) exactly one of the input points,
    //in the same order:
    std::size_t j = 0;
    for ( std::size_t i = 0; i < r.first.size(); ++i ) {
      if ( i )
        REQUIRE( r.first[i] > r.first[i-1] );
      while ( j < x.size() && x[j] != r.first[i] )
        ++j;
      REQUIRE( j < x.size() );
      REQUIRE( y[j] == r.second[i] );
      ++j;
    }
    return r;
  }

  //Indices in the input of the selected points:
  std::vector<std::size_t> indices( const VectD& x, const Res& r )
  {
    std::vector<std::size_t> v;
    for ( double e : r.first )
      v.push_back( static_cast<std::size_t>
                   ( std::lower_bound( x.begin(), x.end(), e ) - x.begin() ) );
    return v;
  }

  //Number of points in a which do not exist in b (both sorted):
  std::size_t nUnmatched( const VectD& a, const VectD& b )
  {
    std::size_t n = 0;
    for ( double v : a )
      n += ( std::binary_search( b.begin(), b.end(), v ) ? 0 : 1 );
    return n;
  }

  //Max error of the linear interpolant through the reduced points, relative to
  //max y (lin), or as max of |ln(f/y)| over points with y above the threshold
  //relative to the max (log):
  struct Errs { double lin, log; };
  Errs interpolationErrors( const VectD& x, const VectD& y, const Res& r,
                            double logthreshold = 1e-8 )
  {
    const double ymax = *std::max_element( y.begin(), y.end() );
    Errs e{ 0.0, 0.0 };
    std::size_t j = 0;
    for ( std::size_t i = 0; i < x.size(); ++i ) {
      while ( j + 2 < r.first.size() && r.first[j+1] < x[i] )
        ++j;
      const double f = ( x[i] - r.first[j] ) / ( r.first[j+1] - r.first[j] );
      const double v = r.second[j] + f * ( r.second[j+1] - r.second[j] );
      e.lin = std::max( e.lin, std::abs( v - y[i] ) / ymax );
      if ( y[i] > logthreshold * ymax && v > 0.0 )
        e.log = std::max( e.log, std::abs( std::log( v / y[i] ) ) );
      else if ( y[i] > logthreshold * ymax )
        e.log = std::max( e.log, 1e10 );
    }
    return e;
  }

  double maxGap( const Res& r )
  {
    double g = 0.0;
    for ( std::size_t i = 1; i < r.first.size(); ++i )
      g = std::max( g, r.first[i] - r.first[i-1] );
    return g;
  }

  std::size_t countInRange( const Res& r, double lo, double hi )
  {
    std::size_t n = 0;
    for ( double v : r.first )
      n += ( v >= lo && v <= hi ) ? 1 : 0;
    return n;
  }

  //A test function with peaks of different widths on a small background and
  //an exponential tail:
  double testFct( double x )
  {
    return ( 0.02
             + 1.0 * gauss( x, 300.0, 12.0 )
             + 0.5 * gauss( x, 520.0, 30.0 )
             + 0.3 * gauss( x, 700.0, 5.0 )
             + 0.6 * std::exp( -x / 150.0 ) );
  }

  void testContractManySizes()
  {
    //All combinations of small input sizes and requested point counts:
    Rng rng( 12345 );
    for ( std::size_t n = 2; n <= 40; ++n ) {
      VectD x = linspaceN( -3.0, 7.5, n );
      VectD y = applyF( x, [](double v){
        return 1.0 + std::sin(2.0*v)*std::sin(2.0*v);
      } );
      for ( std::size_t k = 2; k <= n + 3; ++k )
        checkedReduce( x, y, k );
      //non-uniform x, random non-negative y with exact zeros:
      VectD x2, y2;
      double xv = -1.0;
      for ( std::size_t i = 0; i < n; ++i ) {
        xv += 0.01 + 0.5 * ( 1.0 + rng.u() );
        x2.push_back( xv );
        const double r = rng.u();
        y2.push_back( r < -0.3 ? 0.0 : 10.0 * std::abs( r ) );
      }
      y2[ n / 2 ] += 1.0;//ensure not all zero
      for ( std::size_t k = 2; k <= n + 1; ++k )
        checkedReduce( x2, y2, k );
    }
    std::cout << "contract for many sizes ok" << std::endl;
  }

  void testTargetAtLeastInputSizeIsNoop()
  {
    VectD x = linspaceN( 0.0, 1.0, 50 );
    VectD y = applyF( x, [](double v){ return v*v + 0.1; } );
    for ( std::size_t k : { 50u, 51u, 100u, 100000u } ) {
      auto r = g_reduce( x, y, k, {} );
      REQUIRE( r.first == x );
      REQUIRE( r.second == y );
    }
    std::cout << "no-op when target >= input size ok" << std::endl;
  }

  void testTwoAndThreePoints()
  {
    VectD x = linspaceN( 0.0, 10.0, 101 );
    VectD y = applyF( x, [](double v){ return gauss( v, 3.0, 0.5 ) + 0.01; } );
    auto r2 = checkedReduce( x, y, 2 );
    const VectD expected{ 0.0, 10.0 };
    REQUIRE( r2.first == expected );
    auto r3 = checkedReduce( x, y, 3 );
    REQUIRE( r3.first.size() == 3 );
    REQUIRE( r3.first[1] > 0.0 && r3.first[1] < 10.0 );
    std::cout << "two and three points ok" << std::endl;
  }

  void testConstantFunctionGivesEquidistantPoints()
  {
    //With no curvature at all (in either linear or logarithmic scale), only the
    //constant part of the density remains, and the result must be as
    //equidistant as the input allows. Note that only a constant has this
    //property: for instance a linear ramp has non-zero curvature in ln(y).
    VectD x = linspaceN( -5.0, 5.0, 1001 );
    for ( double c : { 3.0, 1e-30, 1e30 } ) {
      VectD y( x.size(), c );
      auto r = checkedReduce( x, y, 11 );
      auto idx = indices( x, r );
      for ( std::size_t i = 0; i < idx.size(); ++i )
        REQUIRE( idx[i] == 100*i );
    }
    //Same on an input grid with non-equidistant points. Now the requirement
    //is that the points are equidistant in x, as well as possible:
    VectD xg = NC::geomspace( 1.0, 1000.0, 2001 );
    VectD yg( xg.size(), 2.0 );
    auto r = checkedReduce( xg, yg, 11 );
    for ( std::size_t i = 1; i + 1 < r.first.size(); ++i ) {
      const double expected = 1.0 + 999.0 * static_cast<double>(i) / 10.0;
      //closest input point to the expected position:
      auto it = std::lower_bound( xg.begin(), xg.end(), expected );
      double closest = *it;
      if ( it != xg.begin() && std::abs( *std::prev(it) - expected )
           < std::abs( closest - expected ) )
        closest = *std::prev(it);
      REQUIRE( r.first[i] == closest );
    }
    std::cout << "constant function gives equidistant points ok" << std::endl;
  }

  void testExponentialHasNoLogCurvature()
  {
    //For y=exp(-a*x), ln(y) is a straight line and only the linear scale
    //contributes to the density, so points concentrate at small x (where y
    //has the largest curvature) but the result is still sensible:
    VectD x = linspaceN( 0.0, 100.0, 5001 );
    VectD y = applyF( x, [](double v){ return std::exp( -0.2 * v ); } );
    auto r = checkedReduce( x, y, 30 );
    REQUIRE( countInRange( r, 0.0, 20.0 ) >= 12 );
    REQUIRE( maxGap( r ) < 100.0 / ( 0.15 * 29.0 ) + 1.0 );
    std::cout << "exponential function ok" << std::endl;
  }

  void testExactScaleInvariance()
  {
    //Multiplying y by a power of 2 is exact in floating point, and so must give
    //an identical selection. Other factors must give the same selection up to
    //rounding errors:
    VectD x = linspaceN( 0.0, 1000.0, 5000 );
    VectD y = applyF( x, testFct );
    auto ref = checkedReduce( x, y, 60 );
    for ( double scale : { 0.25, 4.0, 1024.0, 1.0/1048576.0,
                           std::ldexp(1.0,200), std::ldexp(1.0,-200) } ) {
      VectD ys( y );
      for ( auto& e : ys )
        e *= scale;
      auto r = checkedReduce( x, ys, 60 );
      REQUIRE( r.first == ref.first );
    }
    for ( double scale : { 3.0, 1e-20, 1e20, 7.77e100 } ) {
      VectD ys( y );
      for ( auto& e : ys )
        e *= scale;
      auto r = checkedReduce( x, ys, 60 );
      REQUIRE( nUnmatched( r.first, ref.first ) <= 1 );
    }
    std::cout << "scale invariance ok" << std::endl;
  }

  void testSymmetry()
  {
    VectD x = linspaceN( -50.0, 50.0, 4001 );
    VectD y = applyF( x, [](double v)
    {
      return ( 0.01 + gauss( v, 0.0, 4.0 )
               + 0.3*gauss( v, 25.0, 2.0 )
               + 0.3*gauss( v, -25.0, 2.0 ) );
    } );
    auto r = checkedReduce( x, y, 41 );
    auto idx = indices( x, r );
    const std::size_t n = x.size();
    for ( std::size_t i = 0; i < idx.size(); ++i ) {
      const std::size_t mirror = n - 1 - idx[ idx.size() - 1 - i ];
      const std::size_t d = mirror > idx[i] ? mirror - idx[i] : idx[i] - mirror;
      REQUIRE( d <= 1 );
    }
    std::cout
      << "symmetric input gives (almost) symmetric selection ok" << std::endl;
  }

  void testPointsConcentrateWhereFunctionVaries()
  {
    //A single Gaussian peak on a flat background. A uniform selection of 41
    //points would place about 2.5 points within +-3 sigma of the peak
    VectD x = linspaceN( 0.0, 1000.0, 5001 );
    VectD y = applyF( x, [](double v)
    {
      return 0.001 + gauss( v, 400.0, 10.0 );
    } );
    auto r = checkedReduce( x, y, 41 );
    REQUIRE( countInRange( r, 370.0, 430.0 ) >= 15 );
    //but the rest should be spread out over the flat parts, not ignored:
    REQUIRE( countInRange( r, 0.0, 300.0 ) >= 1 );
    REQUIRE( countInRange( r, 500.0, 1000.0 ) >= 1 );
    std::cout << "points concentrate where the function varies ok" << std::endl;
  }

  void testMaxGapGuarantee()
  {
    //Regardless of how much the function varies, the equidistant part of the
    //density guarantees no gap much larger than L/(fraction*(N-1)):
    VectD x = linspaceN( 0.0, 1000.0, 20001 );
    VectD y = applyF( x, testFct );
    VectD y2 = applyF( x,
                       [](double v){ return 1e-6 + gauss( v, 500.0, 0.5 ); } );
    for ( double frac : { 0.05, 0.15, 0.3, 0.6, 1.0 } ) {
      for ( std::size_t k : { 20u, 50u, 200u } ) {
        NC::PtReduceCfg cfg;
        cfg.equidistant_fraction = frac;
        for ( auto* yy : { &y, &y2 } ) {
          auto r = checkedReduce( x, *yy, k, cfg );
          const double bound = 1000.0 / ( frac * static_cast<double>(k-1) );
          REQUIRE( maxGap( r ) <= bound + 1.0 );//+1 for input point spacing
        }
      }
    }
    //Without the equidistant part, there is no such guarantee, but it
    //should still work:
    NC::PtReduceCfg cfg0;
    cfg0.equidistant_fraction = 0.0;
    checkedReduce( x, y, 50, cfg0 );
    checkedReduce( x, y2, 50, cfg0 );
    std::cout << "max gap guarantee ok" << std::endl;
  }

  void testEquidistantFractionOneIgnoresFunction()
  {
    VectD x = linspaceN( 0.0, 1000.0, 2001 );
    VectD y1 = applyF( x, testFct );
    VectD y2 = applyF( x, [](double v){ return 1.0 + 0.1*std::sin(v); } );
    NC::PtReduceCfg cfg;
    cfg.equidistant_fraction = 1.0;
    auto r1 = checkedReduce( x, y1, 26, cfg );
    auto r2 = checkedReduce( x, y2, 26, cfg );
    REQUIRE( r1.first == r2.first );
    auto idx = indices( x, r1 );
    for ( std::size_t i = 0; i < idx.size(); ++i )
      REQUIRE( idx[i] == 80*i );
    std::cout << "equidistant_fraction=1 gives equidistant points ok"
              << std::endl;
  }

  void testStepFunctionPutsPointsAtTheStep()
  {
    VectD x = linspaceN( 0.0, 1000.0, 5001 );
    VectD y = applyF( x, [](double v){ return v < 500.0 ? 0.1 : 1.0; } );
    auto r = checkedReduce( x, y, 30 );
    //The step is a single input interval (of width 0.2), and the points in the
    //vicinity of it should resolve it:
    REQUIRE( countInRange( r, 499.0, 501.2 ) >= 2 );
    const double e = interpolationErrors( x, y, r ).lin;
    //Some point must sit on each side of the step, close enough:
    REQUIRE( e < 0.6 );
    std::cout << "step function ok" << std::endl;
  }

  void testZeroRegionsAndAllZero()
  {
    VectD x = linspaceN( 0.0, 100.0, 1001 );
    //Zero except for a bump:
    VectD y = applyF( x, [](double v)
    {
      return ( ( v > 40.0 && v < 50.0 )
               ? std::sin( (v-40.0) * 0.314159265 ) : 0.0 );
    } );
    auto r = checkedReduce( x, y, 25 );
    for ( double v : r.second )
      REQUIRE( std::isfinite( v ) );
    REQUIRE( countInRange( r, 38.0, 52.0 ) >= 6 );
    //All zero (no function at all): equidistant points
    VectD z( x.size(), 0.0 );
    auto rz = checkedReduce( x, z, 11 );
    auto idx = indices( x, rz );
    for ( std::size_t i = 0; i < idx.size(); ++i )
      REQUIRE( idx[i] == 100*i );
    //A single non-zero value:
    VectD one( x.size(), 0.0 );
    one[500] = 1.0;
    checkedReduce( x, one, 15 );
    //Non-zero only at the ends:
    VectD ends( x.size(), 0.0 );
    ends.front() = 1.0;
    ends.back() = 2.0;
    checkedReduce( x, ends, 15 );
    std::cout << "zeros and all-zero input ok" << std::endl;
  }

  void testNonUniformInputGrid()
  {
    //Same smooth function sampled on very different grids should give
    //similar results (concentrating where the function varies, and with
    //similar error), since only the function matters and not the sampling:
    VectD xa = linspaceN( 0.0, 1000.0, 6001 );
    VectD xb = NC::geomspace( 0.5, 1000.0, 6001 );
    xb.front() = 0.0;
    VectD xc;//dense in some places, sparse elsewhere
    for ( double v = 0.0; v < 1000.0;
          v += ( v > 280.0 && v < 320.0 ? 0.05 : 0.5 ) )
      xc.push_back( v );
    xc.push_back( 1000.0 );
    std::vector<double> nearpeak;
    for ( auto* xx : { &xa, &xb, &xc } ) {
      VectD y = applyF( *xx, testFct );
      auto r = checkedReduce( *xx, y, 60 );
      nearpeak.push_back( static_cast<double>( countInRange( r,
                                                             270.0, 330.0 ) ) );
      const auto e = interpolationErrors( *xx, y, r );
      REQUIRE( e.lin < 0.05 );
    }
    for ( double v : nearpeak )
      REQUIRE( std::abs( v - nearpeak.front() ) <= 3.0 );
    std::cout << "non-uniform input grids ok" << std::endl;
  }

  void testTailsAreResolvedThroughLogTerm()
  {
    //A function with a strong peak and a long tail spanning many orders of
    //magnitude. Points in the tail are needed to describe ln(y):
    VectD x = linspaceN( 0.0, 100.0, 10001 );
    VectD y = applyF( x, [](double v)
    {
      return std::exp( -0.02*v*v ) + 1e-12*std::exp( -0.005*v );
    } );
    auto r = checkedReduce( x, y, 60 );
    const auto e = interpolationErrors( x, y, r, 1e-15 );
    REQUIRE( r.first.size() == 60 );
    REQUIRE( countInRange( r, 20.0, 100.0 ) >= 6 );
    REQUIRE( e.log < 10.0 );//sanity: tail values are not wildly off
    //With the log term ignored for values under a high tail_floor, the same
    //data must still work (the tail is considered flat):
    NC::PtReduceCfg cfg;
    cfg.tail_floor = 1e-6;
    auto r2 = checkedReduce( x, y, 60, cfg );
    REQUIRE( countInRange( r2, 0.0, 10.0 ) >= 20 );
    std::cout << "tails and tail_floor ok" << std::endl;
  }

  void testTailFloorMakesSmallValuesIrrelevant()
  {
    //Values far below the tail_floor (relative to the maximum) have no
    //influence on the density. tail_floor now acts as a smooth floor
    //((f^4+tail_floor^4)^(1/4), continuous derivative) rather than a hard
    //clamp (max(tail_floor,f)) -- see docs/claude_session_vdos_fma_reprod.md
    //for why a hard clamp was replaced (it gave sub-floor noise an
    //artificial single-point curvature spike right at the clamp boundary).
    //The quartic (rather than quadratic, i.e. sqrt(f^2+tail_floor^2))
    //combination was needed to keep this test's exact-irrelevance guarantee:
    //a quadratic smooth floor let oscillating sub-floor patterns leak
    //through with just enough amplitude to shift several selected points.
    //The two inputs here are identical where the function is above 1e-6,
    //and differ only below 1e-11 with different patterns, so the result
    //must be the same.
    VectD x = linspaceN( 0.0, 100.0, 4001 );
    NC::PtReduceCfg cfg;
    cfg.tail_floor = 1e-8;
    auto make = [&x]( double amp, double freq ) {
      return applyF( x, [amp,freq](double v){
        const double g = gauss( v, 20.0, 5.0 );
        return g > 1e-6 ? g : amp * ( 1.0 + std::sin( freq * v ) );
      } );
    };
    auto r1 = checkedReduce( x, make( 1e-11, 1.0 ), 40, cfg );
    auto r2 = checkedReduce( x, make( 1e-14, 3.7 ), 40, cfg );
    auto r3 = checkedReduce( x, make( 0.0, 1.0 ), 40, cfg );
    REQUIRE( r1.first == r2.first );
    REQUIRE( r1.first == r3.first );
    std::cout << "tail_floor ok" << std::endl;
  }

  void testRobustToNoiseInY()
  {
    //Relative noise at the level of rounding errors on all y-values must not
    //change the selection (the greedy algorithm can change several points).
    VectD x = linspaceN( 0.0, 1000.0, 8000 );
    VectD y = applyF( x, testFct );
    auto ref = checkedReduce( x, y, 60 );
    Rng rng;
    std::size_t nchanged = 0, maxdiff = 0;
    for ( double sigma : { 1e-16, 1e-15, 1e-14, 1e-13 } ) {
      for ( int trial = 0; trial < 25; ++trial ) {
        VectD yn( y );
        for ( auto& e : yn )
          e *= 1.0 + sigma * rng.u();
        auto r = checkedReduce( x, yn, 60 );
        const auto d = nUnmatched( r.first, ref.first );
        nchanged += ( d > 0 ? 1 : 0 );
        maxdiff = std::max( maxdiff, d );
      }
    }
    //Rare 1-point changes (where a quantile is exactly between two points) are
    //fine, but no more than that:
    REQUIRE( maxdiff <= 1 );
    REQUIRE( nchanged <= 10 );
    std::cout << "robust to noise in y ok" << std::endl;
  }

  void testRobustToRemovalOfSinglePoints()
  {
    //Removing a single input point is what happens when the edge of a tail is
    //moved by a bin due to noise. It must at most move a single output point,
    //and not change the selection of points far away from it at all.
    VectD x = linspaceN( 0.0, 1000.0, 20000 );
    VectD y = applyF( x, testFct );
    auto ref = checkedReduce( x, y, 60 );
    Rng rng( 424242 );
    std::size_t maxdiff = 0, ntested = 0, nchanged = 0;
    for ( int trial = 0; trial < 60; ++trial ) {
      const std::size_t idx = 1 + rng.idx( x.size() - 2 );
      if ( std::binary_search( ref.first.begin(), ref.first.end(), x[idx] ) )
        continue;//removing a selected point legitimately changes the result
      VectD x2( x ), y2( y );
      x2.erase( x2.begin() + static_cast<std::ptrdiff_t>(idx) );
      y2.erase( y2.begin() + static_cast<std::ptrdiff_t>(idx) );
      auto r = checkedReduce( x2, y2, 60 );
      const auto d = nUnmatched( r.first, ref.first );
      maxdiff = std::max( maxdiff, d );
      nchanged += ( d > 0 ? 1 : 0 );
      ++ntested;
    }
    REQUIRE( ntested >= 40 );
    REQUIRE( maxdiff <= 1 );
    REQUIRE( nchanged <= ntested / 3 );
    std::cout << "robust to removal of single input points ok" << std::endl;
  }

  void testMoreRobustThanGreedyAlgorithm()
  {
    //Same experiments as above, comparing directly to the greedy algorithm,
    //on data of the kind for which it is fragile (long fine-grained smooth
    //input, and reduction to few points).
    VectD x = linspaceN( 0.0, 20.0, 20000 );
    VectD y = applyF( x, [](double v){ return 1.0/(1.0+v*v)
          + 0.2/(1.0+4.0*(v-9.0)*(v-9.0)); } );
    auto refEq = checkedReduce( x, y, 30 );
    auto refGr = NC::reducePtsInDistribution( x, y, 30 );
    Rng rng( 999 );
    double sumEq = 0.0, sumGr = 0.0;
    std::size_t n = 0;
    for ( int trial = 0; trial < 40; ++trial ) {
      const std::size_t idx = 1 + rng.idx( x.size() - 2 );
      if ( std::binary_search( refEq.first.begin(),
                               refEq.first.end(), x[idx] )
           || std::binary_search( refGr.first.begin(),
                                  refGr.first.end(), x[idx] ) )
        continue;
      VectD x2( x ), y2( y );
      x2.erase( x2.begin() + static_cast<std::ptrdiff_t>(idx) );
      y2.erase( y2.begin() + static_cast<std::ptrdiff_t>(idx) );
      sumEq += static_cast<double>
        ( nUnmatched( checkedReduce( x2, y2, 30 ).first,
                      refEq.first ) );
      sumGr += static_cast<double>
        ( nUnmatched( NC::reducePtsInDistribution( x2, y2, 30 ).first,
                      refGr.first ) );
      ++n;
    }
    REQUIRE( n >= 25 );
    //On average, the number of changed points is much lower:
    REQUIRE( sumEq < 0.2 * sumGr + 1.0 );
    std::cout << "much more robust than the greedy algorithm ok" << std::endl;
  }

  void testContinuousDependenceOnInput()
  {
    //A smoothly moving peak must move the selected points smoothly: no
    //selected point jumps by more than a few input spacings for a tiny
    //(1e-3 sigma) shift of the peak (much smaller than the input spacing, so
    //some points move to the neighbouring input point, but no more than that)
    VectD x = linspaceN( 0.0, 100.0, 10001 );
    auto sel = [&x]( double mu ) {
      VectD y = applyF( x,
                        [mu](double v){ return 0.001 + gauss( v, mu, 2.0 ); } );
      return checkedReduce( x, y, 40 );
    };
    for ( double mu0 : { 30.0, 47.3, 61.11 } ) {
      auto a = sel( mu0 );
      auto b = sel( mu0 + 0.002 );
      for ( std::size_t i = 0; i < a.first.size(); ++i )
        REQUIRE( std::abs( a.first[i] - b.first[i] ) < 0.03 );
    }
    std::cout << "continuous dependence on the input ok" << std::endl;
  }

  void testQualityComparedToOtherApproaches()
  {
    //The selected points must describe the function well, in comparison to the
    //other algorithms, in both linear and logarithmic scale
    struct Case { const char* name; VectD x; VectD y; std::size_t k; };
    std::vector<Case> cases;
    {
      VectD x = linspaceN( 0.0, 1000.0, 10001 );
      cases.push_back( { "peaks", x, applyF( x, testFct ), 60 } );
      cases.push_back( { "peaks-few", x, applyF( x, testFct ), 25 } );
    }
    {
      VectD x = linspaceN( -20.0, 20.0, 8001 );
      cases.push_back( { "gauss-tail", x,
            applyF( x, [](double v)
            {
              return std::exp(-0.5*v*v) + 1e-9*std::exp(-0.1*std::abs(v));
            } ), 50 } );
    }
    {
      VectD x = linspaceN( 0.0, 30.0, 6001 );
      cases.push_back( { "poisson-like", x,
            applyF( x, [](double v)
            {
              return std::exp( 12.0*std::log( v + 1e-300 ) - v ) / 479001600.0;
            } ), 40 } );
    }
    for ( auto& c : cases ) {
      auto eq = checkedReduce( c.x, c.y, c.k );
      auto gr = NC::reducePtsInDistribution( c.x, c.y, c.k );
      //Uniform selection as a baseline:
      Res uni;
      for ( std::size_t i = 0; i < c.k; ++i ) {
        const std::size_t idx = ( i * ( c.x.size() - 1 ) ) / ( c.k - 1 );
        uni.first.push_back( c.x[idx] );
        uni.second.push_back( c.y[idx] );
      }
      const auto ee = interpolationErrors( c.x, c.y, eq );
      const auto eg = interpolationErrors( c.x, c.y, gr );
      const auto eu = interpolationErrors( c.x, c.y, uni );
      //Not worse than a few times the greedy algorithm (which optimises
      //directly), and clearly better than uniform selection in linear scale:
      REQUIRE( ee.lin < 4.0 * eg.lin + 1e-3 );
      REQUIRE( ee.lin < eu.lin );
      REQUIRE( ee.log < 4.0 * eg.log + 1.0 );
    }
    std::cout << "quality ok" << std::endl;
  }

  void testResultCanBeReducedAgain()
  {
    //Applying the reduction to its own output with the same target size does
    //nothing, and reducing in two steps still gives a good result
    VectD x = linspaceN( 0.0, 1000.0, 10000 );
    VectD y = applyF( x, testFct );
    auto r1 = checkedReduce( x, y, 200 );
    auto r1b = g_reduce( r1.first, r1.second, 200, {} );
    REQUIRE( r1b.first == r1.first );
    auto r2 = checkedReduce( r1.first, r1.second, 50 );
    const auto e = interpolationErrors( x, y, r2 );
    REQUIRE( e.lin < 0.1 );
    std::cout << "repeated reduction ok" << std::endl;
  }

  void testLargeInput()
  {
    VectD x = linspaceN( 0.0, 1000.0, 300000 );
    VectD y = applyF( x, testFct );
    auto r = checkedReduce( x, y, 500 );
    REQUIRE( interpolationErrors( x, y, r ).lin < 0.01 );
    std::cout << "large input ok" << std::endl;
  }

  void testRandomizedContract()
  {
    //Random inputs of random size, with random spikes, zeros and dynamic range:
    Rng rng( 2024 );
    for ( int trial = 0; trial < 300; ++trial ) {
      const std::size_t n = 2 + rng.idx( 400 );
      VectD x;
      double xv = rng.u() * 100.0;
      for ( std::size_t i = 0; i < n; ++i ) {
        xv += 1e-3 + std::abs( rng.u() ) * ( trial % 3 == 0 ? 1e3 : 1.0 );
        x.push_back( xv );
      }
      VectD y;
      for ( std::size_t i = 0; i < n; ++i ) {
        const double r = rng.u();
        double v = std::abs( r ) * std::pow( 10.0, 8.0 * rng.u() );
        if ( r < -0.6 )
          v = 0.0;
        y.push_back( v );
      }
      NC::PtReduceCfg cfg;
      cfg.equidistant_fraction = 0.5 * ( 1.0 + rng.u() );
      cfg.tail_floor = std::pow( 10.0, -30.0 * std::abs( rng.u() ) - 1.0 );
      const std::size_t k = 2 + rng.idx( n + 3 );
      checkedReduce( x, y, k, cfg );
    }
    std::cout << "randomised inputs ok" << std::endl;
  }

  //Full battery of contract/quality tests above, run against whichever
  //algorithm g_reduce is currently bound to (see main()):
  void runSharedTestBattery()
  {
    testContractManySizes();
    testTargetAtLeastInputSizeIsNoop();
    testTwoAndThreePoints();
    testConstantFunctionGivesEquidistantPoints();
    testExponentialHasNoLogCurvature();
    testExactScaleInvariance();
    testSymmetry();
    testPointsConcentrateWhereFunctionVaries();
    testMaxGapGuarantee();
    testEquidistantFractionOneIgnoresFunction();
    testStepFunctionPutsPointsAtTheStep();
    testZeroRegionsAndAllZero();
    testNonUniformInputGrid();
    testTailsAreResolvedThroughLogTerm();
    testTailFloorMakesSmallValuesIrrelevant();
    testRobustToNoiseInY();
    testRobustToRemovalOfSinglePoints();
    testMoreRobustThanGreedyAlgorithm();
    testContinuousDependenceOnInput();
    testQualityComparedToOtherApproaches();
    testResultCanBeReducedAgain();
    testLargeInput();
    testRandomizedContract();
  }

  ////////////////////////////////////////////////////////////////////////////
  // Real (VDOS/SAB-derived) data tests. These load the actual (x,y) arrays  //
  // captured from NCVDOSKnlGrid.cc's setupE0ABGrid call to                  //
  // reducePtsByEquidistribution while investigating a real cross-platform   //
  // sabxs instability (see docs/claude_session_vdos_fma_reprod.md): the     //
  // known-bad case (li_from_li2o_e0grid, where a real CI run showed a       //
  // last-few-ULP difference in the input flipping which discrete point got  //
  // selected), plus six more from other materials/elements so a fix         //
  // targeted at the first cannot silently trade its robustness for theirs.  //
  ////////////////////////////////////////////////////////////////////////////

  struct RealDataset { const char* name; VectD x, y; std::size_t targetN; };

  RealDataset loadRealDataset( const char* name )
  {
    RealDataset ds;
    ds.name = name;
    std::string fn = std::string(name) + ".txt";
    std::ifstream ifs( nctest::find_test_data( "ptreduce", fn.c_str() ) );
    REQUIRE( bool(ifs) );
    std::size_t n;
    ifs >> n >> ds.targetN;
    ds.x.resize(n);
    ds.y.resize(n);
    for ( auto& v : ds.x ) ifs >> v;
    for ( auto& v : ds.y ) ifs >> v;
    REQUIRE( NC::nc_is_grid( ds.x ) );
    return ds;
  }

  std::vector<const char*> realDatasetNames()
  {
    return { "li_from_li2o_e0grid", "al_e0grid", "cu_e0grid", "pb_e0grid",
             "lih_e0grid", "polyethylene_e0grid_c", "polyethylene_e0grid_h" };
  }

  //Indices in x selected by a result (by value, since the result is a subset
  //of x with unique, ordered elements):
  std::vector<std::size_t> selectedIndices( const VectD& x, const VectD& rx )
  {
    std::vector<std::size_t> v;
    std::size_t j = 0;
    for ( double e : rx ) {
      while ( j < x.size() && x[j] != e )
        ++j;
      REQUIRE( j < x.size() );
      v.push_back( j++ );
    }
    return v;
  }

  std::size_t nDiffIndices( const std::vector<std::size_t>& a,
                            const std::vector<std::size_t>& b )
  {
    REQUIRE( a.size() == b.size() );
    std::size_t d = 0;
    for ( std::size_t i = 0; i < a.size(); ++i )
      d += ( a[i] != b[i] ) ? 1 : 0;
    return d;
  }

  double linInterpErrRelToMax( const VectD& x, const VectD& y,
                               const Res& r )
  {
    const double ymax = *std::max_element( y.begin(), y.end() );
    double e = 0.0;
    std::size_t j = 0;
    for ( std::size_t i = 0; i < x.size(); ++i ) {
      while ( j + 2 < r.first.size() && r.first[j+1] < x[i] )
        ++j;
      const double f = ( x[i] - r.first[j] ) / ( r.first[j+1] - r.first[j] );
      const double v = r.second[j] + f * ( r.second[j+1] - r.second[j] );
      e = std::max( e, std::abs( v - y[i] ) / ymax );
    }
    return e;
  }

  //For a given algorithm and dataset, count how many of ntrials independent
  //relative perturbations of y at level sigma change the selection at all
  //(nchanged), and the largest number of selected points that differ from
  //the unperturbed reference in any single trial (maxdiff):
  struct NoiseResult { std::size_t nchanged, maxdiff; };
  NoiseResult noiseRobustness( const ReduceFct& fct, const RealDataset& ds,
                              double sigma, int ntrials = 50,
                              std::uint64_t seed = 12345 )
  {
    auto ref = fct( ds.x, ds.y, ds.targetN, {} );
    auto refIdx = selectedIndices( ds.x, ref.first );
    Rng rng( seed );
    NoiseResult res{ 0, 0 };
    for ( int t = 0; t < ntrials; ++t ) {
      VectD yn( ds.y );
      for ( auto& v : yn )
        v *= ( 1.0 + sigma * rng.u() );
      auto r = fct( ds.x, yn, ds.targetN, {} );
      auto idx = selectedIndices( ds.x, r.first );
      const auto d = nDiffIndices( refIdx, idx );
      res.nchanged += ( d > 0 ? 1 : 0 );
      res.maxdiff = std::max( res.maxdiff, d );
    }
    return res;
  }

  ReduceFct oldFct()
  {
    return []( const VectD& x, const VectD& y, std::size_t k,
              const NC::PtReduceCfg& cfg )
    { return NC::reducePtsByEquidistribution( x, y, k, cfg ); };
  }

  ReduceFct robustFct()
  {
    return []( const VectD& x, const VectD& y, std::size_t k,
              const NC::PtReduceCfg& cfg )
    { return NC::reducePtsByEquidistributionRobust( x, y, k, 3, cfg ); };
  }

  void testRealDataKnownBugCaseIsFixed()
  {
    //This is the exact scenario traced from a real CI run: at a noise level
    //matching what was actually observed cross-platform (~1e-8 relative,
    //see docs/claude_session_vdos_fma_reprod.md), the standard algorithm
    //occasionally (not always -- this is inherently a near-tie) selects a
    //different point than with clean input, while the robust one does not,
    //here and at one order of magnitude higher (1e-7) noise as well:
    auto ds = loadRealDataset( "li_from_li2o_e0grid" );
    auto old_ = oldFct();
    auto rob = robustFct();
    const auto oldAt8 = noiseRobustness( old_, ds, 1e-8 );
    const auto robAt8 = noiseRobustness( rob, ds, 1e-8 );
    REQUIRE( oldAt8.nchanged >= 1 );//the known bug: reproduced here
    REQUIRE( robAt8.nchanged == 0 );//fixed
    const auto oldAt7 = noiseRobustness( old_, ds, 1e-7 );
    const auto robAt7 = noiseRobustness( rob, ds, 1e-7 );
    REQUIRE( oldAt7.nchanged >= 1 );
    REQUIRE( robAt7.maxdiff <= oldAt7.maxdiff );
    std::cout << "known bug case (li_from_li2o_e0grid) fixed ok" << std::endl;
  }

  void testRealDataAcrossMaterialsNoRegression()
  {
    //Same noise-robustness check as above, across several other materials'
    //real data: the robust algorithm must not be meaningfully *less* robust
    //than the standard one on any of them (allowing equality, since several
    //of these are already stable at this noise level with either algorithm),
    //and must not degrade interpolation quality by more than a generous,
    //explicitly bounded factor -- catching a fix for one material's
    //robustness silently regressing another's, either in stability or in
    //quality:
    double sumOldMaxdiff = 0.0, sumRobMaxdiff = 0.0;
    auto old_ = oldFct();
    auto rob = robustFct();
    for ( auto* name : realDatasetNames() ) {
      auto ds = loadRealDataset( name );
      for ( double sigma : { 1e-8, 1e-7 } ) {
        const auto o = noiseRobustness( old_, ds, sigma );
        const auto r = noiseRobustness( rob, ds, sigma );
        sumOldMaxdiff += static_cast<double>( o.maxdiff );
        sumRobMaxdiff += static_cast<double>( r.maxdiff );
      }
      auto refOld = old_( ds.x, ds.y, ds.targetN, {} );
      auto refRob = rob( ds.x, ds.y, ds.targetN, {} );
      const double eOld = linInterpErrRelToMax( ds.x, ds.y, refOld );
      const double eRob = linInterpErrRelToMax( ds.x, ds.y, refRob );
      //Generous bound: quality on real data may vary either way with a
      //wider curvature stencil (it trades a little sharp-feature resolution
      //for noise robustness -- see NCMath.hh's doc comment), but must not
      //collapse:
      REQUIRE( eRob < 4.0 * eOld + 1e-6 );
    }
    //In aggregate (summed over all seven materials and both noise levels),
    //the robust algorithm must be substantially more stable, not just on
    //the one known-bad case:
    std::cout << "  [real data] sum(maxdiff) old=" << sumOldMaxdiff
              << " robust=" << sumRobMaxdiff << std::endl;
    REQUIRE( sumRobMaxdiff < 0.5 * sumOldMaxdiff );
    std::cout << "real data across materials: no regression ok" << std::endl;
  }

}

int main()
{
  g_reduce = oldFct();
  runSharedTestBattery();
  g_reduce = robustFct();
  runSharedTestBattery();

  testRealDataKnownBugCaseIsFixed();
  testRealDataAcrossMaterialsNoRegression();
  return 0;
}

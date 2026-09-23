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
// Tests of VDOS code which is sensitive to floating point contraction (fused //
// multiply-add, FMA), which can make results differ between platforms. See   //
// also app_fmaprobe and app_fmamath.                                         //
//                                                                            //
// Most tests compare results bit-for-bit against reference implementations  //
// in which the rounding of intermediate products is enforced by a volatile   //
// store (see rnd() below). Thus the tests currently define the semantics as  //
// "unfused". If a site is deliberately changed to use an explicit std::fma,  //
// the corresponding reference here must be changed to match. The last tests  //
// print full grids from a complete VDOS expansion, in a precision which      //
// must be reproducible across all platforms.                                 //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/vdos/NCVDOSUtils.hh"
#include "NCrystal/internal/vdos/NCVDOSGn.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/vdos/NCVDOSExpand.hh"
#include "NCrystal/internal/vdos/NCVDOSKnlGrid.hh"
#include "NCrystal/internal/vdos/NCVDOSToScatKnl.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>

namespace NC=NCrystal;
namespace NCV=NCrystal::VDOS;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  //Forces rounding of an intermediate result to double precision, so the
  //compiler can not fuse the operation producing it with a later one:
  inline double rnd( double x ) { volatile double v = x; return v; }

  void requireIdentical( const NC::VectD& a, const NC::VectD& b,
                         const char* what )
  {
    if ( a.size() != b.size() ) {
      std::cout << "ERROR: size mismatch in " << what << std::endl;
      REQUIRE( false );
    }
    for ( auto i : NC::ncrange(a.size()) ) {
      if ( a[i] != b[i] ) {
        std::cout << "ERROR: mismatch in " << what << " at index " << i
                  << ": " << NC::fmt(a[i],"%.17g") << " vs "
                  << NC::fmt(b[i],"%.17g") << std::endl;
        REQUIRE( false );
      }
    }
  }

  void prtgrid( const char* name, const NC::VectD& g )
  {
    std::cout << name << " (" << g.size() << " points):";
    for ( double e : g )
      std::cout << " " << NC::fmt(e,"%.9g");
    std::cout << std::endl;
  }

  NC::VDOSData makeDebyeVDOS()
  {
    //Al-like
    return NC::createVDOSDebye( NC::DebyeTemperature{410.0},
                                NC::Temperature{293.15},
                                NC::SigmaBound{1.5},
                                NC::AtomMass{26.98} );
  }

  void testGnEval()
  {
    //Gn::eval and Gn::evalMany are two implementations of the same
    //interpolation and must give bit-identical results. Also the upper edge of
    //the energy range is calculated in the constructor.
    auto vd = makeDebyeVDOS();
    NC::VDOSEval ve(vd);
    NCV::VDOSGn gn( ve, NCV::VDOSGn::Cfg::Default );
    gn.growMaxOrder( 6 );
    for ( unsigned n = 1; n <= 6; ++n ) {
      NCV::VDOSGn::Order order(n);
      auto er = gn.eRange( order );
      const double bw = gn.binWidth( order );
      const auto& spec = gn.getRawSpectrum( order );
      REQUIRE( er.second
               == NCV::equidistantGridPoint( er.first, bw, spec.size()-1 ) );
      auto grid = NC::linspace( er.first - 3.3*bw, er.second + 2.1*bw, 3001 );
      NC::VectD out, work;
      gn.evalMany( order, grid, out, work );
      REQUIRE( out.size() == grid.size() );
      for ( auto i : NC::ncrange(grid.size()) ) {
        if ( out[i] != gn.eval( order, grid[i] ) ) {
          std::cout << "ERROR: eval and evalMany differ for order " << n
                    << " at E=" << NC::fmt(grid[i],"%.17g") << ": "
                    << NC::fmt(out[i],"%.17g") << " vs "
                    << NC::fmt(gn.eval( order, grid[i] ),"%.17g")
                    << std::endl;
          REQUIRE( false );
        }
      }
    }
    std::cout << "Gn eval/evalMany consistency ok" << std::endl;
  }

  NCV::PWLFct makePWL( double x0, double bw, NC::VectD vals )
  {
    NCV::PWLFct f;
    f.x0 = x0;
    f.binWidth = bw;
    f.dataHolder = std::move(vals);
    f.f = f.dataHolder;
    return f;
  }

  void testPWLMoveSemantics()
  {
    auto f = makePWL( 0.1, 0.2, {1.0,2.0,3.0,4.0} );
    const double* p = f.dataHolder.data();
    REQUIRE( f.f.data() == p );

    //Move construction keeps data and re-points span:
    NCV::PWLFct g( std::move(f) );
    REQUIRE( g.dataHolder.size() == 4 );
    REQUIRE( g.f.data() == g.dataHolder.data() );
    REQUIRE( g.f.size() == 4 );
    REQUIRE( g.x0 == 0.1 && g.binWidth == 0.2 );
    REQUIRE( f.f.empty() && f.dataHolder.empty() );

    //Move assignment:
    NCV::PWLFct h;
    h = std::move(g);
    REQUIRE( h.f.data() == h.dataHolder.data() );
    REQUIRE( h.f.size() == 4 && h.f[3] == 4.0 );
    REQUIRE( g.f.empty() && g.dataHolder.empty() );

    //Non-owning function, referring to external data, keeps referring to it:
    NC::VectD ext{ 1.0, 2.0, 3.0 };
    NCV::PWLFct e;
    e.x0 = 0.0;
    e.binWidth = 1.0;
    e.f = ext;
    NCV::PWLFct e2( std::move(e) );
    REQUIRE( e2.f.data() == ext.data() && e2.dataHolder.empty() );

    //Growing a vector of them (reallocations) keeps everything valid:
    std::vector<NCV::PWLFct> v;
    for ( int i = 0; i < 50; ++i ) {
      v.push_back( makePWL( i, 0.5, {double(i),i+1.0,i+2.0} ) );
      for ( auto& fi : v ) {
        REQUIRE( fi.f.data() == fi.dataHolder.data() );
        REQUIRE( fi.f.size() == 3 );
      }
    }
    for ( int i = 0; i < 50; ++i )
      REQUIRE( v[i].f[0] == double(i) );
    std::cout << "PWLFct move semantics ok" << std::endl;
  }

  NC::VectD refEvalPWLSum( const std::vector<NCV::PWLFct>& fs,
                           const NC::VectD& grid, const NC::VectD& ws )
  {
    //Mirror of the algorithm in evalPWLSum, but with unfused arithmetic:
    NC::VectD out( grid.size(), 0.0 );
    auto xAt = []( const NCV::PWLFct& p, std::size_t i )
    {
      return std::fma( p.binWidth, static_cast<double>(i), p.x0 );
    };
    for ( std::size_t n = 0; n < fs.size(); ++n ) {
      const auto& p = fs[n];
      const auto& f = p.f;
      const std::size_t lastBin = f.size() - 1;
      const double invbw = 1.0 / p.binWidth;
      const double weight = ws.empty() ? 1.0 : ws[n];
      const double xmax = xAt( p, lastBin );
      const double xtol = 1e-9 * p.binWidth;
      std::size_t g = static_cast<std::size_t>
        ( std::lower_bound( grid.begin(), grid.end(), p.x0 - xtol )
          - grid.begin() );
      if ( g == grid.size() || grid[g] > xmax + xtol )
        continue;
      double xLeft = p.x0;
      for ( std::size_t i = 0; i < lastBin && g < grid.size(); ++i ) {
        const double xRight = xAt( p, i + 1 );
        const double y0 = f[i];
        const double y1 = f[i+1];
        const double slope = ( y1 - y0 ) * invbw;
        const double ylo = std::min( y0, y1 );
        const double yhi = std::max( y0, y1 );
        std::size_t end = g;
        while ( end < grid.size() && grid[end] < xRight )
          ++end;
        for ( std::size_t k = g; k < end; ++k ) {
          double v = y0 + rnd( slope * ( grid[k] - xLeft ) );
          v = std::max( ylo, std::min( yhi, v ) );
          out[k] = out[k] + rnd( weight * v );
        }
        g = end;
        xLeft = xRight;
      }
      if ( g < grid.size() && grid[g] <= xmax + xtol ) {
        const double y = f[lastBin];
        std::size_t end = g;
        while ( end < grid.size() && grid[end] <= xmax + xtol )
          ++end;
        for ( std::size_t k = g; k < end; ++k )
          out[k] = out[k] + rnd( weight * y );
      }
    }
    return out;
  }

  void testEvalPWLSum()
  {
    auto genVals = []( unsigned n, unsigned zeroEvery, double a )
    {
      NC::VectD v;
      for ( unsigned i = 0; i < n; ++i )
        v.push_back( ( i % zeroEvery == 0 || i + 1 == n )
                     ? 0.0 : ( 1.0 + i * a ) / ( 2.0 + i % 7 ) );
      return v;
    };
    const double bw = 0.0173;
    std::vector<NCV::PWLFct> fs;
    fs.push_back( makePWL( -0.37, bw, genVals( 40, 5, 0.137 ) ) );
    fs.push_back( makePWL( 0.0031, 0.0091, genVals( 60, 1000000, 0.071 ) ) );
    fs.push_back( makePWL( -0.2, bw, genVals( 45, 4, 0.0313 ) ) );
    const NC::VectD ws{ 0.7, 1.3, 0.11 };

    //Grid: dense, plus the exact node positions of all the functions:
    NC::VectD grid = NC::linspace( -0.5, 0.6, 1537 );
    for ( auto& f : fs )
      for ( std::size_t i = 0; i < f.f.size(); ++i )
        grid.push_back( f.xAt(i) );
    std::sort( grid.begin(), grid.end() );
    grid.erase( std::unique( grid.begin(), grid.end() ), grid.end() );

    auto res = NCV::evalPWLSum( fs, grid, ws );
    requireIdentical( res, refEvalPWLSum( fs, grid, ws ), "evalPWLSum" );
    auto res_nows = NCV::evalPWLSum( fs, grid );
    requireIdentical( res_nows, refEvalPWLSum( fs, grid, {} ),
                      "evalPWLSum (no weights)" );

    //Each function on its own: At all node positions the value must be exactly
    //the node value (not just approximately, and not tiny non-zero numbers at
    //zero-valued nodes, which is what happens if the bin edges are not
    //calculated consistently with the node positions or if the interpolation
    //suffers from cancellation or fusion). Also, values must be non-negative
    //(the node values are) and never above the maximum node value.
    for ( auto& f : fs ) {
      std::vector<NCV::PWLFct> one;
      one.push_back( makePWL( f.x0, f.binWidth,
                              NC::VectD( f.f.begin(), f.f.end() ) ) );
      auto r1 = NCV::evalPWLSum( one, grid );
      const double ymax = *std::max_element( f.f.begin(), f.f.end() );
      for ( double v : r1 )
        REQUIRE( v >= 0.0 && v <= ymax );
      std::size_t nzero = 0;
      for ( std::size_t i = 0; i < f.f.size(); ++i ) {
        const double xn = f.xAt(i);
        const auto k = static_cast<std::size_t>
          ( std::lower_bound( grid.begin(), grid.end(), xn )
            - grid.begin() );
        REQUIRE( k < grid.size() && grid[k] == xn );
        nzero += ( f.f[i] == 0.0 ? 1 : 0 );
        if ( r1[k] != f.f[i] ) {
          std::cout << "ERROR: value " << NC::fmt(r1[k],"%.17g")
                    << " at node " << i << " of piecewise linear function"
                    << " with node value " << NC::fmt(f.f[i],"%.17g")
                    << std::endl;
          REQUIRE( false );
        }
      }
      REQUIRE( nzero >= 1 );
    }
    std::cout << "evalPWLSum ok" << std::endl;
  }

  void testEvalPWLSumEndpoints()
  {
    //A function with non-zero values at its endpoints, evaluated on grids with
    //points which are exactly at, or a few ulps from, the endpoints. This is
    //what happens when the grid points come from nodes of another function
    //with the same endpoint, calculated in a slightly different way. The value
    //at the endpoint is very different from the (zero) value outside, so it
    //must not depend on rounding errors.
    std::vector<NCV::PWLFct> fs;
    fs.push_back( makePWL( -0.37, 0.0173, {0.5,0.9,1.3,1.1,0.7} ) );
    const double xlo = fs.front().x0;
    const double xhi = fs.front().x1();
    for ( int nulp = -3; nulp <= 3; ++nulp ) {
      double a = xlo, b = xhi;
      for ( int i = 0; i < std::abs(nulp); ++i ) {
        a = std::nextafter( a, nulp < 0 ? -NC::kInfinity : NC::kInfinity );
        b = std::nextafter( b, nulp < 0 ? -NC::kInfinity : NC::kInfinity );
      }
      auto r = NCV::evalPWLSum( fs, NC::VectD{ a, b } );
      //Values a few ulps outside/inside the endpoints must be the endpoint
      //values, up to the tiny interpolation error for points inside:
      REQUIRE( std::abs( r[0] - 0.5 ) < 1e-12 );
      REQUIRE( std::abs( r[1] - 0.7 ) < 1e-12 );
    }
    //But not for points clearly outside:
    auto r = NCV::evalPWLSum( fs, NC::VectD{ xlo - 1e-6, xhi + 1e-6 } );
    REQUIRE( r[0] == 0.0 && r[1] == 0.0 );
    std::cout << "evalPWLSum endpoints ok" << std::endl;
  }

  void testPWLNarrowToPos()
  {
    const double bw = 0.0173;
    struct Case { double x0; unsigned n; };
    const Case cases[] = { {-0.37,80}, {-0.0173,10}, {0.00001,30},
                           {0.5,12}, {-5.0,301} };
    for ( auto& c : cases ) {
      NC::VectD v;
      for ( unsigned i = 0; i < c.n; ++i )
        v.push_back( 1.0 + 0.01*i );
      auto f = makePWL( c.x0, bw, v );
      auto r = NCV::pwlNarrowToPos( f );
      const double t = 1e-3 * bw;
      std::size_t i = ( f.x0 < t
                        ? static_cast<std::size_t>( std::ceil((t-f.x0)/bw) )
                        : 0u );
      REQUIRE( r.f.size() == c.n - i );
      REQUIRE( r.f.data() == r.dataHolder.data() );
      REQUIRE( r.x0 == NCV::equidistantGridPoint( c.x0, bw, i ) );
      REQUIRE( r.x0 >= t );
      REQUIRE( r.f.front() == v[i] && r.f.back() == v.back() );
    }
    std::cout << "pwlNarrowToPos ok" << std::endl;
  }

  void testTrimEquidistantGridUpperEdge()
  {
    //NB: EquidistantGrid::x1() uses an explicit std::fma. The function applies
    //a tolerance of 1e-9*binWidth, to keep nodes which mathematically are at
    //xmax regardless of rounding errors:
    NCV::EquidistantGrid g{ -0.37, 0.0173, 200 };
    const std::size_t targets[] = { 2, 3, 17, 100, 199, 200 };
    for ( auto t : targets ) {
      NCV::EquidistantGrid tg = g;
      tg.npts = t;
      const double xmax = tg.x1();
      auto trimmed = [&g]( double x )
      {
        auto r = g;
        NCV::trimEquidistantGridUpperEdge( r, x );
        return r;
      };
      auto r = trimmed( xmax );
      REQUIRE( r.npts == t );
      REQUIRE( r.x0 == g.x0 && r.binWidth == g.binWidth );

      //A few ulps around the node is still at the node:
      double lo = xmax, hi = xmax;
      for ( int i = 0; i < 4; ++i ) {
        lo = std::nextafter( lo, -NC::kInfinity );
        hi = std::nextafter( hi, NC::kInfinity );
        REQUIRE( trimmed( lo ).npts == t );
        REQUIRE( trimmed( hi ).npts == t );
      }
      //But not when clearly below (or above) the node:
      REQUIRE( trimmed( xmax - 1e-6*g.binWidth ).npts == ( t > 2 ? t - 1 : 2 ) );
      REQUIRE( trimmed( xmax + 1e-6*g.binWidth ).npts == t );
      REQUIRE( trimmed( xmax - 0.5*g.binWidth ).npts == ( t > 2 ? t - 1 : 2 ) );
    }
    auto r = g;
    NCV::trimEquidistantGridUpperEdge( r, -100.0 );
    REQUIRE( r.npts == 2 );
    r = g;
    NCV::trimEquidistantGridUpperEdge( r, 100.0 );
    REQUIRE( r.npts == 200 );
    std::cout << "trimEquidistantGridUpperEdge ok" << std::endl;
  }

  void testMakeCommonGrid()
  {
    //Points in the result must be exactly the points we get from the
    //individual grids, which are given by xAt(i):
    std::vector<NCV::EquidistantGrid> grids = {
      { -0.37, 0.0173, 120 }, { 0.0031, 0.0091, 200 }, { -0.2, 0.0173, 45 },
      { 0.4, 0.0517, 30 } };
    NC::VectD cand;
    for ( auto& g : grids ) {
      for ( std::size_t i = 0; i < g.npts; ++i )
        cand.push_back( g.xAt(i) );
    }
    std::sort( cand.begin(), cand.end() );
    auto res = NCV::makeCommonGrid( grids );
    REQUIRE( NC::nc_is_grid( res ) );
    for ( double x : res ) {
      if ( !std::binary_search( cand.begin(), cand.end(), x ) ) {
        std::cout << "ERROR: makeCommonGrid returned the point "
                  << NC::fmt(x,"%.17g")
                  << " which is not exactly a point on one of the input"
          " grids" << std::endl;
        REQUIRE( false );
      }
    }
    REQUIRE( res.front() == cand.front() );
    REQUIRE( res.back() == cand.back() );
    std::cout << "makeCommonGrid ok (" << res.size() << " points)"
              << std::endl;
  }

  void testTrimTailByIntegral()
  {
    //Mirror of the algorithm with unfused arithmetic:
    auto ref = []( NC::VectD x, NC::VectD y, double frac )
    {
      NC::StableSumKahan sum;
      for ( std::size_t i = 1; i < x.size(); ++i )
        sum.add( rnd( ( x[i] - x[i-1] ) * ( y[i-1] + y[i] ) ) );
      const double limit = frac * ( sum.sum() * 0.5 );
      double removed = 0.0;
      while ( x.size() > 2 ) {
        const auto i = x.size() - 2;
        const double dx = x[i+1] - x[i];
        const double area = 0.5 * rnd( dx * ( y[i] + y[i+1] ) );
        if ( removed + area > limit )
          break;
        removed += area;
        x.pop_back();
        y.pop_back();
      }
      return x.size();
    };
    auto x0 = NC::linspace( 0.0, 10.0, 501 );
    NC::VectD y0;
    for ( double x : x0 )
      y0.push_back( 1.0 / ( 1.0 + 4.0*(x-3.0)*(x-3.0) )
                    + 1e-3 / ( 1.0 + x*x ) );
    for ( double frac : { 1e-12, 1e-9, 1e-6, 1e-3, 0.05, 0.5 } ) {
      auto x = x0;
      auto y = y0;
      NCV::trimTailByIntegral( x, y, frac );
      REQUIRE( x.size() == y.size() );
      REQUIRE( x.size() == ref( x0, y0, frac ) );
      std::cout << "trimTailByIntegral(frac=" << NC::fmt(frac,"%g")
                << ") keeps " << x.size() << " points" << std::endl;
    }
  }

  void testPerturbationRobustness()
  {
    //Changing the input by an amount comparable to rounding errors (1e-14
    //relative) should only change the outputs by a similar amount. This is
    //also what happens when the same calculation is carried out on different
    //platforms (e.g. with or without fused multiply-add instructions). Any
    //discrete decision which is sensitive to rounding errors (a comparison
    //with a threshold which a value happens to be exactly at, a value at the
    //edge of a function, a number of grid points, ...) can amplify such a
    //change into a macroscopic one, in particular on the special grids where
    //many values are mathematically exactly on boundaries.
    struct Cfg { double debyeT, temp, mass; int lux; };
    const Cfg cfgs[] = { {410.0,293.15,26.98,2000},
                         {200.0,10.0,207.2,2000},
                         {900.0,600.0,12.0,2001} };
    struct Res { unsigned order; NC::Rectangle sab; NC::VectD a, b, e0; };
    auto run = []( const Cfg& c, double debyeT )
    {
      auto vd = NC::createVDOSDebye( NC::DebyeTemperature{debyeT},
                                     NC::Temperature{c.temp},
                                     NC::SigmaBound{1.5},
                                     NC::AtomMass{c.mass} );
      const NCV::VDOSLux lux( c.lux );
      auto gnexpn = NCV::expandVDOSToGnFcts( vd, lux );
      auto dim = NCV::gridDimFromLux( lux );
      auto ab = NCV::determineAlphaBetaGridFromGn( gnexpn, dim.first,
                                                   dim.second );
      return Res{ gnexpn.Gn.maxOrder().value(), gnexpn.sabRange,
                  ab.first, ab.second,
                  NCV::setupE0ABGrid( gnexpn, 20 ).first };
    };
    auto requireClose = []( const NC::VectD& a, const NC::VectD& b,
                            const char* what )
    {
      REQUIRE( a.size() == b.size() );
      const double scale = NC::ncmax( NC::ncabs(a.front()), NC::ncabs(a.back()) );
      for ( auto i : NC::ncrange(a.size()) ) {
        if ( !( NC::ncabs( a[i] - b[i] ) < 1e-9 * scale ) ) {
          std::cout << "ERROR: " << what << " grid point " << i
                    << " changed from " << NC::fmt(a[i],"%.15g") << " to "
                    << NC::fmt(b[i],"%.15g")
                    << " after a perturbation of the input by ~1e-14"
                    << std::endl;
          REQUIRE( false );
        }
      }
    };
    //Deterministic pseudo-random relative perturbations in [-1e-14,1e-14]:
    std::uint64_t rng = 88172645463325252ULL;
    auto nextPerturbation = [&rng]()
    {
      rng ^= rng << 13;
      rng ^= rng >> 7;
      rng ^= rng << 17;
      return 1e-14 * ( static_cast<double>( rng % 2001 ) / 1000.0 - 1.0 );
    };
    for ( auto& c : cfgs ) {
      const auto base = run( c, c.debyeT );
      for ( int i = 0; i < 5; ++i ) {
        const auto p = run( c, c.debyeT * ( 1.0 + nextPerturbation() ) );
        REQUIRE( p.order == base.order );
        REQUIRE( NC::ncabs( p.sab.y0() - base.sab.y0() )
                 < 1e-9 * NC::ncabs( base.sab.y0() ) );
        REQUIRE( NC::ncabs( p.sab.y1() - base.sab.y1() )
                 < 1e-9 * NC::ncabs( base.sab.y1() ) );
        REQUIRE( NC::ncabs( p.sab.x1() - base.sab.x1() )
                 < 1e-9 * NC::ncabs( base.sab.x1() ) );
        requireClose( base.a, p.a, "alpha" );
        requireClose( base.b, p.b, "beta" );
        requireClose( base.e0, p.e0, "E0" );
      }
    }
    std::cout << "perturbation robustness ok" << std::endl;
  }

  void testExpansionGrids()
  {
    //Full expansion + grid construction for a Debye VDOS. Output must be
    //reproducible across all platforms (to the precision printed):
    auto vd = makeDebyeVDOS();
    const NCV::VDOSLux lux( 2000 );//lowest next-gen level
    auto gnexpn = NCV::expandVDOSToGnFcts( vd, lux );
    std::cout << "Expansion: maxOrder=" << gnexpn.Gn.maxOrder().value()
              << " suggestedEmax=" << NC::fmt(gnexpn.suggestedEmax.dbl(),"%.9g")
              << " alpha2x=" << NC::fmt(gnexpn.alpha2x,"%.9g")
              << " sabRange=[" << NC::fmt(gnexpn.sabRange.x0(),"%.9g")
              << "," << NC::fmt(gnexpn.sabRange.x1(),"%.9g")
              << "]x[" << NC::fmt(gnexpn.sabRange.y0(),"%.9g")
              << "," << NC::fmt(gnexpn.sabRange.y1(),"%.9g")
              << "]" << std::endl;

    auto comb = NCV::getCombinedGnFct( gnexpn );
    REQUIRE( comb.first.size() == comb.second.size() );
    std::cout << "Combined Gn: " << comb.first.size() << " points" << std::endl;
    prtgrid( "Combined Gn beta grid (first 40)",
             NC::VectD( comb.first.begin(), comb.first.begin() + 40 ) );

    auto e0 = NCV::setupE0ABGrid( gnexpn, 20 );
    prtgrid( "E0 grid", e0.first );

    unsigned nalpha, nbeta;
    std::tie( nalpha, nbeta ) = NCV::gridDimFromLux( lux );
    auto ab = NCV::determineAlphaBetaGridFromGn( gnexpn, nalpha, nbeta );
    prtgrid( "Alpha grid", ab.first );
    prtgrid( "Beta grid", ab.second );

    auto kd = NCV::createScatteringKernel( vd, lux );
    double sabsum = 0.0;
    for ( double s : kd.sab )
      sabsum += s;
    std::cout << "Kernel: nalpha=" << kd.alphaGrid.size()
              << " nbeta=" << kd.betaGrid.size()
              << " suggestedEmax=" << NC::fmt(kd.suggestedEmax,"%.9g")
              << " sum(sab)=" << NC::fmt(sabsum,"%.5g") << std::endl;
  }

}

int main()
{
  testGnEval();
  testPWLMoveSemantics();
  testEvalPWLSum();
  testEvalPWLSumEndpoints();
  testPWLNarrowToPos();
  testTrimEquidistantGridUpperEdge();
  testMakeCommonGrid();
  testTrimTailByIntegral();
  testPerturbationRobustness();
  testExpansionGrids();
  return 0;
}

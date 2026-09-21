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
#include <algorithm>

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
               == er.first + rnd( static_cast<double>(spec.size()-1) * bw ) );
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
    for ( std::size_t n = 0; n < fs.size(); ++n ) {
      const auto& p = fs[n];
      const auto& f = p.f;
      const std::size_t lastBin = f.size() - 1;
      const double invbw = 1.0 / p.binWidth;
      const double weight = ws.empty() ? 1.0 : ws[n];
      const double xmax = p.x0 + rnd( static_cast<double>(lastBin)
                                      * p.binWidth );
      std::size_t g = static_cast<std::size_t>
        ( std::lower_bound( grid.begin(), grid.end(), p.x0 )
          - grid.begin() );
      if ( g == grid.size() || grid[g] > xmax )
        continue;
      for ( std::size_t i = 0; i < lastBin && g < grid.size(); ++i ) {
        const double xLeft = p.x0 + rnd( static_cast<double>(i)
                                         * p.binWidth );
        const double xRight = xLeft + p.binWidth;
        const double y0 = f[i];
        const double y1 = f[i+1];
        const double slope = ( y1 - y0 ) * invbw;
        const double intercept = y0 - rnd( slope * xLeft );
        std::size_t end = g;
        while ( end < grid.size() && grid[end] < xRight )
          ++end;
        for ( std::size_t k = g; k < end; ++k )
          out[k] = out[k] + rnd( weight
                                 * ( intercept + rnd( slope * grid[k] ) ) );
        g = end;
      }
      if ( g < grid.size() && grid[g] <= xmax ) {
        const double y = f[lastBin];
        std::size_t end = g;
        while ( end < grid.size() && grid[end] <= xmax )
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
        grid.push_back( f.x0 + rnd( static_cast<double>(i) * f.binWidth ) );
    std::sort( grid.begin(), grid.end() );
    grid.erase( std::unique( grid.begin(), grid.end() ), grid.end() );

    auto res = NCV::evalPWLSum( fs, grid, ws );
    requireIdentical( res, refEvalPWLSum( fs, grid, ws ), "evalPWLSum" );
    auto res_nows = NCV::evalPWLSum( fs, grid );
    requireIdentical( res_nows, refEvalPWLSum( fs, grid, {} ),
                      "evalPWLSum (no weights)" );

    //Nodes where a function is exactly zero must evaluate to exactly zero, not
    //to tiny (possibly negative) numbers from rounding in the interpolation
    //formula. Test with each function on its own. NB: For this we use dyadic
    //x0 and binWidth values, so all node positions and bin boundaries are
    //represented exactly (with general values the bin boundaries computed by
    //accumulation and the node positions computed by multiplication can differ
    //by an ulp, in which case a node can end up in the neighbouring bin and
    //be evaluated with a rounding error there). The products slope*xLeft still
    //require rounding, so fusion is still detected.
    std::vector<NCV::PWLFct> dy;
    dy.push_back( makePWL( -0.375, 0.015625, genVals( 40, 5, 0.137 ) ) );
    dy.push_back( makePWL( 0.0078125, 0.0078125, genVals( 60, 7, 0.071 ) ) );
    dy.push_back( makePWL( -0.25, 0.015625, genVals( 45, 4, 0.0313 ) ) );
    NC::VectD dgrid = NC::linspace( -0.5, 0.6, 1537 );
    for ( auto& f : dy )
      for ( std::size_t i = 0; i < f.f.size(); ++i )
        dgrid.push_back( f.x0 + static_cast<double>(i) * f.binWidth );
    std::sort( dgrid.begin(), dgrid.end() );
    dgrid.erase( std::unique( dgrid.begin(), dgrid.end() ), dgrid.end() );
    requireIdentical( NCV::evalPWLSum( dy, dgrid, ws ),
                      refEvalPWLSum( dy, dgrid, ws ),
                      "evalPWLSum (dyadic)" );
    for ( auto& f : dy ) {
      std::vector<NCV::PWLFct> one;
      one.push_back( makePWL( f.x0, f.binWidth,
                              NC::VectD( f.f.begin(), f.f.end() ) ) );
      auto r1 = NCV::evalPWLSum( one, dgrid );
      std::size_t nzero = 0;
      for ( std::size_t i = 0; i < f.f.size(); ++i ) {
        const double xn = f.x0 + static_cast<double>(i) * f.binWidth;
        const auto k = static_cast<std::size_t>
          ( std::lower_bound( dgrid.begin(), dgrid.end(), xn )
            - dgrid.begin() );
        REQUIRE( k < dgrid.size() && dgrid[k] == xn );
        if ( f.f[i] == 0.0 ) {
          ++nzero;
          if ( r1[k] != 0.0 ) {
            std::cout << "ERROR: non-zero value "
                      << NC::fmt(r1[k],"%.17g")
                      << " at zero-valued node of piecewise linear function"
                      << std::endl;
            REQUIRE( false );
          }
        }
      }
      REQUIRE( nzero >= 1 );
    }
    std::cout << "evalPWLSum ok" << std::endl;
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
      REQUIRE( r.x0 == c.x0 + rnd( static_cast<double>(i) * bw ) );
      REQUIRE( r.x0 >= t );
      REQUIRE( r.f.front() == v[i] && r.f.back() == v.back() );
    }
    std::cout << "pwlNarrowToPos ok" << std::endl;
  }

  void testTrimEquidistantGridUpperEdge()
  {
    //NB: EquidistantGrid::x1() uses an explicit std::fma:
    NCV::EquidistantGrid g{ -0.37, 0.0173, 200 };
    const std::size_t targets[] = { 2, 3, 17, 100, 199, 200 };
    for ( auto t : targets ) {
      NCV::EquidistantGrid tg = g;
      tg.npts = t;
      const double xmax = tg.x1();
      auto r = g;
      NCV::trimEquidistantGridUpperEdge( r, xmax );
      REQUIRE( r.npts == t );
      REQUIRE( r.x0 == g.x0 && r.binWidth == g.binWidth );

      //Just above and below:
      r = g;
      NCV::trimEquidistantGridUpperEdge
        ( r, std::nextafter( xmax, NC::kInfinity ) );
      REQUIRE( r.npts == t );
      r = g;
      NCV::trimEquidistantGridUpperEdge
        ( r, std::nextafter( xmax, -NC::kInfinity ) );
      REQUIRE( r.npts == ( t > 2 ? t - 1 : 2 ) );
      REQUIRE( r.x1() <= std::nextafter( xmax, -NC::kInfinity )
               || r.npts == 2 );
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
    //individual grids. Those are given by x0+binWidth*i (unfused) except the
    //last one which is x1() (with an explicit std::fma):
    std::vector<NCV::EquidistantGrid> grids = {
      { -0.37, 0.0173, 120 }, { 0.0031, 0.0091, 200 }, { -0.2, 0.0173, 45 },
      { 0.4, 0.0517, 30 } };
    NC::VectD cand;
    for ( auto& g : grids ) {
      for ( std::size_t i = 0; i < g.npts; ++i )
        cand.push_back( i + 1 == g.npts
                        ? g.x1()
                        : g.x0 + rnd( g.binWidth * static_cast<double>(i) ) );
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
  testPWLNarrowToPos();
  testTrimEquidistantGridUpperEdge();
  testMakeCommonGrid();
  testTrimTailByIntegral();
  testExpansionGrids();
  return 0;
}

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
// Tests of generic numerical utilities which are sensitive to floating point //
// contraction (fused multiply-add, FMA). Compilers are free to turn a*b+c    //
// into a fused operation, which gives slightly different results on          //
// platforms with (and without) FMA instructions, breaking reproducibility.   //
//                                                                            //
// Most tests here compare results bit-for-bit against reference              //
// implementations in which the rounding of intermediate products is          //
// enforced by a volatile store (see rnd() below). Thus the tests currently   //
// define the semantics as "unfused". If a site is deliberately changed to    //
// use an explicit std::fma, the corresponding reference here must be changed //
// to match. Some tests instead print checksums of full-precision results,    //
// which must be identical on all platforms.                                  //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCFastConvolve.hh"
#include <iostream>
#include <cstring>
#include <cstdint>

namespace NC=NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  //Forces rounding of an intermediate result to double precision, so the
  //compiler can not fuse the operation producing it with a later one:
  inline double rnd( double x ) { volatile double v = x; return v; }

  //Bit-exact comparison (== on doubles is fine here, we never have NaN):
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

  //FNV-1a hash of the bit patterns of a vector of doubles:
  std::string bitHash( const NC::VectD& v )
  {
    std::uint64_t h = 1469598103934665603ULL;
    for ( double d : v ) {
      std::uint64_t u;
      static_assert( sizeof(u) == sizeof(d), "" );
      std::memcpy( &u, &d, sizeof(u) );
      for ( int i = 0; i < 8; ++i ) {
        h ^= ( u >> (8*i) ) & 0xffULL;
        h *= 1099511628211ULL;
      }
    }
    char buf[32];
    std::snprintf( buf, sizeof(buf), "%016llx",
                   static_cast<unsigned long long>(h) );
    return buf;
  }

  void testLinspace()
  {
    struct Rng { double a, b; unsigned n; };
    const Rng rngs[] = { {0.1,0.7,13}, {-3.3,4.1,101}, {1e-3,1e3,777},
                         {-0.37,0.61,2}, {12345.678,12345.9,5000} };
    for ( auto& r : rngs ) {
      auto res = NC::linspace( r.a, r.b, r.n );
      NC::VectD ref;
      const double interval = ( r.b - r.a ) / ( r.n - 1 );
      for ( unsigned i = 0; i + 1 < r.n; ++i )
        ref.push_back( std::fma( static_cast<double>(i), interval, r.a ) );
      ref.push_back( r.b );
      requireIdentical( res, ref, "linspace" );
    }
    std::cout << "linspace ok" << std::endl;
  }

  void testPowspace()
  {
    struct Rng { double a, b; unsigned n; double p; };
    const Rng rngs[] = { {0.005,0.15,14,4.0}, {0.15,2.5,25,1.5},
                         {2.5,317.3,190,2.0}, {0.1,9.7,50,3.0},
                         {0.1,9.7,50,0.5}, {0.3,8.1,60,2.7} };
    for ( auto& r : rngs ) {
      auto res = NC::powspace( r.a, r.b, r.n, r.p );
      const double step = 1.0 / static_cast<double>( r.n - 1 );
      const double delta = r.b - r.a;
      REQUIRE( res.size() == r.n );
      REQUIRE( res.front() == r.a );
      REQUIRE( res.back() == r.b );
      NC::VectD ref( res.size() );
      ref.front() = r.a;
      ref.back() = r.b;
      for ( unsigned i = 1; i + 1 < r.n; ++i ) {
        const double s = static_cast<double>(i) * step;
        if ( r.p == 2.0 )
          ref[i] = std::fma( delta, s*s, r.a );
        else if ( r.p == 3.0 )
          ref[i] = std::fma( delta*(s*s), s, r.a );
        else if ( r.p == 4.0 )
          ref[i] = std::fma( delta, (s*s)*(s*s), r.a );
        else if ( r.p == 1.5 )
          ref[i] = std::fma( delta, s * std::sqrt(s), r.a );
        else if ( r.p == 0.5 )
          ref[i] = std::fma( delta, std::sqrt(s), r.a );
        else
          ref[i] = std::fma( delta, std::pow( s, r.p ), r.a );
      }
      requireIdentical( res, ref, "powspace" );
    }
    std::cout << "powspace ok" << std::endl;
  }

  void testIntervalPos()
  {
    //Uses an explicit std::fma, so the result must be identical on all
    //platforms and independent of contraction:
    const double as[] = { -1.7, 0.0, 0.3, 12.5, 1e-9 };
    const double bs[] = { 0.0, 0.7, 3.1, 12.6, 1.0000001e-9 };
    for ( double a : as ) {
      for ( double b : bs ) {
        if ( !( b >= a ) )
          continue;
        REQUIRE( NC::intervalPos( a, b, 0.0 ) == a );
        REQUIRE( NC::intervalPos( a, b, 1.0 ) == b );
        REQUIRE( NC::intervalPos01( a, b, 0.0 ) == a );
        REQUIRE( NC::intervalPos01( a, b, 1.0 ) == b );
        for ( int i = 0; i <= 100; ++i ) {
          const double r = i * 0.01;
          const double ref = std::fma( r, b, rnd( (1.0-r)*a ) );
          REQUIRE( NC::intervalPos( a, b, r ) == ref );
          const double c = NC::intervalPos01( a, b, r );
          REQUIRE( c >= a && c <= b );
        }
      }
    }
    std::cout << "intervalPos ok" << std::endl;
  }

  void testNclerp()
  {
    //Uses an explicit std::fma, so the result must be identical on all
    //platforms and independent of contraction. Unlike intervalPos, a and b
    //need not be ordered:
    const double as[] = { -1.7, 0.0, 0.3, 12.5, 1e-9, 5.0 };
    const double bs[] = { 0.0, 0.7, 3.1, 12.6, 1.0000001e-9, -5.0 };
    for ( double a : as ) {
      for ( double b : bs ) {
        //Exact at t=0 is provable (0*(b-a) is exactly 0.0 for any finite
        //b-a, so the fma reduces to exactly a); NOT asserting exactness at
        //t=1, which (unlike for intervalPos's fma(t,b,(1-t)*a) form) is not
        //generally guaranteed for fma(t,b-a,a):
        REQUIRE( NC::nclerp( a, b, 0.0 ) == a );
        for ( int i = -20; i <= 120; ++i ) {
          const double t = i * 0.01;
          //b-a is a plain subtraction, not itself at risk of contraction
          //(only a*b+c-shaped expressions are), so no rnd() needed here:
          const double ref = std::fma( t, b - a, a );
          REQUIRE( NC::nclerp( a, b, t ) == ref );
        }
      }
    }
    std::cout << "nclerp ok" << std::endl;
  }

  void testStableSums()
  {
    //The compensation in Kahan/Neumaier summation depends on the addends being
    //rounded doubles. If the product forming an addend is fused into the
    //addition inside add(), the compensation is no longer correct. So adding
    //a product directly must give the same as adding the rounded product:
    NC::StableSum n1, n2, p1, p2;
    NC::StableSumKahan k1, k2;
    for ( int i = 0; i < 5000; ++i ) {
      const double a = 1.0 + i * 1.0e-3;
      const double b = 0.3 + i * 7.1e-4;
      const double c = ( i % 3 ? 1.0 : -1.0 ) * a * 0.37;
      n1.add( a*b );
      n2.add( rnd( a*b ) );
      n1.add( c*b );
      n2.add( rnd( c*b ) );
      p1.addPosVal( a*b );
      p2.addPosVal( rnd( a*b ) );
      k1.add( a*b );
      k2.add( rnd( a*b ) );
      REQUIRE( n1.sum() == n2.sum() );
      REQUIRE( p1.sum() == p2.sum() );
      REQUIRE( k1.sum() == k2.sum() );
    }
    std::cout << "stable sums ok" << std::endl;
  }

  void testFastConvolve()
  {
    //The FFT is designed to be bit-for-bit reproducible across platforms (exact
    //table of phase factors, no calls to trigonometric functions), and so a
    //hash of the bit patterns of the results must be the same everywhere. This
    //detects fused multiply-adds in the FFT butterfly operations.
    NC::VectD a1, a2;
    for ( int i = 0; i < 100; ++i )
      a1.push_back( 1.0 / std::fma( 0.01, static_cast<double>(i*i), 1.0 ) );
    for ( int i = 0; i < 77; ++i )
      a2.push_back( std::fma( 0.1, static_cast<double>( i % 7 + 1 ),
                              1.0 / ( 3.0 + i ) ) );
    NC::FastConvolve fc;
    NC::VectD y, yl;
    fc.convolve( a1, a2, y, 0.0137 );
    fc.convolveLegacy( a1, a2, yl, 0.0137 );
    REQUIRE( y.size() == a1.size() + a2.size() - 1 );
    REQUIRE( yl.size() == y.size() );
    std::cout << "convolve: y[0]=" << NC::fmt(y.front(),"%.10g")
              << " y[50]=" << NC::fmt(y[50],"%.10g")
              << " y[175]=" << NC::fmt(y.back(),"%.10g") << std::endl;
    std::cout << "convolve bit-hash:       " << bitHash(y) << std::endl;
    std::cout << "convolveLegacy bit-hash: " << bitHash(yl) << std::endl;

    //Larger, needs a bigger table of phase factors and more butterfly stages:
    NC::VectD b1, b2;
    for ( int i = 0; i < 1500; ++i )
      b1.push_back( 1.0 / std::fma( 1e-4, static_cast<double>(i*i), 1.0 ) );
    for ( int i = 0; i < 1101; ++i )
      b2.push_back( std::fma( 0.05, static_cast<double>( i % 11 ),
                              1.0 / ( 7.0 + i ) ) );
    fc.convolve( b1, b2, y, 0.0031 );
    std::cout << "convolve (large) bit-hash: " << bitHash(y) << std::endl;
  }

  void testReducePts()
  {
    //The result of the greedy point-reduction depends on tiny differences
    //between importance values of the candidate points, and thus is very
    //sensitive to contraction. The data is deliberately (almost) symmetric,
    //making many candidates nearly tied.
    NC::VectD x, y;
    for ( int i = -100; i <= 100; ++i ) {
      const double xi = i * 0.05;
      x.push_back( xi );
      y.push_back( 1.0 / ( 1.0 + xi*xi )
                   + 0.3 / ( 1.0 + 9.0*(xi-2.0)*(xi-2.0) ) );
    }
    for ( std::size_t target : { 20, 40, 90 } ) {
      auto res = NC::reducePtsInDistribution( x, y, target );
      REQUIRE( res.first.size() == target );
      std::cout << "reducePts(" << target << ") kept indices:";
      for ( double xr : res.first )
        std::cout << " " << static_cast<int>( std::lround( xr/0.05 ) + 100 );
      std::cout << std::endl;
    }
  }

}

int main()
{
  testLinspace();
  testPowspace();
  testIntervalPos();
  testNclerp();
  testStableSums();
  testFastConvolve();
  testReducePts();
  return 0;
}

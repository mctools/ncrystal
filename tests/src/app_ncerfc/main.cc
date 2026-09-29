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

//Validation of the portable ncerf/ncerfc functions: ULP-level agreement
//with high-precision (mpmath-generated, hardwired) reference values,
//wide sweeps against the platform's std::erf/std::erfc, structural
//properties (symmetry, monotonicity across the internal region seams),
//and -- via the exact bit patterns printed into the reference log -- a
//CI-enforced check that results are bit-identical on every platform.
//Timing versus libm is available via NCRYSTAL_NCERFC_TIMING=1 (not part
//of the reference log).

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <chrono>

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  std::uint64_t dblbits( double x )
  {
    std::uint64_t r;
    static_assert( sizeof(r) == sizeof(x), "" );
    std::memcpy( &r, &x, sizeof(r) );
    return r;
  }

  double ulpdiff( double a, double b )
  {
    //Distance in units of ULPs of b (b=reference, assumed finite and
    //normal; a and b assumed same sign region):
    if ( a == b )
      return 0.0;
    const double u = NC::ncabs( std::nextafter( b, a ) - b );
    return NC::ncabs( a - b ) / u;
  }

  void test_mpmath_refvals()
  {
    std::printf("test_mpmath_refvals:\n");
    struct Ref { double x, erfc, erf; };
    //Reference values generated with mpmath at 40 significant digits
    //(rounded here to 22, i.e. exact at double precision). NB: they are
    //evaluated at the exact binary double closest to each x literal --
    //evaluating at the exact decimal instead can differ by >10 ULP for
    //non-representable x, since d(ln(erfc))/dx ~ -2x:
    const Ref refvals[] = {
      { 0.0         , 1.000000000000000000000, 0.0 },
      { 1e-300      , 1.000000000000000000000, 1.128379167095512602172e-300 },
      { -1e-300     , 1.000000000000000000000, -1.128379167095512602172e-300 },
      { 1e-20       , 0.9999999999999999999887, 1.128379167095512512008e-20 },
      { -1e-20      , 1.000000000000000000011, -1.128379167095512512008e-20 },
      { 1.11e-16    , 0.9999999999999998747499, 1.252500875476019003168e-16 },
      { -1.11e-16   , 1.000000000000000125250, -1.252500875476019003168e-16 },
      { 1e-08       , 0.9999999887162083290449, 1.128379167095512559892e-8 },
      { -1e-08      , 1.000000011283791670955, -1.128379167095512559892e-8 },
      { 0.0001      , 0.9998871620836665751251, 0.0001128379163334248748935 },
      { -0.0001     , 1.000112837916333424875, -0.0001128379163334248748935 },
      { 0.01        , 0.9887165844441503828492, 0.01128341555584961715078 },
      { -0.01       , 1.011283415555849617151, -0.01128341555584961715078 },
      { 0.1         , 0.8875370839817151015953, 0.1124629160182848984047 },
      { -0.1        , 1.112462916018284898405, -0.1124629160182848984047 },
      { 0.25        , 0.7236736098317630670149, 0.2763263901682369329851 },
      { -0.25       , 1.276326390168236932985, -0.2763263901682369329851 },
      { 0.4         , 0.5716076449533315235459, 0.4283923550466684764541 },
      { -0.4        , 1.428392355046668476454, -0.4283923550466684764541 },
      { 0.46875     , 0.5073865267820620084118, 0.4926134732179379915882 },
      { -0.46875    , 1.492613473217937991588, -0.4926134732179379915882 },
      { 0.469       , 0.5071601050374203271526, 0.4928398949625796728474 },
      { -0.469      , 1.492839894962579672847, -0.4928398949625796728474 },
      { 0.5         , 0.4795001221869534623173, 0.5204998778130465376827 },
      { -0.5        , 1.520499877813046537683, -0.5204998778130465376827 },
      { 0.75        , 0.2888443663464848684011, 0.7111556336535151315989 },
      { -0.75       , 1.711155633653515131599, -0.7111556336535151315989 },
      { 1.0         , 0.1572992070502851306588, 0.8427007929497148693412 },
      { -1.0        , 1.842700792949714869341, -0.8427007929497148693412 },
      { 1.5         , 0.03389485352468927293302, 0.9661051464753107270670 },
      { -1.5        , 1.966105146475310727067, -0.9661051464753107270670 },
      { 2.0         , 0.004677734981047265837931, 0.9953222650189527341621 },
      { -2.0        , 1.995322265018952734162, -0.9953222650189527341621 },
      { 2.5         , 0.0004069520174449589395642, 0.9995930479825550410604 },
      { -2.5        , 1.999593047982555041060, -0.9995930479825550410604 },
      { 3.0         , 0.00002209049699858544137278, 0.9999779095030014145586 },
      { -3.0        , 1.999977909503001414559, -0.9999779095030014145586 },
      { 3.5         , 7.430983723414127455237e-7, 0.9999992569016276585873 },
      { -3.5        , 1.999999256901627658587, -0.9999992569016276585873 },
      { 3.999       , 1.554474949099499153278e-8, 0.9999999844552505090050 },
      { -3.999      , 1.999999984455250509005, -0.9999999844552505090050 },
      { 4.0         , 1.541725790028001885216e-8, 0.9999999845827420997200 },
      { -4.0        , 1.999999984582742099720, -0.9999999845827420997200 },
      { 4.001       , 1.529078217324873151904e-8, 0.9999999847092178267513 },
      { -4.001      , 1.999999984709217826751, -0.9999999847092178267513 },
      { 5.0         , 1.537459794428034850188e-12, 0.9999999999984625402056 },
      { -5.0        , 1.999999999998462540206, -0.9999999999984625402056 },
      { 6.5         , 3.842148327120647469876e-20, 0.9999999999999999999616 },
      { -6.5        , 1.999999999999999999962, -0.9999999999999999999616 },
      { 8.0         , 1.122429717298292707997e-29, 1.000000000000000000000 },
      { -8.0        , 2.000000000000000000000, -1.000000000000000000000 },
      { 10.0        , 2.088487583762544757001e-45, 1.000000000000000000000 },
      { -10.0       , 2.000000000000000000000, -1.000000000000000000000 },
      { 13.0        , 1.739557315466724521804e-75, 1.000000000000000000000 },
      { -13.0       , 2.000000000000000000000, -1.000000000000000000000 },
      { 17.0        , 1.021228015094260881146e-127, 1.000000000000000000000 },
      { -17.0       , 2.000000000000000000000, -1.000000000000000000000 },
      { 21.0        , 8.032453871022455669021e-194, 1.000000000000000000000 },
      { -21.0       , 2.000000000000000000000, -1.000000000000000000000 },
      { 25.0        , 8.300172571196522752044e-274, 1.000000000000000000000 },
      { -25.0       , 2.000000000000000000000, -1.000000000000000000000 },
      { 26.0        , 5.663192408856142846476e-296, 1.000000000000000000000 },
      { -26.0       , 2.000000000000000000000, -1.000000000000000000000 },
      { 26.5        , 2.210907664263734275929e-307, 1.000000000000000000000 },
      { -26.5       , 2.000000000000000000000, -1.000000000000000000000 },
    };
    double worst_erfc(0.0), worst_erf(0.0);
    for ( auto& r : refvals ) {
      const double u1 = ulpdiff( NC::ncerfc(r.x), r.erfc );
      const double u2 = ulpdiff( NC::ncerf(r.x), r.erf );
      worst_erfc = NC::ncmax( worst_erfc, u1 );
      worst_erf = NC::ncmax( worst_erf, u2 );
      REQUIRE( u1 <= 4.0 );
      REQUIRE( u2 <= 4.0 );
    }
    std::printf("  %u ref points, worst ulp-diff: erfc %.2f, erf %.2f:"
                " OK\n", unsigned(sizeof(refvals)/sizeof(Ref)),
                worst_erfc, worst_erf );
    //Beyond the cutoff the result is flushed to exactly zero:
    for ( double x : { 26.543, 27.0, 100.0, 1e300 } ) {
      REQUIRE( NC::ncerfc(x) == 0.0 );
      REQUIRE( NC::ncerfc(-x) == 2.0 );
      REQUIRE( NC::ncerf(x) == 1.0 && NC::ncerf(-x) == -1.0 );
    }
    REQUIRE( std::isnan( NC::ncerfc( 0.0/0.0 ) ) );
    REQUIRE( std::isnan( NC::ncerf( 0.0/0.0 ) ) );
    std::printf("  cutoff and NaN behaviour: OK\n");
  }

  void test_vs_libm()
  {
    //Wide sweeps against the platform libm. NB: the bar must leave room
    //for the libm's own error (glibc's erf/erfc are correctly rounded
    //nowadays, other platforms are typically within 1-2 ULP):
    std::printf("test_vs_libm:\n");
    double worst(0.0);
    unsigned n(0);
    auto testpt = [&worst,&n]( double x )
    {
      const double u1 = ulpdiff( NC::ncerfc(x), std::erfc(x) );
      worst = NC::ncmax( worst, u1 );
      REQUIRE( u1 <= 8.0 );
      const double referf = std::erf(x);
      if ( NC::ncabs(referf) > 1e-290 ) {
        const double u2 = ulpdiff( NC::ncerf(x), referf );
        worst = NC::ncmax( worst, u2 );
        REQUIRE( u2 <= 8.0 );
      }
      ++n;
    };
    for ( int i = 0; i <= 1200000; ++i ) {
      const double x = -6.0 + i * 1e-5;
      testpt( x );
    }
    for ( int i = 0; i <= 100000; ++i )
      testpt( 4.0 + i * 2.2e-4 );//[4,26]
    for ( int i = 0; i <= 10000; ++i ) {
      const double x = 1e-300 * std::pow( 10.0, i * 0.0299 );//to ~1e-1
      testpt( x );
      testpt( -x );
    }
    //NB: the observed worst ulp-distance depends on the platform's own
    //libm, so it must not enter the (byte-compared) reference log:
    REQUIRE( worst >= 0.0 );
    std::printf("  %u points swept vs std::erf/std::erfc, all within"
                " 8 ulp: OK\n", n );
  }

  void test_properties()
  {
    std::printf("test_properties:\n");
    //Exact symmetries (by construction, but pin them):
    for ( int i = 0; i <= 10000; ++i ) {
      const double x = -27.0 + i * 0.0054;
      REQUIRE( NC::ncerf(-x) == -NC::ncerf(x) );
      const double s = NC::ncerfc(x) + NC::ncerfc(-x);
      REQUIRE( NC::ncabs( s - 2.0 ) < 3e-16 );
    }
    //Monotonicity of erfc (non-increasing), in particular across the
    //internal region seams at 0.46875, 4.0 and 26.543:
    double prev = 2.0;
    for ( int i = 0; i <= 2000000; ++i ) {
      const double x = -27.0 + i * 2.7e-5;
      const double v = NC::ncerfc(x);
      REQUIRE( v <= prev );
      prev = v;
    }
    //erf/erfc consistency where both are well-scaled:
    for ( int i = 0; i <= 100000; ++i ) {
      const double x = i * 6e-5;
      const double d = NC::ncerf(x) + NC::ncerfc(x) - 1.0;
      REQUIRE( NC::ncabs(d) < 3e-16 );
    }
    std::printf("  symmetry, monotonicity and erf+erfc==1: OK\n");
  }

  void print_bitpatterns()
  {
    //The exact bit patterns below are part of the reference log, so any
    //platform producing even a 1-ULP different value will fail the
    //byte-exact CI log comparison -- this is the actual enforcement of
    //the "bit-identical on all platforms" contract:
    std::printf("bit-identity table (uint64 patterns):\n");
    const double xs[] = { 1e-9, 0.001, 0.1, 0.3, 0.46875, 0.469, 0.7,
                          1.0, 1.5, 2.25, 3.0, 3.75, 4.0, 4.5, 6.0,
                          9.0, 14.0, 20.0, 26.0, 26.5 };
    for ( double x : xs ) {
      static_assert( sizeof(unsigned long long)>=8, "" );
      std::printf( "  x=%-8g : erfc(x) %016llx erfc(-x) %016llx"
                   " erf(x) %016llx\n", x,
                   (unsigned long long)dblbits( NC::ncerfc(x) ),
                   (unsigned long long)dblbits( NC::ncerfc(-x) ),
                   (unsigned long long)dblbits( NC::ncerf(x) ) );
    }
  }

  void timing()
  {
    if ( !NC::ncgetenv_bool("NCERFC_TIMING") )
      return;
    struct BRange { const char* lbl; double lo, hi; };
    for ( auto br : { BRange{"x in [0,0.46]",0.0,0.46},
                      BRange{"x in [0.47,4]",0.47,4.0},
                      BRange{"x in [4,26]",4.0,26.0} } ) {
      for ( int mode = 0; mode < 2; ++mode ) {
        const int nn = 10000000;
        const double dx = ( br.hi - br.lo ) / nn;
        double x = br.lo;
        NC::StableSum sum;
        auto t0 = std::chrono::steady_clock::now();
        for ( int i = 0; i < nn; ++i ) {
          sum.add( mode ? NC::ncerfc(x) : std::erfc(x) );
          x += dx;
        }
        auto dt = std::chrono::duration<double,std::nano>(
          std::chrono::steady_clock::now() - t0 ).count();
        std::printf("  TIMING %-14s %s : %5.1f ns/call (checksum"
                    " %g)\n", br.lbl, mode ? "ncerfc   " : "std::erfc",
                    dt / nn, sum.sum() );
      }
    }
  }

}

int main()
{
  test_mpmath_refvals();
  test_vs_libm();
  test_properties();
  print_bitpatterns();
  timing();
  std::printf("All tests passed.\n");
  return 0;
}

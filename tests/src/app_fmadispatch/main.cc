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
// Self-test for NCRYSTAL_FMADISPATCH_ATTR (see docs/devel_fma_attribute.md). //
// Applies it to a FastConvolve-shaped hot loop (matching the jr=a*c-b*d,     //
// ji=a*d+b*c spectral-multiply step in NCFastConvolve.cc), written with      //
// explicit std::fma throughout as required, and checks three things:         //
//                                                                            //
//  1. Unconditionally: the dispatched function gives bit-for-bit identical   //
//     results to the same code, never multiversioned -- the attribute must   //
//     never change what is computed, only how fast it runs.                  //
//  2. Only when NCRYSTAL_FMADISPATCH_ENABLED on this build: the dispatched   //
//     function is actually measurably faster. This is the point of the       //
//     test: a future regression (e.g. a toolchain/linking change that        //
//     silently stops the ifunc resolver from engaging, quietly falling back  //
//     to the slow generic path) must be caught here, rather than silently    //
//     running without its benefit forever after.                             //
//  3. On platforms/compilers we are confident about, NCRYSTAL_FMADISPATCH_   //
//     ENABLED actually has the value it should. This is an INDEPENDENT       //
//     cross-check, using a crude compiler/platform-predefine guess of our    //
//     own (deliberately not reusing any CMake/simplebuild logic), so a bug   //
//     that made the real detection (the CMake probe, or simplebuild's        //
//     assumption -- see ncrystal_fmadispatch.cmake and sbgen/main.py)        //
//     silently disagree with what the platform obviously supports gets       //
//     caught here too, instead of us quietly running blind forever.          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/core/NCException.hh"
#include <iostream>
#include <chrono>

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

////////////////////////////////////////////////////////////////////////////
// Our own, independent guess (see point 3 above) at whether this platform //
// should have the technique enabled or disabled. Intentionally simple and //
// NOT reused from anywhere else in the project: its entire value is in    //
// being a second opinion. Confident cases only -- anything else (32-bit   //
// x86, other OSes, in-between compiler versions, ...) is left             //
// unclassified and simply not checked:                                    //
////////////////////////////////////////////////////////////////////////////

#if ( defined(__x86_64__) || defined(__i386__) )                        \
  && ( defined(__linux__) || defined(__APPLE__) )                       \
  && ( ( defined(__clang__) && __clang_major__ >= 14 )                  \
       || ( defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 6 ) )
#  define NCTEST_FMADISPATCH_EXPECT_ENABLED 1
#else
#  define NCTEST_FMADISPATCH_EXPECT_ENABLED 0
#endif

#if defined(__aarch64__) || defined(__arm__)    \
  || defined(_M_ARM64) || defined(_M_ARM)       \
  || defined(_WIN32) || defined(_WIN64)
#  define NCTEST_FMADISPATCH_EXPECT_DISABLED 1
#else
#  define NCTEST_FMADISPATCH_EXPECT_DISABLED 0
#endif

//These two must never both claim a confident classification at once (that
//would mean the two #if conditions above overlap, a bug in this test
//itself). Both being 0 (the unclassified, "don't know" case) is fine and
//expected by design -- e.g. an old-but-not-ancient x86/Linux compiler, or
//an architecture that is neither x86 nor ARM:
static_assert( !( NCTEST_FMADISPATCH_EXPECT_ENABLED && NCTEST_FMADISPATCH_EXPECT_DISABLED ),
              "EXPECT_ENABLED and EXPECT_DISABLED must not both be true" );

//These are compile-time constants, so check them with static_assert rather
//than a runtime REQUIRE -- a wrong value fails the build outright:
#if NCTEST_FMADISPATCH_EXPECT_ENABLED
static_assert( NCRYSTAL_FMADISPATCH_ENABLED == 1,
              "NCRYSTAL_FMADISPATCH_ENABLED should be 1 on this platform/compiler"
              " (x86, GCC/Clang, Linux/macOS) -- see doc/devel_fma_attribute.md" );
#elif NCTEST_FMADISPATCH_EXPECT_DISABLED
static_assert( NCRYSTAL_FMADISPATCH_ENABLED == 0,
              "NCRYSTAL_FMADISPATCH_ENABLED should be 0 on this platform"
              " (ARM or Windows) -- see doc/devel_fma_attribute.md" );
#endif

namespace {

  NCRYSTAL_FMADISPATCH_ATTR
  void cmulDispatch( double* re, double* im,
                     const double* cre, const double* cim, std::size_t n )
  {
    for ( auto i : NC::ncrange(n) ) {
      double a = re[i], b = im[i], c = cre[i], d = cim[i];
      re[i] = std::fma( a, c, -(b*d) );
      im[i] = std::fma( a, d, b*c );
    }
  }

  //Identical code, never multiversioned, so it always runs exactly as
  //written, on every platform:
  void cmulReference( double* re, double* im,
                      const double* cre, const double* cim, std::size_t n )
  {
    for ( auto i : NC::ncrange(n) ) {
      double a = re[i], b = im[i], c = cre[i], d = cim[i];
      re[i] = std::fma( a, c, -(b*d) );
      im[i] = std::fma( a, d, b*c );
    }
  }

#if NCRYSTAL_FMADISPATCH_ENABLED
  //Best-of-ntrials timing, in ns/element. Guarded by the same #if as its
  //only use below: when dispatch is disabled, this would otherwise be an
  //unused function, which fails the project's -Werror build:
  double timeIt( void(*f)(double*,double*,const double*,const double*,std::size_t),
                 const NC::VectD& re0, const NC::VectD& im0,
                 const NC::VectD& cre, const NC::VectD& cim,
                 int ntrials, int nrep )
  {
    const std::size_t n = re0.size();
    double best = 1e300;
    for ( auto trial : NC::ncrange(ntrials) ) {
      (void)trial;
      NC::VectD re(re0), im(im0);
      const auto t0 = std::chrono::steady_clock::now();
      for ( auto r : NC::ncrange(nrep) ) {
        (void)r;
        f( re.data(), im.data(), cre.data(), cim.data(), n );
      }
      const auto t1 = std::chrono::steady_clock::now();
      volatile double sink = re[0] + im[0];
      (void)sink;
      const double ns = std::chrono::duration<double,std::nano>(t1-t0).count();
      best = std::min( best, ns / ( static_cast<double>(nrep) * static_cast<double>(n) ) );
    }
    return best;
  }
#endif

}

int main()
{
  //3. Cross-check NCRYSTAL_FMADISPATCH_ENABLED against our own, independent
  //guess, on the platforms/compilers we are confident about. Already
  //enforced above by static_assert (a wrong value fails the build, not this
  //line); this just reports which case applied:
#if NCTEST_FMADISPATCH_EXPECT_ENABLED
  std::cout << "Sanity check: this looks like a platform/compiler where"
    " NCRYSTAL_FMADISPATCH should be ENABLED (x86, GCC/Clang, Linux/macOS),"
    " and it is." << std::endl;
#elif NCTEST_FMADISPATCH_EXPECT_DISABLED
  std::cout << "Sanity check: this looks like a platform where"
    " NCRYSTAL_FMADISPATCH should be DISABLED (ARM or Windows), and it is."
            << std::endl;
#else
  std::cout << "Sanity check: platform/compiler not confidently classified"
    " by this test; skipping the enabled/disabled expectation check."
            << std::endl;
#endif

  const std::size_t n = 32768;
  NC::VectD re0(n), im0(n), cre(n), cim(n);
  for ( auto i : NC::ncrange(n) ) {
    const double x = static_cast<double>(i);
    re0[i] = 1.0/(1.0+0.001*x);
    im0[i] = 0.3*std::sin(0.01*x);
    cre[i] = std::cos(0.02*x);
    cim[i] = std::sin(0.02*x);
  }

  //1. Correctness, unconditionally:
  NC::VectD reA(re0), imA(im0), reB(re0), imB(im0);
  cmulDispatch( reA.data(), imA.data(), cre.data(), cim.data(), n );
  cmulReference( reB.data(), imB.data(), cre.data(), cim.data(), n );
  const bool correct = ( reA == reB && imA == imB );
  std::cout << "correctness check (dispatched == reference, bit for bit): "
            << ( correct ? "PASS" : "FAIL" ) << std::endl;
  REQUIRE( correct );

#if NCRYSTAL_FMADISPATCH_ENABLED
  //2. This build claims dispatch is available: verify it is not merely
  //compiling but actually being used. The threshold is set well below the
  //~5-15x typically measured (examples/fmadispatch_investigate), to stay
  //robust on a loaded/virtualised CI machine while still failing hard on
  //"no dispatch happened at all" (which would show up as ~1x):
  const int ntrials = 8, nrep = 50;
  const double tDispatch = timeIt( cmulDispatch, re0, im0, cre, cim, ntrials, nrep );
  const double tReference = timeIt( cmulReference, re0, im0, cre, cim, ntrials, nrep );
  const double speedup = tReference / tDispatch;
  std::cout << "cmulDispatch  : " << tDispatch << " ns/element" << std::endl;
  std::cout << "cmulReference : " << tReference << " ns/element" << std::endl;
#ifdef NDEBUG
  std::cout << "speedup (reference/dispatch): " << speedup
            << "x (required: >1.5x)" << std::endl;
  REQUIRE( speedup > 1.5 );
#else
  //Unoptimised (e.g. simplebuild debug cache) builds do not give
  //target_clones clones the optimisation needed to show their benefit;
  //not enforcing a minimum here:
  std::cout << "speedup (reference/dispatch): " << speedup
            << "x (unoptimised/debug build: not enforcing a minimum)"
            << std::endl;
#endif
#else
  std::cout << "NCRYSTAL_FMADISPATCH_ENABLED=0 on this build:"
    " skipping speedup check." << std::endl;
#endif

  std::cout << "app_fmadispatch: all checks passed" << std::endl;
  return 0;
}

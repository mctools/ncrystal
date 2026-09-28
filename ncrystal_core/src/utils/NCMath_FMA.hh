#ifndef NCrystal_Math_FMA_hh
#define NCrystal_Math_FMA_hh

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

#include "NCrystal/internal/utils/NCFMADispatch.hh"
#include "NCrystal/internal/utils/NCMath.hh"

//The ordinary and MSVC/x86-only Windows-fast variants of NCMath.cc's
//expm1_reducedarg_taylor14/detail_stable_expm1/detail_stable_exp -- see
//NCFMADispatch.hh for the macros used below and the full mechanism/
//rationale, and NCMath_WINFMA.cc for the other half.
//
//expm1_reducedarg_taylor14 is called only from within this same file/TU
//(by detail_stable_expm1), so unlike the other two below it needs no
//WINFMA_DECLARE/FORWARD pair of its own: it is written once, unconditionally,
//as a plain NCRYSTAL_FMADISPATCH_ATTR anonymous-namespace helper, and simply
//gets compiled under whatever flags its *caller's* translation unit uses --
//once NCMath_WINFMA.cc dispatches into it via its own /arch:AVX2 build,
//this helper (also compiled into that same TU) already benefits, with no
//separate twin needed.
//
//detail_stable_expm1/detail_stable_exp use NCRYSTAL_FMADISPATCH_DECLARATOR_C
//(not the plain DECLARATOR): both are already extern "C" +
//NCRYSTAL_APPLY_C_NAMESPACE-wrapped in their ordinary form too, for the
//unrelated (Apple/Mach-O) rule-5 reason documented at their original
//definition site -- see docs/devel_fma_attribute.md. Their mutual/
//recursive calls (by detail_stable_exp, of itself and of
//detail_stable_expm1) are written using the plain, un-suffixed extern "C"
//name in both branches (shared, unconditional body text) -- from the
///arch:AVX2-compiled twin, such a call still reaches the *ordinary* entry
//point (the only place that exact symbol is ever defined), but that entry
//point's own WINFMA_FORWARD check immediately re-dispatches back to the
//fast twin, so the call still ends up on the fast path, just via one extra
//indirection.

namespace NCRYSTAL_NAMESPACE {
  namespace {

    //expm1(r)=e^r-1 via a 14th order Taylor expansion (terms r^1/1! to
    //r^14/14!), Horner's method with std::fma throughout for reproducibility
    //(same technique as safe_xcothx in NCVDOSEval.cc). Only accurate for
    //the reduced range |r|<=ln(2)/2~0.3466 produced by the Cody-Waite
    //reduction in detail_stable_expm1 below: truncation error there is
    //~4e-18 relative (first omitted term, r^15/15!), comfortably below
    //double precision. Coefficients are exact rationals rounded to the
    //nearest double, so identical on every platform by construction:
    NCRYSTAL_FMADISPATCH_ATTR
    double expm1_reducedarg_taylor14( double r )
    {
      constexpr double c1 = 1.0;
      constexpr double c2 = 1.0/2;
      constexpr double c3 = 1.0/6;
      constexpr double c4 = 1.0/24;
      constexpr double c5 = 1.0/120;
      constexpr double c6 = 1.0/720;
      constexpr double c7 = 1.0/5040;
      constexpr double c8 = 1.0/40320;
      constexpr double c9 = 1.0/362880;
      constexpr double c10 = 1.0/3628800;
      constexpr double c11 = 1.0/39916800;
      constexpr double c12 = 1.0/479001600;
      constexpr double c13 = 1.0/6227020800.0;
      constexpr double c14 = 1.0/87178291200.0;
      double p = c14;
      p = std::fma( r, p, c13 );
      p = std::fma( r, p, c12 );
      p = std::fma( r, p, c11 );
      p = std::fma( r, p, c10 );
      p = std::fma( r, p, c9 );
      p = std::fma( r, p, c8 );
      p = std::fma( r, p, c7 );
      p = std::fma( r, p, c6 );
      p = std::fma( r, p, c5 );
      p = std::fma( r, p, c4 );
      p = std::fma( r, p, c3 );
      p = std::fma( r, p, c2 );
      p = std::fma( r, p, c1 );
      return r*p;
    }

  }

  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, stable_expm1,
    ( double x )
  )
  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, stable_exp,
    ( double x )
  )

  //Forward declarations of the plain (non-_winfma) extern "C" names:
  //detail_stable_exp's body below calls itself recursively and calls
  //detail_stable_expm1, always by this same, un-suffixed name (shared,
  //unconditional body text -- see NCFMADispatch.hh). In the ordinary
  //branch that name is simply its own (about to be defined) declarator,
  //but in the /arch:AVX2-compiled NCMath_WINFMA.cc translation unit only
  //the *_winfma-suffixed twins are actually defined there, so without an
  //explicit declaration here, that TU would never see a declaration of
  //the plain names at all -- a compile error (C++ has no implicit
  //declarations), not just an extra indirection. Harmless in the ordinary
  //branch too (an extern "C" declaration ahead of a matching definition
  //is always fine):
  extern "C" double NCRYSTAL_APPLY_C_NAMESPACE(detail_stable_expm1)( double x );
  extern "C" double NCRYSTAL_APPLY_C_NAMESPACE(detail_stable_exp)( double x );

  //Cody-Waite argument reduction constants for splitting x = n*ln(2) + r
  //with |r|<=ln(2)/2: ln2_hi carries only the top ~32 bits of ln(2) (all
  //lower bits exactly 0), so n*ln2_hi is exact for any |n| relevant here,
  //and the two-step fma reduction below recovers r to near full double
  //precision even though x and n*ln(2) can be far larger than r itself.
  //Standard constants (as used by e.g. fdlibm's exp()):
  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_stable_expm1,stable_expm1)
  ( double x )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(stable_expm1,(x));
    constexpr double ln2_hi = 6.93147180369123816490e-01;
    constexpr double ln2_lo = 1.90821492927058770002e-10;
    constexpr double invln2 = 1.44269504088896338700e+00;
    //Avoids std::exp/std::log/std::expm1/std::log1p entirely (unlike the
    //Newton-Raphson-on-libm approach this replaced), so the result is
    //identical on every platform, not just very close. An earlier version of
    //this function range-reduced by repeated exact halving down to a small
    //threshold, then reconstructed via the doubling identity
    //expm1(2y)=t*(t+2) applied once per halving -- correct, and still used
    //in stable_sinh/stable_tanh, but each doubling step roughly *doubles*
    //the relative error of the running result, so ~15 doublings (needed to
    //cover the physically relevant range) compounds a tiny base error into
    //several ULP, failing app_stablemathfct's 8-ULP bar for |x|>~10.
    //Replaced with the standard single-shot Cody-Waite reduction instead
    //(x=n*ln(2)+r, |r|<=ln(2)/2, both std::fma steps below), which has no
    //such compounding since there is only one reduction step:
    if ( !std::isfinite(x) )
      return ncisnan(x) ? x : ( x < 0.0 ? -1.0 : x );//nan; -inf->-1; +inf->+inf
    //Short-circuit comfortably before exp(x) itself would overflow/underflow
    //to inf/-1 to full double precision anyway: besides being pointless work,
    //this keeps n below safely within int range for the static_cast below
    //(without it, e.g. x=1e10 would overflow that cast -- undefined
    //behaviour -- long before exp(x) itself would legitimately overflow):
    if ( x >= 710.0 )
      return kInfinity;
    if ( x <= -40.0 )
      return -1.0;
    const double n = std::round( x * invln2 );
    double r = std::fma( -n, ln2_hi, x );
    r = std::fma( -n, ln2_lo, r );
    const double expm1_r = expm1_reducedarg_taylor14( r );
    //expm1(x) = 2^n*exp(r) - 1 = 2^n*(1+expm1(r)) - 1 = 2^n*expm1(r) +
    //(2^n-1). std::ldexp is exact (pure exponent adjustment, no rounding),
    //and 2^n-1 is likewise exact (subtracting an exactly representable value
    //no larger than 1 in magnitude from 1 never needs more mantissa bits than
    //already available), so the only rounding in the whole reconstruction is
    //the single final std::fma:
    const double pow2n = std::ldexp( 1.0, static_cast<int>(n) );
    return std::fma( pow2n, expm1_r, pow2n - 1.0 );
  }

  //Calls the *detail* expm1 directly (not the public wrapper), to stay on
  //the C-linkage, dispatched fast path throughout, including the recursive
  //self-call for negative x:
  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_stable_exp,stable_exp)
  ( double x )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(stable_exp,(x));
    //exp(x)=1+expm1(x): unlike the "obvious" exp(x)-1 (which cancels
    //catastrophically for small x, the reason expm1 exists in the first
    //place), this direction has no cancellation at all for x>=0 -- it is a
    //plain addition of two same-signed quantities, introducing no extra
    //error beyond the single final rounding.
    //
    //But this identity alone is unsafe for negative x: once exp(x) drops
    //below ulp(1)~1.1e-16, expm1(x) correctly rounds to exactly -1.0 (nothing
    //wrong with expm1 there -- that IS its correctly-rounded value), yet
    //"1.0+(-1.0)" is then exactly 0.0 regardless of how much smaller than
    //that the true exp(x) still is (confirmed empirically: 1+expm1(-199.3)
    //rounds to 0.0 vs the true ~1.4e-87 -- a total loss of precision, and
    //this holds for any implementation of expm1, not just this one).
    //Sidestepped with the standard reciprocal trick already used by
    //exp_approx/exp_negarg_approx in NCMath.cc: for negative x, compute
    //exp(-x) instead (a large positive argument, handled exactly by the case
    //above) and reciprocate -- division doesn't cancel, so it stays accurate
    //all the way down to underflow:
    if ( x < 0.0 )
      return 1.0 / NCRYSTAL_APPLY_C_NAMESPACE(detail_stable_exp)( -x );
    return 1.0 + NCRYSTAL_APPLY_C_NAMESPACE(detail_stable_expm1)( x );
  }

}

#endif

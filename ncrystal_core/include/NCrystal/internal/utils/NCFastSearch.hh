#ifndef NCrystal_FastSearch_hh
#define NCrystal_FastSearch_hh

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

#include "NCrystal/internal/utils/NCSpan.hh"

//Hints to the optimiser that the following branch condition is
//data-dependent/unpredictable, encouraging it to lower the branch as a
//conditional move (cmov on x86) instead of a mispredictable conditional
//jump. Only Clang exposes this directly (__builtin_unpredictable); GCC
//already reliably emits cmov for the plain condition without any hint, and
//other compilers fall back to the plain (correct, just not necessarily as
//fast) condition. Verified via disassembly (both GCC 15 and Clang 21,
//x86-64) that this combination gets both compilers to emit cmov for
//fastLowerBoundIdx/fastUpperBoundIdx below.
#ifdef __clang__
#  define NCRYSTAL_UNPREDICTABLE(x) __builtin_unpredictable(x)
#else
#  define NCRYSTAL_UNPREDICTABLE(x) (x)
#endif

namespace NCRYSTAL_NAMESPACE {

  //Drop-in replacements for std::lower_bound/std::upper_bound on a sorted
  //array of double, returning the same index that
  //(std::lower_bound/upper_bound(first,first+n,val)-first) would, but
  //written so that the inner comparison compiles to a branchless
  //conditional move rather than an ordinary conditional branch.
  //
  //Same precondition and behaviour as the standard library functions: first
  //must be sorted in non-decreasing order (ties/duplicate values are fine,
  //e.g. a cumulative-sum array with some zero-weight entries), n>=1.
  //
  //Motivation: while a binary search's O(log2 N) comparison count is
  //unbeatable in principle, on real hardware its actual cost is usually
  //dominated by branch *mispredictions* rather than raw comparisons or even
  //cache misses: each comparison's outcome is data-dependent and close to
  //50/50 for essentially-random queries, so the CPU's branch predictor
  //does no better than a coin flip, paying a full misprediction penalty
  //(~15-20 cycles) roughly half the time. A conditional move has no such
  //penalty. This stops mattering once the array is so large that memory
  //latency dominates everything else, but measurements (see below) show
  //the crossover is gradual, not a cliff: even well beyond that point the
  //branchless version was never slower, just decreasingly faster.
  //
  //Measured (nocheck/fastsearch/bench_search.cc, not part of the repo) on
  //log-spaced double arrays queried at random positions, on both GCC 15
  //and Clang 21 (x86-64, -O3): a robust 3-5x faster than
  //std::lower_bound/upper_bound for array sizes from a handful of
  //elements up to ~1e5-3e5 (a few hundred KB to a couple MB -- i.e.
  //comfortably within a typical machine's L2/L3 cache), fading to
  //~1.1-2.4x by 1e6-1e7 elements (tens of MB, exceeding L3) and ~1.2x by
  //3e7 (hundreds of MB, firmly into main-memory-latency territory) -- but
  //not observed to cross over into being slower at any size tried.
  //Results verified identical to the standard library's across millions
  //of random queries, arrays containing duplicate/tied values, and many
  //edge cases (see tests/src/app_fastsearch).
  std::size_t fastLowerBoundIdx( const double* first, std::size_t n, double val ) noexcept;
  std::size_t fastUpperBoundIdx( const double* first, std::size_t n, double val ) noexcept;
  inline std::size_t fastLowerBoundIdx( Span<const double> s, double val ) noexcept
  {
    return fastLowerBoundIdx( s.data(), s.size(), val );
  }
  inline std::size_t fastUpperBoundIdx( Span<const double> s, double val ) noexcept
  {
    return fastUpperBoundIdx( s.data(), s.size(), val );
  }

}

////////////////////////////
// Inline implementations //
////////////////////////////

inline std::size_t NCrystal::fastLowerBoundIdx( const double* a, std::size_t n,
                                                double val ) noexcept
{
  std::size_t lo = 0;
  std::size_t len = n;
  while ( len > 1 ) {
    std::size_t half = len/2;
    bool advance = ( a[lo+half-1] < val );
    lo += NCRYSTAL_UNPREDICTABLE(advance) ? half : 0;
    len -= half;
  }
  if ( len && NCRYSTAL_UNPREDICTABLE( a[lo] < val ) )
    ++lo;
  return lo;
}

inline std::size_t NCrystal::fastUpperBoundIdx( const double* a, std::size_t n,
                                                double val ) noexcept
{
  std::size_t lo = 0;
  std::size_t len = n;
  while ( len > 1 ) {
    std::size_t half = len/2;
    bool advance = !( val < a[lo+half-1] );//a[lo+half-1] <= val
    lo += NCRYSTAL_UNPREDICTABLE(advance) ? half : 0;
    len -= half;
  }
  if ( len && NCRYSTAL_UNPREDICTABLE( !( val < a[lo] ) ) )
    ++lo;
  return lo;
}

#undef NCRYSTAL_UNPREDICTABLE

#endif

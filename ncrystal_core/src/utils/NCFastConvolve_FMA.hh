#ifndef NCrystal_FastConvolve_FMA_hh
#define NCrystal_FastConvolve_FMA_hh

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

//Both the ordinary (NCRYSTAL_FMADISPATCH_ATTR) and the MSVC/x86-only
//Windows-fast (NCRYSTAL_WIN_FMA) variants of NCFastConvolve.cc's two hottest
//inner loops, kept side by side -- see NCFMADispatch.hh for the macros used
//below and the full mechanism/rationale, and NCFastConvolve_WINFMA.cc for the
//other half. Each function is only actually written once: the parameter list
//and body are plain, unconditional source, so there is nothing for the two
//variants to drift out of sync on.

namespace NCRYSTAL_NAMESPACE {
  //The anonymous namespace below gives the ordinary (non-extern-"C") variant
  //real internal, per-TU-private linkage -- not needed (and not applied) for
  //the extern "C" NCRYSTAL_WIN_FMA variant, which already gets external
  //linkage from "C" language linkage regardless of namespace nesting (see
  //docs/devel_fma_attribute.md), same as NCMath.cc's detail_stable_expm1/
  //detail_stable_exp:
#ifndef NCRYSTAL_WIN_FMA
  namespace {
#endif

    NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
      fastConvolveSpectralMultiply,
      ( double* ncrestrict data1, const double* ncrestrict data2, std::size_t n )
    )
    NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
      fastConvolveButterflyRun,
      ( double* ncrestrict data_j, double* ncrestrict data_sympos,
        const double* ncrestrict wtable,
        std::ptrdiff_t wtable_stride, bool is_forward, int count )
    )

    //Pointwise complex multiply of two interleaved-[re,im,re,im,...] arrays, in
    //place into the first (data1 *= data2), replacing std::complex<
    //double>::operator*= (whose exact formula/precision is otherwise up to the
    //standard library implementation -- an extra, uncontrolled source of
    //cross-platform variance, beyond just FMA contraction). Uses explicit
    //std::fma and NCRYSTAL_FMADISPATCH_ATTR. data1/data2 are ncrestrict since
    //they are always two distinct std::complex<double> buffers at the call
    //site (see NCFastConvolve.cc's convolve()); verified (not just assumed)
    //below:
    NCRYSTAL_FMADISPATCH_DECLARATOR(fastConvolveSpectralMultiply)
    ( double* ncrestrict data1, const double* ncrestrict data2, std::size_t n )
    {
      NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(fastConvolveSpectralMultiply,(data1,data2,n));
      nc_assert( buffersDisjoint( data1, 2*n, data2, 2*n ) );
      for ( std::size_t i = 0; i < n; ++i ) {
        const double a = data1[2*i], b = data1[2*i+1];
        const double c = data2[2*i], d = data2[2*i+1];
        data1[2*i]   = std::fma( a, c, -(b*d) );
        data1[2*i+1] = std::fma( a, d, b*c );
      }
    }

    //One stage of the FFT butterfly, applied to a run of count consecutive
    //j-values (see the caller in fft<is_forward> in NCFastConvolve.cc for how
    //a stage decomposes into such runs): reads/writes count complex numbers
    //each at data_j and data_sympos (both interleaved [re,im,...], with a
    //fixed stride of 2 doubles between consecutive elements of the run), and
    //reads count complex numbers from wtable (also interleaved, but strided
    //by wtable_stride doubles between consecutive elements of the run, not
    //necessarily 2). data_j/data_sympos/wtable are ncrestrict: at the call
    //site, data_j/data_sympos are always the two disjoint halves of one
    //butterfly stage's index range (data_sympos = data_j - i1*2, run length
    //i1*2, so [data_sympos,data_j) and [data_j,data_j+i1*2) are adjacent, not
    //overlapping), and wtable is always a separate, cached W-table array:
    NCRYSTAL_FMADISPATCH_DECLARATOR(fastConvolveButterflyRun)
    ( double* ncrestrict data_j, double* ncrestrict data_sympos,
      const double* ncrestrict wtable, std::ptrdiff_t wtable_stride,
      bool is_forward, int count )
    {
      NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(
        fastConvolveButterflyRun,
        (data_j,data_sympos,wtable,wtable_stride,is_forward,count) );
#ifndef NDEBUG
      nc_assert( count >= 0 && wtable_stride > 0 );
      nc_assert( buffersDisjoint( data_j, std::size_t(2*count),
                                  data_sympos, std::size_t(2*count) ) );
      //wtable's touched extent spans (count-1) strides plus one final pair:
      const std::size_t wtable_extent
        = count == 0 ? 0 : std::size_t( (count-1)*wtable_stride + 2 );
      nc_assert( buffersDisjoint( data_j, std::size_t(2*count),
                                  wtable, wtable_extent ) );
      nc_assert( buffersDisjoint( data_sympos, std::size_t(2*count),
                                  wtable, wtable_extent ) );
#endif
      for ( int k = 0; k < count; ++k ) {
        const double a = data_j[0], b = data_j[1];
        const double c = wtable[0];
        const double d = ( is_forward ? wtable[1] : -wtable[1] );
        const double jr = std::fma( a, c, -(b*d) );
        const double ji = std::fma( a, d, b*c );
        const double sr = data_sympos[0], si = data_sympos[1];
        data_j[0] = sr - jr;
        data_j[1] = si - ji;
        data_sympos[0] = sr + jr;
        data_sympos[1] = si + ji;
        data_j += 2;
        data_sympos += 2;
        wtable += wtable_stride;
      }
    }

#ifndef NCRYSTAL_WIN_FMA
  }
#endif
}

#endif

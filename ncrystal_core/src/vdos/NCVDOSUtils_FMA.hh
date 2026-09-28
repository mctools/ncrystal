#ifndef NCrystal_VDOSUtils_FMA_hh
#define NCrystal_VDOSUtils_FMA_hh

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

//The ordinary and MSVC/x86-only Windows-fast variants of
//NCVDOSUtils.cc's pwlSumAccumulateSegment -- see NCFMADispatch.hh for the
//macros used below and the full mechanism/rationale, and
//NCVDOSUtils_WINFMA.cc for the other half.

namespace NCRYSTAL_NAMESPACE {
#ifndef NCRYSTAL_WIN_FMA
  namespace {
#endif

    NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
      void, pwlSumAccumulateSegment,
      ( double* ncrestrict out, const double* ncrestrict gridPtr,
        std::size_t g, std::size_t end,
        double y0, double slope, double xLeft,
        double ylo, double yhi, double weight )
    )

    //out/gridPtr are ncrestrict: at the call site (evalPWLSum in
    //NCVDOSUtils.cc) they are always a freshly allocated output VectD and
    //the caller-provided input grid Span, which can never be the same
    //underlying allocation:
    NCRYSTAL_FMADISPATCH_DECLARATOR(void,pwlSumAccumulateSegment)
    ( double* ncrestrict out, const double* ncrestrict gridPtr,
      std::size_t g, std::size_t end,
      double y0, double slope, double xLeft,
      double ylo, double yhi, double weight )
    {
      NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(
        pwlSumAccumulateSegment,
        (out,gridPtr,g,end,y0,slope,xLeft,ylo,yhi,weight) );
      nc_assert( end >= g );
      nc_assert( buffersDisjoint( out+g, end-g, gridPtr+g, end-g ) );
      for ( std::size_t k = g; k < end; ++k ) {
        const double v = ncclamp( std::fma(slope,gridPtr[k]-xLeft,y0),
                                  ylo, yhi );
        out[k] = std::fma( weight, v, out[k] );
      }
    }

#ifndef NCRYSTAL_WIN_FMA
  }
#endif
}

#endif

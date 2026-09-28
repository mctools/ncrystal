#ifndef NCrystal_VDOSGn_FMA_hh
#define NCrystal_VDOSGn_FMA_hh

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

//The ordinary and MSVC/x86-only Windows-fast variants of NCVDOSGn.cc's
//vdosGnInterpolateDensityRun -- see NCFMADispatch.hh for the macros used
//below and the full mechanism/rationale, and NCVDOSGn_WINFMA.cc for the
//other half.

namespace NCRYSTAL_NAMESPACE {
#ifndef NCRYSTAL_WIN_FMA
  namespace {
#endif

    NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
      void, vdosGnInterpolateDensityRun,
      ( double* ncrestrict out, const double* ncrestrict f,
        const double* ncrestrict ix_as_dbl,
        const double* ncrestrict spec, std::size_t count )
    )

    //Batched form of interpolateDensity's nclerp(spec[ix],spec[ix+1],f)
    //call in NCVDOSGn.cc; must match it bit for bit (see the debug
    //assertion in interpolateDensityMany there). See
    //doc/devel_fma_attribute.md for the audit this requires. out/f/
    //ix_as_dbl/spec are ncrestrict: at the call site they are always op/
    //buf_f/buf_ix (disjoint sub-ranges of out's own buffer and workbuf's
    //two halves, verified by the caller's own no_overlap asserts plus the
    //buf_ix=buf_f+n offset) and m_spec.data() (a wholly separate,
    //long-lived member array); out/f/ix_as_dbl's mutual non-overlap (all
    //length count) is re-verified here, spec is excluded since its length
    //is not available to this function:
    NCRYSTAL_FMADISPATCH_DECLARATOR(void,vdosGnInterpolateDensityRun)
    ( double* ncrestrict out, const double* ncrestrict f,
      const double* ncrestrict ix_as_dbl,//indices, as double
      const double* ncrestrict spec, std::size_t count )
    {
      NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(
        vdosGnInterpolateDensityRun, (out,f,ix_as_dbl,spec,count) );
      nc_assert( buffersDisjoint( out, count, f, count ) );
      nc_assert( buffersDisjoint( out, count, ix_as_dbl, count ) );
      nc_assert( buffersDisjoint( f, count, ix_as_dbl, count ) );
      for ( std::size_t k = 0; k < count; ++k ) {
        const std::size_t ix = static_cast<std::size_t>(ix_as_dbl[k]);
        out[k] = nclerp( spec[ix], spec[ix+1], f[k] );
      }
    }

#ifndef NCRYSTAL_WIN_FMA
  }
#endif
}

#endif

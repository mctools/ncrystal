#ifndef NCrystal_FillHKL_FMA_hh
#define NCrystal_FillHKL_FMA_hh

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
#include "NCrystal/internal/utils/NCVector.hh"

//The per-hkl numerical kernels of NCFillHKL.cc (with their MSVC/x86-only
//Windows-fast twins in NCFillHKL_WINFMA.cc), written with every
//multiply-add shape as an explicit std::fma: the plain reciprocal-space
//transforms, phase dot-products and |F|^2 assembly were otherwise
//contracted differently by different compilers/targets (observed as
//last-bit F2 differences flipping knife-edge counts in app_fillhkl on
//the aarch64 CI legs), and this pipeline must give bit-identical
//results everywhere. See NCFMADispatch.hh for the macros and
//docs/devel_fma_attribute.md for the rules followed here.

namespace NCRYSTAL_NAMESPACE {

  namespace {
#if defined(_MSC_VER) && !defined(__clang__)
#  define NCFILLHKLFMA_ALWAYS_INLINE __forceinline
#else
#  define NCFILLHKLFMA_ALWAYS_INLINE inline __attribute__((always_inline))
#endif

    NCFILLHKLFMA_ALWAYS_INLINE
    double fillhkl_dot_fma( const Vector& a, const Vector& b )
    {
      return std::fma( a.x(), b.x(),
                       std::fma( a.y(), b.y(), a.z() * b.z() ) );
    }
  }

  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, fillhkl_ksq,
    ( const double* m9, double h, double k, double l )
  )
  NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(
    double, fillhkl_fsq,
    ( const Vector* hkl, const double* factors,
      const Vector* const* posarrs, const std::size_t* nposarr,
      std::size_t nspecies )
  )

  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_fillhkl_ksq,fillhkl_ksq)
  ( const double* m9, double h, double k, double l )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(fillhkl_ksq,(m9,h,k,l));
    //(rec_lat*hkl).mag2() with rec_lat given as its 9 row-major
    //elements:
    const double r0 = std::fma( m9[0], h, std::fma( m9[1], k, m9[2]*l ) );
    const double r1 = std::fma( m9[3], h, std::fma( m9[4], k, m9[5]*l ) );
    const double r2 = std::fma( m9[6], h, std::fma( m9[7], k, m9[8]*l ) );
    return std::fma( r0, r0, std::fma( r1, r1, r2*r2 ) );
  }

  NCRYSTAL_FMADISPATCH_DECLARATOR_C(double,detail_fillhkl_fsq,fillhkl_fsq)
  ( const Vector* hkl, const double* factors,
    const Vector* const* posarrs, const std::size_t* nposarr,
    std::size_t nspecies )
  {
    NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(fillhkl_fsq,
                                        (hkl,factors,posarrs,
                                         nposarr,nspecies));
    //|F|^2 = real^2+imag^2 with real+i*imag = sum_i factor_i *
    //sum_j exp(2pi*i*hkl.dot(pos_ij)). Inner (atomic position) sums use
    //StableSum (additions only, so contraction-proof) of the canonical
    //sincos_2pix values; the small outer (species) sum and the final
    //assembly are explicit fma:
    double re(0.0), im(0.0);
    for ( std::size_t i = 0; i < nspecies; ++i ) {
      const double factor = factors[i];
      if ( !factor )
        continue;
      StableSum cpsum, spsum;
      const Vector* pos = posarrs[i];
      const std::size_t npos = nposarr[i];
      for ( std::size_t j = 0; j < npos; ++j ) {
        const auto spcp = sincos_2pix( fillhkl_dot_fma( *hkl, pos[j] ) );
        cpsum.add( spcp.cos );
        spsum.add( spcp.sin );
      }
      re = std::fma( factor, cpsum.sum(), re );
      im = std::fma( factor, spsum.sum(), im );
    }
    return std::fma( re, re, im*im );
  }

}

#endif

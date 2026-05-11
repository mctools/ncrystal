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

#include "NCrystal/internal/sab/NCSABSurveyor.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;
namespace NCS = NCrystal::SABUtils;

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {
    namespace {

      bool rectIntersectsIdentityLine(double x0, double y0,
                                      double x1, double y1) {
        //Branchless efficient check whether or not the axis-aligned rectangle
        //with corners at (x0,y0) and (x1,y1) intersects the line y=x.
        nc_assert(x1>x0);
        nc_assert(y1>y0);
        return std::min<double>(x1, y1) >= std::max<double>(y0, x0);
      }

    }
  }
}

NCS::SABSurveyor::SABSurveyor( const SABData& sd )
  : SABSurveyor( sd.alphaGrid(), sd.betaGrid() )
{
}

NCS::SABSurveyor::SABSurveyor( const VectD& alphaGrid,
                               const VectD& betaGrid )
{
  const std::size_t na_sizet = alphaGrid.size();
  const std::size_t nb_sizet = betaGrid.size();

  //Cell index must ultimately fit in an unsigned 16bit integer (for combined
  //packing as an 32bit unsigned integer).
  if ( !( na_sizet <= 65000 && nb_sizet < 65000 ) )
    NCRYSTAL_THROW2(BadInput,"SAB grid too large (max size is 65000x65000)");

  //Verify that we have a suitable grid (nc_is_grid only in dbg builds):
  nc_assert_always( na_sizet >= 2 && nb_sizet >= 2 );
  nc_assert( nc_is_grid(alphaGrid) );
  nc_assert( nc_is_grid(betaGrid) );
  nc_assert_always( alphaGrid.front()>=0.0 );

  using idx_t = std::uint_fast32_t;
  static_assert(std::is_same<idx_t,cellidx_t>::value,"");
  const idx_t nalpha = static_cast<idx_t>(na_sizet);
  const idx_t nbeta  = static_cast<idx_t>(nb_sizet);
  const auto ncells_sizet = (na_sizet-1)*(nb_sizet-1);
  m_touch.reserve( ncells_sizet );
  m_cover.reserve( ncells_sizet );

  //The first energy value E at which a given point in (alpha,beta) space is
  //accessible is what we must find for all grid points in order to answer
  //questions about which cells are touched or covered by different
  //phase-spaces. It is given by: E/kT = (alpha-beta)^2/(4alpha). A cell will be
  //considered touched if at least one of its corners is touched, OR it
  //intersects the line alpha=beta, OR its right edge intersects the line
  //alpha=-beta. A cell will be considered covered if all of its corners are
  //touched (but the corner at highest alpha and beta does not need to be
  //checked).

  //We need a buffer of size nb_sizet, normally on the stack but using the heap
  //if needed:
  std::unique_ptr<double[]> buf_heap;
  constexpr std::size_t nsmall = 800;
  double buf_stack[nsmall];
  double* prevalpha_ebuf_begin;
  if ( nb_sizet <= nsmall ) {
    prevalpha_ebuf_begin = buf_stack;
  } else {
    buf_heap = ncmake_unique_array<double>(nb_sizet);
    prevalpha_ebuf_begin = &buf_heap[0];
  }

  //Note that cells include two alpha or beta grid points, so we treat ialpha=0
  //and ibeta=0 specially, and just iterate from ialpha=1 and ibeta=1.

  //Deal with ia=0 row, placing contents into the prevalpha_ebuf:
  {
    const double alpha0 = alphaGrid.front();
    nc_assert_always(alpha0>=0);
    double * it = prevalpha_ebuf_begin;
    if ( alpha0 == 0 ) {
      //special treatment of alpha0=0. In this case, only beta=0 is ever reached
      //before infinite energy.
      for ( auto bval : betaGrid )
        *it++ = ( bval == 0.0 ? 0.0 : kInfinity );
    } else {
      const double inv4a = 0.25 / alpha0;
      for ( auto bval : betaGrid )
        *it++ = ncsquare(alpha0-bval)*inv4a;
    }
    nc_assert(it == prevalpha_ebuf_begin+nb_sizet);
  }

  double aval_prev = alphaGrid.front();
  for ( idx_t ia = 1 ; ia < nalpha; ++ia ) {
    const double aval = vectAt(alphaGrid,ia);
    nc_assert(aval>0.0);
    const double inv4a = 0.25 / aval;
    const idx_t packidx_ib0 = ia << 16;
    double bval_prev = betaGrid.front();
    double e_prevb = ncsquare(aval-bval_prev)*inv4a;
    double * itPAEBUF = prevalpha_ebuf_begin;
    for ( idx_t ib = 1 ; ib < nbeta; ++ib ) {
      const double bval = vectAt(betaGrid,ib);
      const idx_t packidx = packidx_ib0 | ib;
      const double e = ncsquare(aval-bval)*inv4a;
      const double e_prevab = *itPAEBUF;
      const double e_preva = *std::next(itPAEBUF);
      double etouch;
     if ( bval < 0.0 && valueInInterval( aval_prev, aval, -bval ) )
        etouch = -bval;
      else
        etouch = ( rectIntersectsIdentityLine(bval_prev, aval_prev,
                                              bval, aval)
                   ? 0.0
                   : ncmin( e,e_prevb,e_prevab,e_preva) );
      const double ecover = ncmax( e_prevb,e_prevab,e_preva);
      m_touch.emplace_back( etouch, packidx );
      m_cover.emplace_back( ecover, packidx );
      *itPAEBUF++ = e_prevb;
      bval_prev = bval;
      e_prevb = e;
    }
    nc_assert(std::next(itPAEBUF) == prevalpha_ebuf_begin+nb_sizet);
    *itPAEBUF = e_prevb;
    aval_prev = aval;
  }

  nc_assert( m_touch.size() == ncells_sizet );
  nc_assert( m_cover.size() == ncells_sizet );
  std::sort( m_touch.begin(), m_touch.end() );
  std::sort( m_cover.begin(), m_cover.end() );
}

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
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;
namespace NCS = NCrystal::SABUtils;

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {
    namespace {

      inline bool rectIntersectsIdentityLine( double x0, double y0,
                                              double x1, double y1 ) {
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
  m_data.reserve( ncells_sizet );
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
    const idx_t packidx_ib0 = (ia-1) << 16;
    double bval_prev = betaGrid.front();
    double e_prevb = ncsquare(aval-bval_prev)*inv4a;
    double * itPAEBUF = prevalpha_ebuf_begin;
    for ( idx_t ib = 1 ; ib < nbeta; ++ib ) {
      const double bval = vectAt(betaGrid,ib);
      const idx_t packidx = packidx_ib0 | (ib-1);
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
      m_data.emplace_back( etouch, ecover, packidx );

      *itPAEBUF++ = e_prevb;
      bval_prev = bval;
      e_prevb = e;
    }
    nc_assert(std::next(itPAEBUF) == prevalpha_ebuf_begin+nb_sizet);
    *itPAEBUF = e_prevb;
    aval_prev = aval;
  }

  nc_assert( m_data.size() == ncells_sizet );
  std::sort( m_data.begin(), m_data.end() );

  //fixme: the following lines seem innocent, but they cost >70ms or so!!!!
  m_touch.reserve( ncells_sizet );
  for ( auto& cc : m_data )
    m_touch.emplace_back( cc.e_touch, cc.cellidx );
  nc_assert( m_touch.size() == ncells_sizet );
  std::sort( m_touch.begin(), m_touch.end() );

  m_cover.reserve( ncells_sizet );
  for ( auto& cc : m_data )
    m_cover.emplace_back( cc.e_cover, cc.cellidx );
  nc_assert( m_cover.size() == ncells_sizet );
  std::sort( m_cover.begin(), m_cover.end() );
}

NCS::SABCellSurvey::SABCellSurvey( double alpha1, double alpha2,
                                   double beta1, double beta2,
                                   double E_div_kT
                                   //fixme^^^: need typesafe E_div_kT
                                   )
{
  //fixme: special-early return for the few cases we are likely to encounter
  //mostly!
  nc_assert( !ncisnan(alpha1) );
  nc_assert( !ncisnan(alpha2) );
  nc_assert( !ncisnan(beta1) );
  nc_assert( !ncisnan(beta2) );
  nc_assert( !ncisnan(E_div_kT) );
  nc_assert( !ncisinf(alpha1) );
  nc_assert( !ncisinf(alpha2) );
  nc_assert( !ncisinf(beta1) );
  nc_assert( !ncisinf(beta2) );
  nc_assert( alpha1>=0.0 );
  nc_assert( alpha2>alpha1 );
  nc_assert( beta2>beta1 );
  nc_assert( E_div_kT > 0.0 );//fixme: we could allow also 0?
  nc_assert( !ncisinf(E_div_kT) );//fixme: we could allow inf?

  double a1 = alpha1;
  double a2 = alpha2;
  const double b1 = beta1;
  const double b2 = beta2;
  const double e = E_div_kT;

  if (b2 <= -e)
    return;//no overlap
  const double twoe = e+e;
  const double tmp = 2*std::sqrt(e*(b2+e));
  const double ap2 = twoe+b2+tmp;//alpha^+(b2)
  if ( a1 >= ap2 )
    return;//no overlap
  double am2 = twoe+b2-tmp;//alpha^-(b2)

  double am1(-1.0), ap1(-1.0);
  if ( b1 >= -e ) {
    const double tmp2 = 2*std::sqrt(e*(b1+e));
    const double twoe_plus_b1 = twoe+b1;
    am1 = twoe_plus_b1-tmp2;//alpha^-(b1)
    ap1 = twoe_plus_b1+tmp2;//alpha^+(b1)
  }

  //snap along alpha:
  a2 = ncmin(a2,ap2);
  if ( b2 < 0 ) {
    a1 = ncmax(a1,am2);
  } else if ( b1 > 0 ) {
    a1 = ncmax(a1,am1);
  }
  if (!(a2 > a1))
    return;//no overlap

  struct Interval final {
    Interval( double aa1, double aa2, bool bb )
      : a1(aa1), a2(aa2), bounded(bb) { nc_assert( a2 > a1 ); }
    double a1, a2;
    bool bounded;
  };
  //fixme: Can we do it without smallvector?
  SmallVector<Interval,3> intervals_lower, intervals_upper;

  //Look at lower bounds:
  if ( b1 <= -e ) {
    intervals_lower.emplace_back(a1,a2,true);//bounded by beta^-(alpha)
  } else {
    double aa1 = a1;
    nc_assert(am1!=-1.0);
    const double ncmin_am1_a2 = ncmin(am1,a2);
    if ( aa1 < ncmin_am1_a2 ) {
      intervals_lower.emplace_back(aa1,
                                   ncmin_am1_a2,true);//bounded by beta^-(alpha)
      aa1 = ncmin_am1_a2;
    }
    const double ncmin_ap1_a2 = ncmin(ap1,a2);
    if ( aa1 < ncmin_ap1_a2 ) {
      intervals_lower.emplace_back(aa1,ncmin_ap1_a2,false);//bounded by b1 edge
      aa1 = ncmin_ap1_a2;
    }
    if ( aa1 < a2 )
      intervals_lower.emplace_back(aa1,a2,true);//bounded by beta^-(alpha)
  }
  //Look at upper bounds:
  if ( b2 <= 0.0 )
    am2 = 0.0;
  {
    double aa1 = a1;
    const double ncmin_am2_a2 = ncmin(am2,a2);
    if ( aa1 < ncmin_am2_a2 ) {
      intervals_upper.emplace_back(aa1,ncmin_am2_a2,true);//bounded by beta^+(alpha)
      aa1 = ncmin_am2_a2;
    }
    if ( aa1 < a2 )
      intervals_upper.emplace_back(aa1,a2,false);//bounded by b2 edge
  }

  //Combine intervals:
  // (FIXME: Special case intervals_lower.size()==1 and intervals_upper.size() == 1?
  nc_assert( !intervals_lower.empty() && !intervals_upper.empty() );
  nc_assert( intervals_lower.front().a1 == intervals_upper.front().a1 );
  nc_assert( intervals_lower.back().a2 == intervals_upper.back().a2 );
  double au = intervals_upper.back().a2, al;

  while ( !intervals_lower.empty() ) {
    bool bm_bounded =  intervals_lower.back().bounded;
    bool bp_bounded =  intervals_upper.back().bounded;
    nc_assert_always(!intervals_upper.empty());
    if ( intervals_lower.back().a1 >= intervals_upper.back().a1 ) {
      al = intervals_lower.back().a1;
      intervals_lower.pop_back();
    } else {
      al = intervals_upper.back().a1;
      intervals_upper.pop_back();
    }
    if ( au > al ) {
      m_regions.emplace_back();
      nc_assert(au > al);
      m_regions.back().alpha_low = al;
      m_regions.back().alpha_up = au;
      m_regions.back().is_bounded_by_betaminus = bm_bounded;
      m_regions.back().is_bounded_by_betaplus = bp_bounded;
    }
    au = al;
  }
  nc_assert( intervals_upper.size()==1 );
  const auto nsize = m_regions.size();
  nc_assert( nsize <= nmax_regions );

  if ( nsize == 0 )
    return;

  //Testing with various integration algs indicates that we benefit from more
  //careful splitting around the phasespace endpoints. This should happen very
  //rarely. Note that we have incremented nmax_regions by 2 due to this!
  std::size_t isplit = nsize;
  for ( auto i : ncrange( nsize ) ) {
    auto& r = m_regions[i];
    if ( r.is_bounded_by_betaminus &&
         e > r.alpha_low && e < r.alpha_up )
      isplit = i;
  }
  if ( isplit < nsize ) {
    //split at endpoint at beta=-e
    RegionList newregions;
    for ( auto i : ncrange( nsize ) ) {
      auto& r = m_regions[i];
      newregions.push_back(r);
      if ( i == isplit ) {
        newregions.back().alpha_low = e;
        newregions.push_back(r);
        newregions.back().alpha_up = e;
      }
    }
    std::swap(newregions,m_regions);
  }

  //split at endpoint near alpha=0:
  const double aal = m_regions.back().alpha_low;
  const double aau = m_regions.back().alpha_up;
  if ( aal < 0.01*aau ) {
    double asplit = 0.05*aau;
    if ( aal < 0.001*aau )
      asplit = 0.005*aau;
    if ( aal < 0.0001*aau )
      asplit = 0.0005*aau;
    m_regions.push_back(m_regions.back());
    m_regions.back().alpha_up = asplit;
    m_regions.at(m_regions.size()-2).alpha_low = asplit;
  }

}

void NCS::SABCellSurvey::toJSON( std::ostream& os ) const
{
  os << ("{\"regions_format\":[\"alphalow\",\"alphaup\","
         "\"is_bounded_by_beta-\",\"is_bounded_by_beta+\"],\"regions\":[");
  for( auto i : ncrange( m_regions.size() ) ) {
    auto& r = m_regions.at(i);
    if ( i )
      os << ',';
    os << '[';
    streamJSON(os,r.alpha_low);
    os<<',';
    streamJSON(os,r.alpha_up);
    os<<',';
    streamJSON(os,r.is_bounded_by_betaminus);
    os<<',';
    streamJSON(os,r.is_bounded_by_betaplus);
    os<<']';
  }
  os<<"]}";
}

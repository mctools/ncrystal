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

#define NCRYSTAL_USE_CELLINFO_RADIX_SORT

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {
    namespace {

      using CellInfo = SABSurveyor::CellInfo;

#ifdef NCRYSTAL_USE_CELLINFO_RADIX_SORT

      template <typename FuncDblKey, typename FuncTieBreak>
      void radixSortByDouble( std::unique_ptr<CellInfo[]>& dataArray,
                              std::size_t n,
                              FuncDblKey key_extract,
                              FuncTieBreak tie_breaker) {
        //Specialised radix sort, intended to sort structs where a single double
        //field is the primary sorting key (with tie-breaker for other
        //fields. For simplicity we just hardwirde the CellInfo type, but it
        //could be templated.
        //
        //The double fields are assumed to be non-negative finite values, or
        //+infinity. Negative zero is NOT allowed.

        using T = CellInfo;

        if (n < 2)
          return;

        std::unique_ptr<T[]> buffer = ncmake_unique_array_noinit<T>(n);
        T* src = dataArray.get();
#ifndef NDEBUG
        //Check validity (needed for our std::memcpy trick to work, see below):
        for ( std::size_t i = 0; i < n; ++i ) {
          double val = key_extract(src[i]);
          nc_assert_always( val >= 0.0 && !std::signbit(val) );
        }
#endif
        T* dst = buffer.get();

        constexpr int BITS_PER_PASS = 16;//NB: BUCKET=2**16 is large, but we
                                         //allocate on the heap to be safe.
                                         //We could reduce to 8 if needed, but
                                         //that was slightly less efficient in
                                         //profiling.
        constexpr int TOTAL_PASSES  = ( sizeof(double) * 8 ) / BITS_PER_PASS;
        constexpr int BUCKETS       = 1 << BITS_PER_PASS;
        constexpr int MASK          = BUCKETS - 1;
        auto tmpbuf = ncmake_unique_array_noinit<std::size_t>(2*BUCKETS);

        std::size_t * counts = tmpbuf.get();
        std::size_t * offsets = counts + BUCKETS;
        auto resetCounts = [&counts]() { std::fill_n(counts, BUCKETS,
                                                     std::size_t{0}); };
        static_assert( TOTAL_PASSES == 4, "" );
        static_assert( BUCKETS == 65536, "" );

        for (int pass = 0; pass < TOTAL_PASSES; ++pass) {
          const int shift = pass * BITS_PER_PASS;
          auto radix = [shift,key_extract,src](std::size_t i)
          {
            double v = key_extract(src[i]);
            //The trick is to interpret the double bits as an uint64_t, and then
            //extract BITS_PER_PASS in each pass. With our constraints of
            //non-negative (no -0.0 but +inf allowed) integers, this should have
            //the same sorting with standard double implementations:
            static_assert(std::numeric_limits<double>::is_iec559,"");
            static_assert(sizeof(double) == sizeof(std::uint64_t), "");
            std::uint64_t bits;
            std::memcpy(&bits, &v, sizeof(double));
            int res = (bits >> shift) & MASK;
            nc_assert( res >= 0 && res < BUCKETS );
            return res;
          };
          resetCounts();
          for (std::size_t i = 0; i < n; ++i)
            ++counts[radix(i)];
          offsets[0] = 0;
          for (int i = 1; i < BUCKETS; ++i)
            offsets[i] = offsets[i - 1] + counts[i - 1];
          for (std::size_t i = 0; i < n; ++i)
            dst[offsets[radix(i)]++] = src[i];
          //every second pass we use the original array as the new target
          //buffer for the next step:
          std::swap(src, dst);
        }

        //Even number of swaps should have left the results in dataArray.data():
        static_assert( TOTAL_PASSES%2 == 0, "" );
        nc_assert( src == dataArray.get() );

        //Release buffer:
        dst = nullptr;
        buffer = nullptr;

        //Radix sort done on primary double key, now resolve ties:
        double key0 = key_extract(src[0]);
        const double nm1 = n-1;
        for (std::size_t i = 1, i0 = 0; i <= nm1; ++i) {
          double key = key_extract(src[i]);
          if ( i == nm1 || key != key0 ) {
            if ( i - i0 > 1 ) {
              std::sort( src + i0, src + i, tie_breaker);
            }
            i0 = i;
            key0 = key;
          }
        }
      }
#endif// NCRYSTAL_USE_CELLINFO_RADIX_SORT

      inline void sortCellInfo( std::unique_ptr<CellInfo[]>& data,
                                std::size_t n )
      {
#ifdef NCRYSTAL_USE_CELLINFO_RADIX_SORT
        //Note: The comparisons functions in the next statement must be carefully
        //aligned with the CellInfo::operator<(..) implementation!
        radixSortByDouble( data, n,
                           []( CellInfo& c ) -> double& { return c.e_touch; },
                           []( const CellInfo& a, const CellInfo& b) -> bool {
                             return  ( ( a.e_cover != b.e_cover )
                                       ? ( a.e_cover < b.e_cover )
                                       : ( a.cellidx.val < b.cellidx.val ) );
                           } );
#else
        std::sort( data.get(), data.get()+n );
#endif
        //As a sanity check, std::is_sorted should yield true if everything was
        //implemented consistently:
        nc_assert_always( std::is_sorted(data.get(), data.get()+n) );//fixme _always
      }

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
  SABIdx::verifyGridSizes(na_sizet,nb_sizet);

  //Verify that we have a suitable grid (nc_is_grid only in dbg builds):
  nc_assert( na_sizet >= 2 && nb_sizet >= 2 );
  nc_assert( nc_is_grid(alphaGrid) );
  nc_assert( nc_is_grid(betaGrid) );
  nc_assert_always( alphaGrid.front()>=0.0 );

  using idx_t = cellidx_t::index_t;
  static_assert(std::is_same<idx_t,SABIdx::PackedIndex::index_t>::value,"");
  const idx_t nalpha = static_cast<idx_t>(na_sizet);
  const idx_t nbeta  = static_cast<idx_t>(nb_sizet);
  const auto ncells_sizet = (na_sizet-1)*(nb_sizet-1);
  auto cellinfo_array = ncmake_unique_array_noinit<CellInfo>(ncells_sizet);
  auto itData = cellinfo_array.get();

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
    buf_heap = ncmake_unique_array_noinit<double>(nb_sizet);
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
    static_assert( std::is_same<idx_t,std::uint32_t>::value, "" );
    const idx_t packidx_ib0 = (ia-1) << 16;
    nc_assert( packidx_ib0
               == SABIdx::PackedIndex::fromAlphaIdxBetaIdx(ia-1,0).val );
    double bval_prev = betaGrid.front();
    double e_prevb = ncsquare(aval-bval_prev)*inv4a;
    double * itPAEBUF = prevalpha_ebuf_begin;
    for ( idx_t ib = 1 ; ib < nbeta; ++ib ) {
      const double bval = vectAt(betaGrid,ib);
      const idx_t packidx = packidx_ib0 | (ib-1);
      nc_assert( packidx
                 == SABIdx::PackedIndex::fromAlphaIdxBetaIdx( ia-1,ib-1).val );
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
      itData->e_touch = etouch;
      itData->e_cover = ecover;
      itData->cellidx = SABIdx::PackedIndex{packidx};
      ++itData;
      *itPAEBUF++ = e_prevb;
      bval_prev = bval;
      e_prevb = e;
    }
    nc_assert(std::next(itPAEBUF) == prevalpha_ebuf_begin+nb_sizet);
    *itPAEBUF = e_prevb;
    aval_prev = aval;
  }
  nc_assert( itData == cellinfo_array.get() + ncells_sizet );

  //Finally, we must sort the cells. Profiling shows this to be a very important
  //bottleneck for material initialisation time, and this is improved by the
  //usage of a custom radix sort instead of a simple std::sort invocation.
  sortCellInfo( cellinfo_array, ncells_sizet );

  //All done:
  m_dataHolder = std::move(cellinfo_array);
  m_dataSpan = Span<const CellInfo>( m_dataHolder.get(),
                                     m_dataHolder.get() + ncells_sizet );
}

NCS::SABCellSurvey::SABCellSurvey( double alpha1, double alpha2,
                                   double beta1, double beta2,
                                   double E_div_kT
                                   //fixme^^^: need typesafe E_div_kT
                                   )
{
  //fixme: special-early return for the few cases we are likely to encounter
  //mostly?
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
    Interval( double aa1, double aa2, bool bb ) ncnoexceptndebug
      : a1(aa1), a2(aa2), bounded(bb) { nc_assert( a2 > a1 ); }
    double a1, a2;
    bool bounded;
  };

  TinyVector<Interval,3> intervals_lower, intervals_upper;

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
  double au = intervals_upper.back().a2;
  double al;

  while ( !intervals_lower.empty() ) {
    bool bm_bounded =  intervals_lower.back().bounded;
    bool bp_bounded =  intervals_upper.back().bounded;
    nc_assert(!intervals_upper.empty());
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

static_assert(std::is_trivially_default_constructible
              <NC::SABIdx::PackedIndex>::value,"");

#ifndef NCrystal_SABCellSample_hh
#define NCrystal_SABCellSample_hh

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

#include "NCrystal/internal/sab/NCSABCellInteg.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {

    /////////////////////////////////////
    // Utilities for SAB Cell sampling //
    /////////////////////////////////////

    class LogLinDistSampler final {
    public:

      //////////////////////////////////////////////////////////////////////////
      // Class which can be used to sample alpha values according to the usual
      // log-linear alpha interpolation scheme of S(alpha,beta) tables.
      //
      // Technically this means sampling a value of x on the interval [a,b]
      // according to a function f(x) which is log-linear (linear in log(f(x)))
      // when both f(a)>0 and f(b)>0, and otherwise linear in f(x). The
      // constructor also supports force_linlin=true, which can be used to for
      // the interpolation to always be linear in f(x). This is needed when
      // sampling correctly over parts of the alpha-range of an S(alpha,beta)
      // cell where fa or fb is zero at the original corners.
      //
      // For efficiency, pre-cached values of log(f(a)) and log(f(b)) must be
      // provided, but they will be completely ignored in case of linear
      // interpolation.

      LogLinDistSampler( double a, double fa, double logfa,
                         double b, double fb, double logfb,
                         bool force_linlin = false );

      double sample( RNG& ) const;

    private:
      bool m_islinlin;
      double m_cache1, m_cache2, m_cache3, m_cache4;
      double m_a, m_b;
    };

    class FullCellSampler final : private NoCopyMove  {
    public:

      //////////////////////////////////////////////////////////////////////////
      //Sample (alpha,beta) point over a full cell. Note that for efficiency,
      //only a pointer to the CellData object is copied by the constructor, so
      //the .sampleAlphaBeta(..)  method should not be called after the lifetime
      //of the CellData object has ended. If only a single sampling is needed,
      //the static .sampleOneAlphaBeta(..)  method can be used, avoiding any
      //lifetime concerns.

      FullCellSampler( const CellData* );

      PairDD sampleAlphaBeta( RNG& rng );

      static PairDD sampleOneAlphaBeta( const CellData& c, RNG& rng )
      {
        return FullCellSampler(&c).sampleAlphaBeta(rng);
      }

    private:
      Optional<LogLinDistSampler> m_samplers[2];
      const CellData* m_cellptr;
      double m_prob1;
    };

    class BoundedCellSampler final : private NoCopyMove {

      //////////////////////////////////////////////////////////////////////////
      //Sample (alpha,beta) in the intersection between a given cell and the
      //kinematic bounds given by the phasespace available to a particular
      //neutron E/kT.

    public:

      // To construct a BoundedCellSampler object, expensive preprocessing is
      // required. To facilitate better caching and initialisation performance,
      // such info is kept in a ProcessedSurveyInfo object which can be created
      // ahead of time and is optimised for storage in a contiguous data
      // structure elsewhere. The initialisation is specific to a given neutron
      // energy, and requires that a SABCellSurvey is already available. All
      // parameters provided must be identical to those with which the
      // SABCellSurvey was constructed.

      struct ProcessedSurveyInfo {
        //Overlay values and (if set), alpha-range in which samples must fall.
        //Note: overlay values are calculated as S(bmiddle)*bwidth*(b2-b1) + a
        //      bit of safety.
        float overlay1, overlay2;
        float alpha_low;//-1 if not set
        float alpha_up;//-1 if not set
      };
      static ProcessedSurveyInfo
      processSurveyInfo( const SABCellSurvey&,
                         double alpha1, double alpha2,
                         double beta1, double beta2,
                         double E_div_kT );

      //Initialise. The caller must ensure that all arguments are consistent
      //(e.g. ProcessedSurveyInfo must be created from the same CellData and
      //E_div_kT values as those passed along). The probability_b1_edge must be
      //given by W1/(W1+W2) where W1 is the alpha integral along the beta1 edge
      //of the cell and W2 the one along the beta2 edge:
      BoundedCellSampler( const CellData& cell,
                          const ProcessedSurveyInfo& oi,
                          double E_div_kT,
                          double probability_b1_edge );

      //Actual sampling method. For diagnostics purposes, this includes a count
      //of the internal number of tries (only the last of these was accepted).
      struct Result {
        double alpha;
        double beta;
        std::uint_fast64_t ntries;
      };
      Result sampleAlphaBeta( RNG& );

    private:
      struct EdgeData {
        bool islinlin;
        double s_low, s_up, lns_low, lns_up;
      };
      CellData m_cell;
      Optional<EdgeData> m_edge[2];
      double m_overlay[2];
      double m_4e;
      double m_prob1;
      double m_restrict_a1;
      double m_restrict_a2;
      void initEdge( EdgeData&, unsigned ) const;
    };
  }
}

////////////////////////////
// Inline implementations //
////////////////////////////

inline double NCrystal::SABUtils::LogLinDistSampler::sample( RNG& rng ) const
{
  double val;
  if ( m_islinlin ) {
    if ( m_cache2 > m_cache1 && rng.generate()*m_cache2 > m_cache1 ) {
      val = m_cache3 + m_cache4 * std::sqrt(rng.generate());
    } else {
      val = m_cache3 + m_cache4 * rng.generate();
    }
  } else {
    nclikely val = m_a + m_cache1 * randKPowX(m_cache2,m_cache3,rng.generate());
  }
  return ncclamp( val, m_a, m_b );
}

inline NCrystal::SABUtils::BoundedCellSampler::
BoundedCellSampler( const CellData& cell,
               const ProcessedSurveyInfo& psi,
               double E_div_kT,
               double probability_b1_edge )
  : m_cell(cell),
    m_4e(4.0*E_div_kT),
    m_prob1(probability_b1_edge),
    m_restrict_a1(static_cast<double>(psi.alpha_low)),
    m_restrict_a2(static_cast<double>(psi.alpha_up))
{
  nc_assert( m_prob1 >= 0 && m_prob1 <= 1 );
  nc_assert( E_div_kT > 0 && std::isfinite(E_div_kT) );
  m_overlay[0] = static_cast<double>(psi.overlay1);
  m_overlay[1] = static_cast<double>(psi.overlay2);
  if ( !( m_restrict_a1>0.0 && m_restrict_a1 < cell.a2 ) )
    m_restrict_a1 = cell.a1;
  if ( !(m_restrict_a2 > m_restrict_a1 && m_restrict_a2 < cell.a2) )
    m_restrict_a2 = cell.a2;

  nc_assert( m_overlay[0] > 0.0 );
  nc_assert( m_overlay[1] > 0.0 );
  nc_assert( std::isfinite(m_overlay[0]) );
  nc_assert( std::isfinite(m_overlay[1]) );
}

#endif

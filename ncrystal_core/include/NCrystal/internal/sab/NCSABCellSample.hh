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

    class LogLinDistSampler {
    public:

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

    class FullCellSampler : private NoCopyMove  {
    public:

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
      val = m_cache3 + m_cache4 * std::sqrt(rng.generate());//fixme: sqrt(R) xvs minmax(R,R)
    } else {
      val = m_cache3 + m_cache4 * rng.generate();
    }
  } else {
    nclikely val = m_a + m_cache1 * randKPowX(m_cache2,m_cache3,rng.generate());
  }
  return ncclamp( val, m_a, m_b );
}

#endif

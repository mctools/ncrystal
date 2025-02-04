#ifndef NCrystal_SABSamplerModels_hh
#define NCrystal_SABSamplerModels_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
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

#include "NCrystal/internal/sab/NCSABSampler.hh"
#include "NCrystal/internal/utils/NCPointwiseDist.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace SAB {

    class SABSamplerAtE_Alg1 : public SABSamplerAtE {
      //A sampler which implements Algorithm1 of the Algorithm 1. of the
      //sampling paper (https://doi.org/10.1016/j.jcp.2018.11.043).
      //
      //On top of that, for NCrystal v3.1.0, a special treatment of the first
      //beta-bin was added, since the piece-wise linear assumption was too crude
      //there.

    public:
      PairDD sampleAlphaBeta(double ekin_div_kT, RNG&) const final;

      struct CommonCache {
        const std::shared_ptr<const SABData> data;
        const VectD logsab, alphaintegrals_cumul;
      };
      class AlphaSampleInfo  {
        //Class able to sample alpha for a given energy and beta-value.
      public:
        struct SAPoint {
          double alpha = 0, sval = 0, logsval = 0;
          unsigned alpha_idx = 0;//the grid idx by which the front/back tail is bounded.
        };
        SAPoint pt_front, pt_back;
        double prob_front = 0;//1.0 means narrow, 2.0 means 0 cross-section at value, sample linearly in [pt_front.alpha,pt_back.alpha]
        double prob_notback = 0;//prob_front+prob_middle
      };

      SABSamplerAtE_Alg1( std::shared_ptr<const CommonCache>,
                          VectD&& betaVals,
                          VectD&& betaWeights,
                          std::vector<AlphaSampleInfo>&&,
                          std::size_t ibetaOffset,
                          double firstBinKinematicEndpointValue = 1.0 );

    private:
      // Sample alpha from F(alpha|beta_j,Ei) (line 7-8 of Alg. 1 in the sampling
      // paper). NB: this needs to work with a single random number, the
      // percentile, for purposes of interpolating between two beta-rows:
      double sampleAlpha(std::size_t ibeta, double rand_percentile) const;

      //Data:
      std::shared_ptr<const CommonCache> m_common;
      PointwiseDist m_betaSampler;
      std::vector<AlphaSampleInfo> m_alphaSamplerInfos;
      std::size_t m_ibetaOffset;
      double m_firstBinKinematicEndpointValue;//if <=0, lowest beta value in
                                              //m_betaSampler has been moved
                                              //down -1/3*firstbinwidth to fake
                                              //the weight of a sqrt-like curve
                                              //rather than a trapezoidal curve
                                              //in the first bin.
    };

    class SABSamplerAtE_NoScatter : public SABSamplerAtE {
      //Special technical sampler which doesn't actually scatter (i.e. returns
      //alpha=beta=0). For usage of edge-cases with vanishing cross-section.
    public:
      PairDD sampleAlphaBeta(double, RNG&) const final { return {0.0,0.0}; }
    };

#if 0
    class SABSamplerAtE_Gridded : public SABSamplerAtE {
      //A sampler which improves the rejection-sampling described in
      //https://doi.org/10.1016/j.jcp.2018.11.043 by expanding the kinematically
      //accessible domain at Ei, D, slightly so it aligns with grid
      //cell-boundaries. It then samples (alpha,beta) on this expanded D by
      //first picking a grid-cell according to the integral of S(alpha,beta)
      //over that cell, and secondly sampling (alpha,beta) within that
      //cell. Finally, if the chosen (alpha,beta) algorithm is not within the
      //kinematically accessible region of a given neutron ekin, the point is
      //rejected and sampling starts over. This sampler should provide more
      //mathematically sound result, without dubious interpolation artifacts,
      //but becomes inefficient for neutrons at very small energies, where the
      //kinematically accessible region becomes much narrower than the grid
      //cells.
      //...
    };
#endif
  }
}

#endif

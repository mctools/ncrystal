#ifndef NCrystal_SABKBCellSmpl_hh
#define NCrystal_SABKBCellSmpl_hh

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

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/phys_utils/NCKinUtils.hh"
#include "NCrystal/interfaces/NCSABData.hh"
#include "NCrystal/internal/sab/NCScatKnlData.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {

    //Utilities for sampling within a cell, respecting kinematic bounds.

    class ParabolicBandBoxSampler final {
    public:

      //Helper class which has a quick initialisation and which can efficiently
      //sample (x,y) uniformly in the 2D plane within the intersection between a
      //box (0<=x0<x<x1, 0<=y0<y<y1) and the shape given by x>=0 and
      //(1-sqrt(x))^2<y<(1+sqrt(x))^2:

    public:
      ParabolicBandBoxSampler( double x0, double y0,
                               double x1, double y1 );
      bool canSample() const { return m_data[1] != -1.0; }
      PairDD sample( RNG& rng ) const
      {
        if ( fullyInside() )
          return sampleFullBox(rng);
        else
          return sampleCrossing(rng);
      }

      //For debugging:
      void toJSON(std::ostream&) const;

    private:
      double m_data[4];
      bool fullyInside() const { return m_data[2] > 0.0; }
      PairDD sampleCrossingParallelogram( RNG& rng ) const;
      PairDD sampleCrossing( RNG& rng ) const
      {
        if ( m_data[3]<0.0 )
          return sampleCrossingParallelogram(rng);
        while (true) {
          auto res = sampleFullBox(rng);
          if ( ncsquare(res.second-(1.0+res.first)) < 4.0 * res.first )
            return res;
        }
      }
      PairDD sampleFullBox( RNG& rng ) const
      {
        nc_assert(canSample());
        const double x0(m_data[0]), y0(m_data[1]), y1(m_data[3]);
        const double x1(ncabs(m_data[2]));
        PairDD res;
        res.first = ncmin(x1,x0 + rng.generate()*(x1-x0));
        res.second = ncmin(y1,y0 + rng.generate()*(y1-y0));
        return res;
      }
    };


     class SABKBCellSampler final {

       // For a given cell defined by an alpha and beta range, this class can be
       // used to efficiently sample an (alpha,beta) range uniformly in the
       // intersection between the cell and the phase-space accessible to a
       // neutron of a specific E/kT value.
       //
       // This class is intended to have a relatively fast initialisation time,
       // and a rather efficient sampling. Meaning that it can be initialised
       // and used on-demand during samplings. Although it is best used where
       // this class only initialised once and then used for multiple samplings
       // to handle the S-value variations across the cell. It could in
       // principle be cached if needed be though, if the 32bytes of storage is
       // deemed OK.

     public:

       SABKBCellSampler( double ekin_div_kT,
                         double beta0, double alpha0,
                         double beta1, double alpha1 )
         : SABKBCellSampler( internal_kT_div_ekin_t{},
                             [ekin_div_kT]()
                             {
                               nc_assert(ekin_div_kT>0.0);
                               nc_assert(std::isfinite(ekin_div_kT));
                               const double kT_div_ekin = 1.0/ekin_div_kT;
                               nc_assert(kT_div_ekin>0.0);
                               nc_assert(std::isfinite(kT_div_ekin));
                               return kT_div_ekin;
                             }(),
                             beta0, alpha0, beta1, alpha1 )
       {
#ifndef NDEBUG
         m_dbg_e = ekin_div_kT;
#endif
       }


       //Note that in order to preserve memory, the sampling requires
       //ekin_div_kT as an argument. This MUST be the same as was used in the
       //constructor.
       PairDD sampleAlphaBeta( RNG& rng, double ekin_div_kT )
       {
         nc_assert(m_dbg_e == ekin_div_kT);
         PairDD xy = m_pbbsampler.sample(rng);
         return { ekin_div_kT*xy.second, ekin_div_kT*(xy.first-1.0) };
       }
     private:
       ParabolicBandBoxSampler m_pbbsampler;
#ifndef NDEBUG
       double m_dbg_e = 0.0;
#endif
       struct internal_kT_div_ekin_t {};
       SABKBCellSampler( internal_kT_div_ekin_t,
                         double kT_div_ekin,
                         double beta0, double alpha0,
                         double beta1, double alpha1 )
         : m_pbbsampler( 1.0 + beta0*kT_div_ekin, alpha0*kT_div_ekin,
                         1.0 + beta1*kT_div_ekin, alpha1*kT_div_ekin )
       {
       }


     };

  }
}

#endif

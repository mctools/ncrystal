#ifndef NCrystal_MMC_Utils_hh
#define NCrystal_MMC_Utils_hh

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/minimc/NCMMC_Basket.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace Utils {

      //Internal utilities that are useful when implementing MMC simulation
      //engines.

      void calcProbTransm( NumberDensity nd, std::size_t N,
                           const double * ncrestrict xs_or_nullptr,
                           const double * ncrestrict dist,
                           double * ncrestrict out );

      void propagate( NeutronBasket& b,
                      const double* ncrestrict dists );

      void propagateAndAttenuate( NeutronBasket& b,
                                  NumberDensity nd,
                                  const double* ncrestrict dists,
                                  const double* ncrestrict xsvals = nullptr );

      void sampleRandDists( RNG& rng, NumberDensity nd,
                            const double * ncrestrict dists,
                            const double * ncrestrict xsvals,
                            std::size_t N,
                            double * ncrestrict tgt );

      void scatterGivenMu( RNG& rng,
                           NeutronBasket& b,
                           double * ncrestrict mu_vals );

    }
  }
}

#endif

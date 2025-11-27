#ifndef NCrystal_MMC_Utils_hh
#define NCrystal_MMC_Utils_hh

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/minimc/NCMMC_Basket.hh"

namespace NCRYSTAL_NAMESPACE {

  class MatCfg;

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


      //Fixme: Move the following explanation somewhere more appropriate and
      //maintainable, and refer to it there?
      //
      // The decodeScenario function parses a MiniMC quick simulation scenario
      // string, according to the syntax:
      //
      //  "ENERGY [pencil] [on [THICKNESS] [sphere|slab]] [COUNT times]"
      //
      //  The meaning of the various parameters (in uppercase above) and
      //  keywords (in lowercase above) is explained in the following.
      //
      //  If given, COUNT is the initial number of neutrons to input into the
      //  simulation. Otherwise, the default value is 1e6 for isotropic
      //  materials, and 1e5 for anisotropic materials.
      //
      //  ENERGY is the monochromatic beam energy like "1.8Aa", "25meV" or
      //  "0.1eV". Special units "BT" means the bragg threshold of the material
      //  (or 4.0Aa in case material does not have one), and "kT" means a
      //  kinetic energy equal to Boltzmann's constant times the material
      //  temperature. Using either special unit will cause the result to be
      //  rounded to 6 significant digits.
      //
      //  "pencil" is an optional keyword related to the beam profile (see below).
      //
      //  THICKNESS is the material thickness like "1mm", "2m", "0.4cm", or
      //  "2.5mfp". The unit "mfp" corresponds to the mean-free-path length for a
      //  neutron scattering interaction in the material (rounded to 6 significant
      //  digits). Default THICKNESS is "1mfp".
      //
      // The keywords "sphere" or "slab" can be used to select the sample geometry
      // (default is a sphere).
      //
      // For spherical geometry, the beam profile will by default be taken to be a
      // beam with a uniform circular profile, of the same radius as the
      // sphere. However, if the keyword "pencil" is provided, a pencil beam hitting
      // the sphere centrally is used instead,
      //
      // As a special case, an empty scenario string is interpreted in the same way
      // as a scenario string with contents "0.8BT".
      //
      // For flexibility and usage from the cmdline, colons (:) and underscores (_)
      // can be used as whitespace. Additionally, all repeated whitespace (tabs,
      // newlines, etc.) is converted into a single space before parsing, and
      // trailing or leading whitespace is trimmed away.
      //

      struct ScenarioDecoded {
        std::string geomcfg;
        std::string srccfg;
        std::string enginecfg;
      };
      ScenarioDecoded decodeScenario( const MatCfg&, const char* scenario );

    }
  }
}

#endif

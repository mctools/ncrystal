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

#include "NCrystal/internal/minimc/NCMMC_Defs.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  class MatCfg;

  namespace MiniMC {

    namespace Utils {

      //Internal utilities that are useful when implementing MMC simulation
      //engines.

      //Transmission probability. Handles xs=inf always, and handles dist=inf
      //correctly if geom_is_unbounded=true.
      //
      // The result will be:
      //  dist = inf : P = 0.0 (NOTE: Requires geom_is_unbounded=true!)
      //  dist = 0   : P = 1.0
      //  0<dist<inf : P = exp(-k*xs*dist) (k=numdens*100.0)
      //
      //Note that for (xs,dist)=(inf,0) we return P=0.0 and not P=1.0, and this
      //choice is to make tallying of unbounded geometries easier.
      void calcProbTransm( NumberDensity nd, std::size_t N,
                           bool geom_is_unbounded,
                           const double * ncrestrict xs_or_nullptr,
                           const double * ncrestrict dist,
                           double * ncrestrict out );

      //Moves the neutron forward. If geom_is_unbounded, this also works as
      //expected for dist=inf (e.g. moving a neutron with any pos and
      //dir=(0,-1/sqrt2,1/sqrt2) forward dist=inf, will result in a neutron with
      //position (0,-inf,inf). Notably this will NOT give (nan,-inf,inf).
      void propagate( NeutronBasket& b,
                      bool geom_is_unbounded,
                      const double* ncrestrict dists );

      //Combines the propagate and calcProbTransm functions in order to move the
      //neutrons forward. Notably if dist=0 this will leave the neutron
      //unchanged and if dist=inf (requires geom_is_unbounded=true) it will be
      //left with weight=0.
      void propagateAndAttenuate( NeutronBasket& b,
                                  NumberDensity nd,
                                  bool geom_is_unbounded,
                                  const double* ncrestrict dists,
                                  const double* ncrestrict xsvals = nullptr );

      void sampleRandDists( RNG& rng, NumberDensity nd,
                            const double * ncrestrict dists,
                            const double * ncrestrict xsvals,
                            std::size_t N,
                            double * ncrestrict tgt );

      //Helper function suitable for implementing distToVolumeExit (see
      //NCMMC_Geom.hh) for an infinite slab bounded by planes at
      //x = +-slab_halfthickness.
      void distToSlabExit( const double * ncrestrict x,
                           const double * ncrestrict ux,
                           double * ncrestrict out_dist,
                           std::size_t n, double slab_halfthickness );

      //Helper function suitable for implementing distToVolumeEntry (see
      //NCMMC_Geom.hh) for an infinite slab bounded by planes at
      //x = +-slab_halfthickness.
      void distToSlabEntry( const double * ncrestrict x,
                            const double * ncrestrict ux,
                            double * ncrestrict out_dist,
                            std::size_t n, double slab_halfthickness );

      //Format value with most appropriate unit. The fmtstr is used as in
      //NCFmt.hh:
      void fmtBestUnit( std::ostream& os, Length,
                        const char* fmtstr = nullptr );

      //returns macroscocopic XS, or attenuation coefficient, in units of
      //1/m. NB, using nc_as_const for pre C++17:
      inline constexpr double macroXS( NumberDensity nd, CrossSect xs ) noexcept
      {
        return 100.0 * nc_as_const(nd).dbl() * nc_as_const(xs).dbl();
      }

      inline double fast_sqrt_clippos( double x ) noexcept {
        //Using this in a loop does not actually allow for vectorization,
        //https://gcc.gnu.org/bugzilla/show_bug.cgi?id=91645 . However, perhaps
        //this will be solved in GCC 13/14 or is already solved in clang/msvc?
        //Further investigations needed.
#if defined(__clang__) || defined(__GNUC__)
        return __builtin_sqrt(__builtin_fmax(0.0,x));
#else
        return std::sqrt(ncmax(0.0,x));
#endif
      }

      // Parse MiniMC scenario string.
      // More info: https://github.com/mctools/ncrystal/wiki/minimc_scenario
      //            ncrystal minimc --doc=scenario
      struct ScenarioDecoded {
        std::string geomcfg;
        std::string srccfg;
      };
      ScenarioDecoded decodeScenario( const MatCfg&, StrView scenario );

    }
  }
}

#endif

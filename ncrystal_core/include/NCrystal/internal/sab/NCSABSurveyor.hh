#ifndef NCrystal_SABSurveyor_hh
#define NCrystal_SABSurveyor_hh

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

#include "NCrystal/interfaces/NCSABData.hh"
#include "NCrystal/internal/utils/NCSpan.hh"
#include "NCrystal/internal/utils/NCTinyVector.hh"
#include "NCrystal/internal/sab/NCSABIdx.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {

    class SABSurveyor final : private MoveOnly {

      // Class which is used to investigate how a given SAB alpha,beta-grid
      // layout intersects the phase-space curves of various E/kT values.
      //
      // Specifically, it provides a list of cells ordered by when the energy is
      // high enough that the corresponding neutron phasespace is large enough to
      // respectively touch or cover them.
      //
      // This is intended to serve the basis of further processing in order to
      // provide cross sections or samplings.

    public:
      SABSurveyor( const VectD& alphaGrid,
                   const VectD& betaGrid );
      SABSurveyor( const SABData& );//convenience

      using cellidx_t = SABIdx::PackedIndex;

      struct CellInfo final {
        double e_touch;
        double e_cover;
        cellidx_t cellidx;
        bool operator<(const CellInfo& o) const noexcept {
          if ( e_touch != o.e_touch )
            return e_touch < o.e_touch;
          if ( e_cover != o.e_cover )
            return e_cover < o.e_cover;
          return cellidx.val < o.cellidx.val;
        }
      };

      Span<const CellInfo> data() const noexcept { return m_dataSpan; }

    private:
      Span<const CellInfo> m_dataSpan;
      std::unique_ptr<CellInfo[]> m_dataHolder;
    };

    class SABCellSurvey final : private NoCopyMove {

      // Class which is used to locate and classify the various regions within a
      // specific S(alpha,beta) cell which should be integrated separately to
      // find the total contribution to the integral of S(alpha,beta) within a
      // particular neutron phasespace curve - and to both ensure that the
      // integrand is smooth enough for Romberg integration within each region,
      // as well as making it possible to find any rectangular regions if
      // present. The subsequent integration is not performed by this class, as
      // it could depend on S(alpha,beta) interpolation schemes and choices of
      // numerical quadrature algorithms.
      //
      // The resulting region list is sorted, so the regions at upper alpha
      // values appear FIRST in the list, and regions are always consecutive
      // (i.e. alpha_low of the i'th region equals alpha_up of the (i+1)th
      // region).
    public:
      SABCellSurvey( double alpha1, double alpha2,
                     double beta1, double beta2,
                     double E_div_kT );

      struct Region final {
        double alpha_low;
        double alpha_up;
        bool is_bounded_by_betaminus;
        bool is_bounded_by_betaplus;
      };
      static constexpr unsigned nmax_regions = 6;
      using RegionList = TinyVector<Region,nmax_regions>;
      const RegionList& regions() const { return m_regions; }
      void toJSON( std::ostream& ) const;

    private:
      RegionList m_regions;
    };

  }
}

#endif

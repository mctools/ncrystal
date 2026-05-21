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
      SABSurveyor( const VectD& alphaGrid, const VectD& betaGrid );
      SABSurveyor( const SABData& );//convenience

      //Packed cell index (unpack via unpackCellIdx below):
      using cellidx_t = std::uint_fast32_t;

      //Each entry in a CellList is (E/kT,cell index):
      using CellList = std::vector<std::pair<double,cellidx_t>>;

      //Get sorted list of the minimum energy (E/kT) needed before the
      //phase-space reaches ("touches") a given cell:
      const CellList& getTouchList() const noexcept { return m_touch; }

      //Get sorted list of the minimum energy (E/kT) needed before the
      //phase-space completely covers a given cell:
      const CellList& getCoverList() const noexcept { return m_cover; }

      //Unpack cell idx to (ialpha,ibeta):
      template<class TUInt = unsigned>
      static std::pair<TUInt,TUInt> unpackCellIdx( cellidx_t ci )
      {
        static_assert( std::numeric_limits<TUInt>::max()
                       >= std::numeric_limits<std::uint16_t>::max(), "" );
        constexpr cellidx_t mask = 0xFFFFu;
        return { TUInt( ci >> 16 ), TUInt( ci & mask ) };
      }

      struct CellInfo final {
        double e_touch;
        double e_cover;
        cellidx_t cellidx;
        CellInfo( double et, double ec, cellidx_t ci ) noexcept
          : e_touch(et), e_cover(ec), cellidx(ci) {}
        bool operator<(const CellInfo& o) const noexcept {
          if ( e_touch != o.e_touch )
            return e_touch < o.e_touch;
          if ( e_cover != o.e_cover )
            return e_cover < o.e_cover;
          return cellidx < o.cellidx;
        }
      };
      const std::vector<CellInfo>& data() const noexcept { return m_data; }


    private:
      CellList m_touch, m_cover;
      std::vector<CellInfo> m_data;//fixme this better? (if so, discard others)
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
      //fixme: something else than smallvector? A fixed array + Span interface?
      static constexpr unsigned nmax_regions = 4;
      using RegionList = SmallVector<Region,nmax_regions>;
      const RegionList& regions() { return m_regions; }
      void toJSON( std::ostream& ) const;
    private:
      RegionList m_regions;
    };

  }
}

#endif

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

    private:
      CellList m_touch, m_cover;
    };
  }
}

#endif

#ifndef NCrystal_MMC_Slab_hh
#define NCrystal_MMC_Slab_hh

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

#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/minimc/NCMMC_Basket.hh"
#include "NCrystal/internal/minimc/NCMMC_Utils.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    class Slab {
      double m_dz;
    public:
      Slab( Length dz )
        : m_dz(dz.get())
      {
        static_assert( Length::meter == 1.0, "" );
        nc_assert_always( dz.dbl() > 0.0 );
        nc_assert_always( dz.dbl() < 1e199 );
      }

      void toJSONDecodedItems(std::ostream& os) const
      {
        os << "\"name\":\"slab\",\"dz\":";
        streamJSON(os,m_dz);
      }

      void toCfgString(std::ostream& os) const
      {
        os << "slab;dz="<<fmt(m_dz);
      }

      bool pointIsInside( const Vector& v ) const
      {
        return ncabs(v.z()) <= m_dz;
      }

      void distToVolumeEntry( const NeutronBasket& nb,
                              Span<double> tgt,
                              std::size_t offset ) const
      {
        nc_assert( tgt.size() >= nb.nused);
        nc_assert( offset < nb.nused);
        Utils::distToSlabEntry( nb.z + offset, nb.uz + offset,
                                tgt.data() + offset, nb.nused - offset,
                                m_dz );
      }

      void distToVolumeExit( const NeutronBasket& nb,
                             Span<double> tgt ) const
      {
        nc_assert( tgt.size() >= nb.nused);
        Utils::distToSlabExit( nb.z, nb.uz, tgt.data(), nb.nused, m_dz );
      }
    };

  }
}

#endif

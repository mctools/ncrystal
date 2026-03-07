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

#include "NCMMC_BasketSrcFiller.hh"

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace detail {
      namespace {
        void propagateDistanceImpl( double * ncrestrict x,
                                    double * ncrestrict y,
                                    double * ncrestrict z,
                                    const double * ncrestrict ux,
                                    const double * ncrestrict uy,
                                    const double * ncrestrict uz,
                                    const double * ncrestrict dist,
                                    std::size_t n ) ncnoexceptndebug
        {
          for ( std::size_t i = 0; i < n; ++i )
            x[i] += dist[i] * ux[i];
          for ( std::size_t i = 0; i < n; ++i )
            y[i] += dist[i] * uy[i];
          for ( std::size_t i = 0; i < n; ++i )
            z[i] += dist[i] * uz[i];
        }
      }
    }
  }
}

void NC::MiniMC::detail::propagateDistance( NeutronBasket& nb,
                                            std::size_t basket_offset,
                                            const BasketValBufDbl& distances ) ncnoexceptndebug
{
  nc_assert( basket_offset <= nb.size());
  auto& f = nb.fields;
  propagateDistanceImpl( f.x.data + basket_offset,
                         f.y.data + basket_offset,
                         f.z.data + basket_offset,
                         f.ux.data + basket_offset,
                         f.uy.data + basket_offset,
                         f.uz.data + basket_offset,
                         distances.data + basket_offset,
                         nb.size() - basket_offset );
}

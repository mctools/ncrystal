#ifndef NCrystal_MMC_BasketUtils_hh
#define NCrystal_MMC_BasketUtils_hh

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

#include <cstring>
#include "NCrystal/internal/minimc/NCMMC_Defs.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace BasketUtils {

      //type-safe memcopy:
      template<class TData>
      inline void memcpydata( TData* ncrestrict dst,
                              const TData* ncrestrict src,
                              std::size_t n ) ncnoexceptndebug
      {
        assert(n>0);
        nc_assert( n > 0 );
        nc_assert( dst!=nullptr );
        nc_assert( src!=nullptr );
        nc_assert( ( dst >= src+n ) || ( src >= dst+n ) );//no overlap
        std::memcpy( (void*)dst, (void*)src, n*sizeof(TData));
      }

      inline void basketfields_set( NeutronBasketFields& dst,
                                    const NeutronBasketFields& src,
                                    std::size_t i_o,
                                    std::size_t i )
      {
        dst.x.data[i] = src.x[i_o];
        dst.y.data[i] = src.y[i_o];
        dst.z.data[i] = src.z[i_o];
        dst.ux.data[i] = src.ux[i_o];
        dst.uy.data[i] = src.uy[i_o];
        dst.uz.data[i] = src.uz[i_o];
        dst.w.data[i] = src.w[i_o];
        dst.ekin.data[i] = src.ekin[i_o];
      }

      inline void basketfields_setrange( NeutronBasketFields& dst,
                                         const NeutronBasketFields& src,
                                         std::size_t i_o,
                                         std::size_t i,
                                         std::size_t n ) ncnoexceptndebug
      {
        nc_assert( &dst != &src );//memcpy, not memmov
        nc_assert( n>=1 );
        nc_assert( i + n <= basket_N );
        nc_assert( i_o + n <= basket_N );
        memcpydata<double>( dst.x.data + i, src.x.data + i_o, n);
        memcpydata<double>( dst.y.data + i, src.y.data + i_o, n);
        memcpydata<double>( dst.z.data + i, src.z.data + i_o, n);
        memcpydata<double>( dst.ux.data + i, src.ux.data + i_o, n);
        memcpydata<double>( dst.uy.data + i, src.uy.data + i_o, n);
        memcpydata<double>( dst.uz.data + i, src.uz.data + i_o, n);
        memcpydata<double>( dst.w.data + i, src.w.data + i_o, n);
        memcpydata<double>( dst.ekin.data + i, src.ekin.data + i_o, n);
      }

      inline void
      basketfields_validateIfDbg( const NeutronBasketFields& f,
                                  std::size_t nused ) ncnoexceptndebug
      {
        (void)f;
        (void)nused;
#ifndef NDEBUG
        nc_assert( nused <= basket_N );
        for ( std::size_t i = 0; i < nused; ++i ) {
          nc_assert( std::isfinite(f.x[i])
                     && std::isfinite(f.y[i]) && std::isfinite(f.z[i]) );
          nc_assert( std::isfinite(f.ux[i])
                     && std::isfinite(f.uy[i]) && std::isfinite(f.uz[i]) );
          nc_assert( std::isfinite(f.w[i]) && std::isfinite(f.ekin[i]) );
          nc_assert( std::fabs(f.ux[i]*f.ux[i]+f.uy[i]*f.uy[i]
                               +f.uz[i]*f.uz[i]-1.0)<1e-9 );
        }
#endif
      }

      inline void basket_validateIfDbg( const NeutronBasket& b )
      {
        basketfields_validateIfDbg(b.fields,b.nused);
      }

      inline NeutronEnergy&
      ekin_obj( NeutronBasketFields& f, std::size_t i ) ncnoexceptndebug
      {
        nc_assert( i < basket_N );
        static_assert(sizeof(NeutronEnergy)==sizeof(double),"");
        return *reinterpret_cast<NeutronEnergy*>( &f.ekin.data[i] );
      }

      inline const NeutronEnergy&
      ekin_obj( const NeutronBasketFields& f, std::size_t i ) ncnoexceptndebug
      {
        nc_assert( i < basket_N );
        static_assert(sizeof(NeutronEnergy)==sizeof(double),"");
        return *reinterpret_cast<const NeutronEnergy*>( &f.ekin.data[i] );
      }

      inline NeutronEnergy&
      ekin_obj( NeutronBasket& nb, std::size_t i ) ncnoexceptndebug
      {
        nc_assert( i < nb.size() );
        return ekin_obj( nb.fields, i );
      }

      inline const NeutronEnergy&
      ekin_obj( const NeutronBasket& nb, std::size_t i ) ncnoexceptndebug
      {
        nc_assert( i < nb.size() );
        return ekin_obj( nb.fields, i );
      }

      ncnodiscard17 inline NeutronDirection
      dir_obj( const NeutronBasketFields& f, std::size_t i ) ncnoexceptndebug
      {
        nc_assert( i < basket_N );
        return NeutronDirection{ f.ux[i], f.uy[i], f.uz[i] };
      }

      ncnodiscard17 inline NeutronDirection
      dir_obj( const NeutronBasket& nb, std::size_t i ) ncnoexceptndebug
      {
        nc_assert( i < nb.size() );
        return dir_obj( nb.fields, i );
      }

    }
  }
}

#endif

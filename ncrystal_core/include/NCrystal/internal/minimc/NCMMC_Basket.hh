#ifndef NCrystal_MMC_Basket_hh
#define NCrystal_MMC_Basket_hh

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

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    //For efficiency, handle larger number of neutrons at once, with each field
    //in a separate array. This has several advantages, mainly:
    //
    // 1) Call overhead (e.g. with virtual functions) neglible when treating
    //    many particles at once.
    // 2) Easier to achieve simd vectorisation, hot caches, etc.

    namespace detail {
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
    }

    struct NeutronBasketFields final : NoCopyMove {
      BasketValBufDbl x, y, z, ux, uy, uz, ekin, w;
      //fixme: standalone functions in basketutils instead? The we can move
      //NeutronBasketFields+NeutronBasket to Defs.hh

      void set( const NeutronBasketFields& o,
                std::size_t i_o,
                std::size_t i )//fixme switch arg order
      {
        x[i] = o.x[i_o];
        y[i] = o.y[i_o];
        z[i] = o.z[i_o];
        ux[i] = o.ux[i_o];
        uy[i] = o.uy[i_o];
        uz[i] = o.uz[i_o];
        w[i] = o.w[i_o];
        ekin[i] = o.ekin[i_o];
      }

      void setrange( const NeutronBasketFields& o,
                     std::size_t i_o,
                     std::size_t i,
                     std::size_t n ) ncnoexceptndebug
      {
        nc_assert( this != &o );//memcpy, not memmov
        nc_assert( n>=1 );
        nc_assert( i + n <= basket_N );
        nc_assert( i_o + n <= basket_N );
        detail::memcpydata<double>( x.data + i, o.x.data + i_o, n);
        detail::memcpydata<double>( y.data + i, o.y.data + i_o, n);
        detail::memcpydata<double>( z.data + i, o.z.data + i_o, n);
        detail::memcpydata<double>( ux.data + i, o.ux.data + i_o, n);
        detail::memcpydata<double>( uy.data + i, o.uy.data + i_o, n);
        detail::memcpydata<double>( uz.data + i, o.uz.data + i_o, n);
        detail::memcpydata<double>( w.data + i, o.w.data + i_o, n);
        detail::memcpydata<double>( ekin.data + i, o.ekin.data + i_o, n);
      }


      void validateIfDbg( std::size_t nused ) const ncnoexceptndebug
      {
#ifndef NDEBUG
        nc_assert( nused <= basket_N );
        for ( std::size_t i = 0; i < nused; ++i ) {
          nc_assert( std::isfinite(x[i])
                     && std::isfinite(y[i]) && std::isfinite(z[i]) );
          nc_assert( std::isfinite(ux[i])
                     && std::isfinite(uy[i]) && std::isfinite(uz[i]) );
          nc_assert( std::isfinite(w[i]) && std::isfinite(ekin[i]) );
          nc_assert( ncabs(ux[i]*ux[i]+uy[i]*uy[i]+uz[i]*uz[i]-1.0)<1e-9 );
        }
#endif
      }
    };

    struct NeutronBasket final : NoCopyMove {
      NeutronBasketFields fields;
      std::size_t nused = 0;

      void validateIfDbg() const ncnoexceptndebug
      {
        fields.validateIfDbg(nused);
      }
      ncnodiscard17 NeutronEnergy& ekin_obj( std::size_t i ) noexcept
      {
        static_assert(sizeof(NeutronEnergy)==sizeof(double),"");
        return *reinterpret_cast<NeutronEnergy*>( &fields.ekin.data[i] );
      }
      ncnodiscard17 const NeutronEnergy& ekin_obj( std::size_t i ) const noexcept
      {
        static_assert(sizeof(NeutronEnergy)==sizeof(double),"");
        return *reinterpret_cast<const NeutronEnergy*>( &fields.ekin.data[i] );
      }

      ncnodiscard17 NeutronDirection dir_obj( std::size_t i ) const noexcept
      {
        return NeutronDirection{ fields.ux[i], fields.uy[i], fields.uz[i] };
      }

      constexpr bool full() const noexcept { return nused==basket_N; }
      constexpr bool empty() const noexcept { return nused==0; }
      constexpr std::size_t size() const noexcept { return nused; }

      void copyEntryFromOther( const NeutronBasket& o,
                               std::size_t i_o,
                               std::size_t i ) ncnoexceptndebug
      {
        fields.set(o.fields,i_o,i);
      }

      void appendEntriesFromOther( const NeutronBasket& other,
                                   std::size_t i_o,
                                   std::size_t n ) ncnoexceptndebug
      {
        nc_assert( nused + n <= basket_N );
        nc_assert( i_o + n <= other.size() );
        fields.setrange( other.fields, i_o, nused, n );
        nused += n;

      }

    };
  }
}

#endif

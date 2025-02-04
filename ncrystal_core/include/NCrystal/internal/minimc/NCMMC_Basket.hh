#ifndef NCrystal_MMC_Basket_hh
#define NCrystal_MMC_Basket_hh

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

#include "NCrystal/internal/minimc/NCMMC_Defs.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    namespace detail {
      template<class TData>
      inline void memcpydata( TData* ncrestrict dst,
                              const TData* ncrestrict src,
                              std::size_t n ) noexcept
      {
        std::memcpy( (void*)dst, (void*)src, n*sizeof(TData));
      }
    }

    //For efficiency, handle larger number of neutrons at once, with each field
    //in a separate array. This has several advantages:
    //
    // 1) Call overhead (e.g. with virtual functions) neglible when treating
    //    many particles at once.
    // 2) Easier to achieve simd vectorisation, hot caches, etc.
    // 3) Easy to export the particles in a basket to a Python layer, since we
    //    can simply create a Numpy array-view of each data array.

    struct NeutronBasket : NoCopyMove {
      static constexpr auto N = basket_N;
      double x[N];
      double y[N];
      double z[N];
      double ux[N];
      double uy[N];
      double uz[N];
      double w[N];
      double ekin[N];
      //TODO: Time as well? (remember to update it in all propagations if so)
      std::size_t nused = 0;

      ncnodiscard17 NeutronEnergy& ekin_obj( std::size_t i ) noexcept
      {
        static_assert(sizeof(NeutronEnergy)==sizeof(double),"");
        return *reinterpret_cast<NeutronEnergy*>( &ekin[i] );
      }
      ncnodiscard17 const NeutronEnergy& ekin_obj( std::size_t i ) const noexcept
      {
        static_assert(sizeof(NeutronEnergy)==sizeof(double),"");
        return *reinterpret_cast<const NeutronEnergy*>( &ekin[i] );
      }

      ncnodiscard17 NeutronDirection dir_obj( std::size_t i ) const noexcept { return NeutronDirection{ ux[i], uy[i], uz[i] }; }

      constexpr bool full() const noexcept { return nused==N; }
      constexpr bool empty() const noexcept { return nused==0; }
      constexpr std::size_t size() const noexcept { return nused; }

      void copyEntryFromOther( const NeutronBasket& o, std::size_t i_o, std::size_t i ) ncnoexceptndebug
      {
        nc_assert( this != &o );
        x[i] = o.x[i_o];
        y[i] = o.y[i_o];
        z[i] = o.z[i_o];
        ux[i] = o.ux[i_o];
        uy[i] = o.uy[i_o];
        uz[i] = o.uz[i_o];
        w[i] = o.w[i_o];
        ekin[i] = o.ekin[i_o];
      }

      void appendEntriesFromOther( const NeutronBasket& o, std::size_t i_o, std::size_t n ) ncnoexceptndebug
      {
        nc_assert( this != &o );
        nc_assert( n>=1 );
        nc_assert( this->size() + n <= basket_N );
        nc_assert( i_o + n <= o.size() );
        const std::size_t p = this->size();
        this->nused += n;
        detail::memcpydata<double>( this->x + p, o.x + i_o, n);
        detail::memcpydata<double>( this->y + p, o.y + i_o, n);
        detail::memcpydata<double>( this->z + p, o.z + i_o, n);
        detail::memcpydata<double>( this->ux + p, o.ux + i_o, n);
        detail::memcpydata<double>( this->uy + p, o.uy + i_o, n);
        detail::memcpydata<double>( this->uz + p, o.uz + i_o, n);
        detail::memcpydata<double>( this->w + p, o.w + i_o, n);
        detail::memcpydata<double>( this->ekin + p, o.ekin + i_o, n);
      }

    };

    //We might need extra fields during a simulation, for e.g. caches
    //(i.e. cross sections) or statistics collection
    //(i.e. number-of-scatterings). Those caches will be templated (TCache) and
    //a basket of those will be kept with the normal neutron basket in a
    //CachedNeutronBasket object. To ensure we add no significant memory usage
    //if TCache has no data members, we use the empty-base-object optimisation
    //to add the cache data.

    template<class TCache>
    class CachedNeutronBasket {
    public:
      static_assert( std::is_standard_layout<TCache>::value, "" );
      CachedNeutronBasket() noexcept {}
      NeutronBasket neutrons;
      TCache cache;

      constexpr bool full() const noexcept { return neutrons.full(); }
      constexpr bool empty() const noexcept { return neutrons.empty(); }
      constexpr std::size_t size() const noexcept { return neutrons.size(); }

      void copyEntryFromOther( const CachedNeutronBasket& o, std::size_t i_o, std::size_t i ) ncnoexceptndebug
      {
        nc_assert( i_o < o.size() );
        nc_assert( i < this->size() );
        neutrons.copyEntryFromOther( o.neutrons, i_o, i );
        if (!std::is_empty<TCache>::value)
          this->cache.copyEntryFromOther( o.cache, i_o, i );
      }

      void copyEntry( std::size_t i_tgt, std::size_t i_src ) ncnoexceptndebug
      {
        this->copyEntryFromOther( *this, i_src, i_tgt );
      }

      std::size_t appendEntryFromOther( const CachedNeutronBasket& o, std::size_t i_o ) ncnoexceptndebug
      {
        nc_assert(!this->full());
        std::size_t i = neutrons.nused++;
        this->copyEntryFromOther( o, i_o, i );
        return i;
      }

      void appendEntriesFromOther( const CachedNeutronBasket& o, std::size_t i_o, std::size_t n ) ncnoexceptndebug
      {
        nc_assert( this != &o );
        nc_assert( n>=1 );
        nc_assert( this->size() + n <= NeutronBasket::N );
        nc_assert( i_o + n <= o.size() );

        const std::size_t p = this->size();
        neutrons.appendEntriesFromOther( o.neutrons, i_o, n );//nb: call updates this->neutrons.nused value.
        if (!std::is_empty<TCache>::value)
          this->cache.copyEntriesFromOther( o.cache, p, i_o, n );
      }


    };

  }
}

#endif

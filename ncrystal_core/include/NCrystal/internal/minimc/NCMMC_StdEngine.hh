#ifndef NCrystal_MMC_StdEngine_hh
#define NCrystal_MMC_StdEngine_hh

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

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// A simulation engine which employs just a few "tricks" for variance-        //
// reduction:                                                                 //
//                                                                            //
//   * After an elastic scattering in an isotropic material, the cross        //
//     section is unchanged.                                                  //
//                                                                            //
//   * We always force scatterings and always emit transmitted neutrons       //
//     unconditionally (thus allowing more vectorisation!), but we use        //
//     russian roulette before adding the scattered states to the pending     //
//     stacks for further processing (otherwise the models would blow up and  //
//     spend infinite time on unprobable paths).                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/minimc/NCMMC_RunSim.hh"
#include "NCrystal/internal/minimc/NCMMC_Basket.hh"
#include "NCrystal/internal/minimc/NCMMC_BasketMgr.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    struct DPCacheData {
      using data_t = DPCacheData;
      using idx_t = std::size_t;
      int nscat[NeutronBasket::N];
      bool sawinelas[NeutronBasket::N];
      double scatxsval[NeutronBasket::N];//if we already know scat-xs (<0.0 means unknown)

      void markAsMissedTarget( std::size_t i ) noexcept { this->nscat[i] = -1; }

      void markScatteredElastic( std::size_t i ) noexcept
      {
        ++this->nscat[i];
      }

      void markScatteredInelastic( std::size_t i ) noexcept {
        ++this->nscat[i];
        this->sawinelas[i] = true;
      }

      void init( std::size_t i ) noexcept
      {
        nscat[i] = 0;
        sawinelas[i] = false;
        scatxsval[i] = -1.0;
      }

      void copyEntryFromOther( const data_t& o, std::size_t i_o, std::size_t i ) noexcept
      {
        nscat[i] = o.nscat[i_o];
        sawinelas[i] = o.sawinelas[i_o];
        scatxsval[i] = o.scatxsval[i_o];

      }
      void copyEntriesFromOther( const data_t& o, std::size_t i,
                                 std::size_t i_o, std::size_t n ) ncnoexceptndebug
      {
        nc_assert( this != &o );
        nc_assert( i+n <= NeutronBasket::N );
        nc_assert( i_o+n <= NeutronBasket::N );
        detail::memcpydata<int>( nscat + i,o.nscat + i_o,n );
        detail::memcpydata<bool>( sawinelas + i, o.sawinelas + i_o, n );
        detail::memcpydata<double>( scatxsval + i, o.scatxsval + i_o, n );
      }

    };
    static_assert( std::is_standard_layout<DPCacheData>::value, "" );

    class StdEngine final {
    public:
      using Cache = DPCacheData;
      using basket_t = CachedNeutronBasket<Cache>;
      using basketmgr_t = BasketMgr<basket_t>;
      using basket_holder_t = typename basketmgr_t::basket_holder_t;
      static constexpr auto local_basket_max_poolsize = 4;
      using heapmempool_t = HeapMemPool<basket_holder_t::heapmem_t,
                                        (local_basket_max_poolsize==0
                                         ?1:local_basket_max_poolsize)>;

      using matdef_t = MatDef;
    private:
      StdEngineOptions m_opt;
      double m_opt_roulette_survivor_boost;
      matdef_t m_mat;
      CachePtr m_sct_cacheptr;
      CachePtr m_abs_cacheptr;

      heapmempool_t m_mempool;
      basket_holder_t allocateBasket( basketmgr_t& mgr )
      {
        if ( m_mempool.empty() )
          return mgr.allocateBasket();
        return { m_mempool.allocate() };
      }
      void deallocateBasket( basketmgr_t& mgr, basket_holder_t&& bh )
      {
        if ( m_mempool.size() == local_basket_max_poolsize )
          mgr.deallocateBasket( std::move(bh) );
        else
          m_mempool.deallocate( bh.stealMemory() );
      }

      //We need a few buffers (careful with memory usage!):
      double m_buf_disttoexit[basket_N];
      double m_buf_xs_abs[basket_N];
      double m_buf_ptransm[basket_N];
      double m_buf_disttoscat[basket_N];

    public:

      StdEngine( matdef_t md, StdEngineOptions opts = {} );

      shared_obj<StdEngine> clone_so()
      {
        return makeSO<StdEngine>( m_mat, m_opt );
      }

      using resultfct_t = std::function<void(const basket_t&)>;
      void advanceSimulation( RNG& rng,
                              const Geometry& geom,
                              basket_holder_t&& inbasket_holder,
                              basketmgr_t& mgr,
                              const resultfct_t& resultFct );
    };

  }
}

#endif

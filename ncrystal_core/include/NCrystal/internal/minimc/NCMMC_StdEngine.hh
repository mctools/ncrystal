#ifndef NCrystal_MMC_StdEngine_hh
#define NCrystal_MMC_StdEngine_hh

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
#include "NCrystal/internal/minimc/NCMMC_UBView.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    struct DPCacheData {
      using data_t = DPCacheData;
      using idx_t = std::size_t;
      BasketValBufInt nscat;
      BasketValBufBool sawinelas;
      BasketValBufDbl scatxsval;//if already know scat-xs (<0.0 means unknown)

      //fixme: Check which of these we really need?
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
        nc_assert( i+n <= basket_N );
        nc_assert( i_o+n <= basket_N );
        detail::memcpydata<int>( nscat.data + i,o.nscat.data + i_o,n );
        detail::memcpydata<bool>( sawinelas.data + i,
                                  o.sawinelas.data + i_o, n );
        detail::memcpydata<double>( scatxsval.data + i,
                                    o.scatxsval.data + i_o, n );
      }

    };

    static_assert( std::is_standard_layout<DPCacheData>::value, "" );

    class BasketView_StdEngine final : public BasketView {
      //Fixme: We should migrate to a different view mechanism.
      const CachedNeutronBasket<DPCacheData> * m_b;
    public:
      BasketView_StdEngine( const CachedNeutronBasket<DPCacheData> * b )
        : m_b(b) {}
      const NeutronBasket& neutrons() const override { return m_b->neutrons; }
      const BasketValBufInt * nscat() const override { return &m_b->cache.nscat; }
      const BasketValBufBool * sawinelas() const override { return &m_b->cache.sawinelas; }
      const NeutronBasketFields* neutrons_original() const override { return nullptr; }
    };

    class StdEngine final {
    public:
      using cache_t = DPCacheData;
      using basket_t = CachedNeutronBasket<cache_t>;
      using basketmgr_t = BasketMgr<basket_t>;
      using basket_holder_t = typename basketmgr_t::basket_holder_t;
      static constexpr int local_basket_max_poolsize = 4;
      using heapmem_t = typename basket_holder_t::heapmem_t;
      using heapmempool_t = HeapMemPool<heapmem_t,local_basket_max_poolsize>;
      using matdef_t = MatDef;
      using basket_view_t = BasketView_StdEngine;
    private:
      EngineOpts m_opt;
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

      StdEngine( matdef_t md, const EngineOpts& opts = {} );

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

    class SimEngine : NoCopyMove {
    public:
      //Fixme: this interface should go somewhere else.
      // A simulation engine is an engine which is able to move the simulation
      // forward one step at a time, processing a pending basket of neutrons and
      // either discard neutrons, pass them on for additional simulation steps,
      // or serve them up for tallies. This is all done via the basket manager
      // and result callback function, and for efficiency everything is done in
      // terms of baskets of neutrons.
      //
      // The class implementing this will most likely accept options including
      // for geometry, material and variance reduction in its constructor.

      virtual ~SimEngine() = default;

      //Produce a clone (most likely to have a different SimEngine object for
      //each thread in a multithreaded simulation):
      virtual shared_obj<SimEngine> clone() const = 0;

      //Advance the simulation one step. This does not need to be multi-thread
      //safe.
      using TallyFct = std::function<void(const UniversalBasket&)>;
      virtual void step( UniversalBasket,
                         RNG&,
                         UniversalBasketMgr&,
                         const TallyFct& tallyfct ) = 0;
    };

    //Create a std simulation engine through this factory function:
    shared_obj<SimEngine> createStdSimEngine( GeometryPtr,
                                              MatDef,
                                              const EngineOpts& opts = {} );


  }
}

#endif

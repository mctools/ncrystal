#ifndef NCrystal_MMC_SimMgr_hh
#define NCrystal_MMC_SimMgr_hh

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

#include "NCrystal/internal/minimc/NCMMC_UBView.hh"//fixme rename
#include "NCrystal/internal/minimc/NCMMC_Tally.hh"
#include "NCrystal/internal/minimc/NCMMC_EngineOpts.hh"
#include "NCrystal/internal/minimc/NCMMC_CBMgr.hh"
#ifndef NCRYSTAL_DISABLE_THREADS
#  include <thread>
#  include <condition_variable>
#endif

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    //Simulation manager class, responsible for managing worker threads, etc.  .

    //fixme: pimpl everything?

    class SimMgr final : NoCopyMove {
    public:

      SimMgr( const EngineOpts& eopts,
              shared_obj<SimEngine> engine,
              BasketManagementPair basketmgrs,
              shared_obj<TallyMgr> tallymgr,
              Optional<CB::CBMgrInput> callback = NullOpt )
        : m_engineOpts(eopts),
          m_bmgr( std::move( basketmgrs.first ) ),
          m_bprovider( std::move( basketmgrs.second ) ),
          m_engine( std::move(engine) ),
          m_tallymgr( std::move(tallymgr) )
      {
        if ( callback.has_value() )
          m_cbmgr.emplace( std::move(callback.value()) );
      }

      struct LaunchSimReturnVal {
        ParticleCountSum miss;
        ParticleCountSum tallied;
      };
      //returns particles missing geometry
      LaunchSimReturnVal launchSimulations( ThreadCount nthreads,
                                            std::uint64_t seed )
      {
#ifndef NCRYSTAL_DISABLE_THREADS
        return launchSimulationsImpl(nthreads, seed);
#else
        if ( nthreads.get() > 1 )
          NCRYSTAL_WARN("NCrystal installation does not support threads."
                        " Running simulation in single thread and not "
                        "the requested "<<nthreads<<" threads");
        return launchSimulationsImpl( seed );
#endif
      }

    private:
      EngineOpts m_engineOpts;
      shared_obj<UniversalBasketMgr> m_bmgr;
      shared_obj<InputBasketProvider> m_bprovider;
      shared_obj<SimEngine> m_engine;
      shared_obj<TallyMgr> m_tallymgr;
      Optional<CB::CBMgr> m_cbmgr;
      // GeometryPtr m_geom;
      // EngineOpts m_engineOpts;
      // shared_obj<basketmgr_t> m_basketmgr;
      // shared_obj<basket_srcfiller_t> m_srcfiller;
      // shared_obj<TSim> m_sim;
      // shared_obj<TallyMgr> m_tallymgr;
#ifndef NCRYSTAL_DISABLE_THREADS
      SmallVector<std::thread,64> m_workers;
      struct CommonThreadWaitingInfo {
        ThreadCount nthreads = ThreadCount{0};
        std::atomic<ThreadCount::value_type> nthreads_inactive = {0};
        std::mutex mutex;
        std::condition_variable condvar;
      };

      static ParticleCountSum
      doWork( RNG& rng,
              SimEngine& engine,
              UniversalBasketMgr& basketManager,
              InputBasketProvider& basketProvider,
              std::function<void(const UniversalBasket&)>& result_fct,
              CommonThreadWaitingInfo& common,
              bool do_ignoreMiss )
      {
        UniversalBasket basket;
        std::function<void(const UniversalBasket&)> result_fct_srcmiss;
        if ( !do_ignoreMiss )
          result_fct_srcmiss = result_fct;
        ParticleCountSum missStats;

        while( true ) {
          //Get pending basket. If there are no more pending, it might
          //indicate that we are done. But before concluding that, we do first
          //retry after a while, since other threads might have something to
          //provide.
          auto tryFillBasket =  [&basket,
                                 &basketManager,
                                 &basketProvider,
                                 &rng,
                                 &result_fct_srcmiss,
                                 &common,
                                 &missStats]
          {
            if ( !basket.valid() )
              basket = basketProvider.getInputBasket( rng,
                                                      result_fct_srcmiss,
                                                      missStats );
            nc_assert(!basket.valid()||basket.size()>0);
            return basket.valid();
          };

           auto allThreadsInactive = [&common]()
           {
             return common.nthreads_inactive.load() == common.nthreads.get();
           };

          bool we_are_marked_inactive = false;
          if (!tryFillBasket()) {
            //Did not get anything. Mark ourselves inactive and check right away
            //if this is because all threads are now done with all src and
            //pending particles:
            auto nthreads_inactive = common.nthreads_inactive.fetch_add(1) + 1;
            we_are_marked_inactive = true;
            if ( nthreads_inactive == common.nthreads.get() ) {
              if (!tryFillBasket()) {
                //Truly no more src particles!
                std::unique_lock<std::mutex> lock(common.mutex);
                common.condvar.notify_all();
                return missStats;
              }
            }
          }

          while( !basket.valid() ) {
            nc_assert_always( we_are_marked_inactive );
            //Ok, there was no pending baskets available and we have marked
            //ourselves as inactive, but there might be some other active
            //threads. Let us wait for them to potentially give us something
            //(via the common conditions variable):

            std::unique_lock<std::mutex> lock(common.mutex);
            common.condvar.wait( lock, [&tryFillBasket,&allThreadsInactive]()
            {
              return tryFillBasket() || allThreadsInactive();
            });
            if ( !basket.valid()
                 && allThreadsInactive()
                 && !tryFillBasket() ) {
              //make sure other waiting threads might notice
              common.condvar.notify_all();
              return missStats;//Truly no more src particles!
            }
            if ( basket.valid() && we_are_marked_inactive ) {
              common.nthreads_inactive.fetch_sub(1);
              we_are_marked_inactive = false;
            }
          }

          nc_assert_always(!we_are_marked_inactive);
          nc_assert_always(basket.valid());
          nc_assert_always(basket.size()>0);
          //Got something, and we are marked as active. Process it:
          engine.step( std::move(basket), rng, basketManager, result_fct );

          //Probably our call to advanceSimulation has produced some pending
          //baskets. However, before we start notifying other threads, we can
          //have a look ourselves (that way we can avoid messing with the
          //common.mutex/condvar). Of course, if we find nothing ourselves,
          //there won't be anything to notify other threads about.
          if ( tryFillBasket() ) {
            std::unique_lock<std::mutex> lock(common.mutex);
            common.condvar.notify_one();
          }
        }
      }

      LaunchSimReturnVal launchSimulationsImpl( ThreadCount nthreads,
                                                std::uint64_t seed )
      {
        if ( nthreads.indicatesAutoDetect() )
          nthreads = ThreadCount{ std::thread::hardware_concurrency() };
        if ( nthreads.get() < 1 )
          nthreads = ThreadCount{1};
        nc_assert_always( nthreads.get() >= 1 && nthreads.get() < 9999 );
      // EngineOpts m_engineOpts;
      // shared_obj<UniversalBasketMgr> m_bmgr;
      // shared_obj<InputBasketProvider> m_bprovider;
      // shared_obj<SimEngine> m_engine;
      // shared_obj<TallyMgr> m_tallymgr;
        auto bprovider_copy = m_bprovider.getsp();
        auto engine_copy = m_engine.getsp();
        // auto sf_copy = m_srcfiller.getsp();
        // auto sim_copy = m_sim.getsp();
        auto tallymgr_copy = m_tallymgr.getsp();
        CommonThreadWaitingInfo common;
        common.nthreads = nthreads;
        const bool do_ignoreMiss = ( m_engineOpts.ignoreMiss
                                     == EngineOpts::IgnoreMiss::YES );
        std::vector<ParticleCountSum> missStats;
        std::vector<ParticleCountSum> ntalliedStats;
        missStats.resize( nthreads.get() );
        ntalliedStats.resize( nthreads.get() );
        auto rng_next = createBuiltinRNG( seed );
        auto nthrval = nthreads.get();
        nc_assert_always( nthrval >=1 && nthrval <= 99999 );
        for ( auto i : ncrange(nthrval) ) {
          if ( i )
            rng_next = rng_next->createJumped();
          auto bmgr_copy = m_bmgr->cloneMgrForThread().getsp();
          std::shared_ptr<RNGStream> rng = rng_next;
          ParticleCountSum& thread_missStats = missStats.at(i);
          ParticleCountSum& thread_ntalliedStats = ntalliedStats.at(i);
          CB::CBMgr* cbmgrptr = nullptr;
          if ( m_cbmgr.has_value() )
            cbmgrptr = &m_cbmgr.value();
          m_workers.push_back(std::thread([rng,
                                           bmgr_copy,
                                           bprovider_copy,
                                           engine_copy,
                                           tallymgr_copy,
                                           &common,do_ignoreMiss,
                                           &thread_missStats,
                                           &thread_ntalliedStats,
                                           cbmgrptr]()
          {
            NCRYSTAL_DEBUGMMCMSG( "In thread "<<std::this_thread::get_id()
                                  <<" RNG @ "<<(void*)rng.get()
                                  <<" (first val gen: "<<rng->generate()<<")" );
            auto theengine_so = engine_copy->clone();
            SimEngine& theengine = *theengine_so.get();
            TallyPtr tally_so = tallymgr_copy->getIndependentTallyPtr();
            UniversalBasketMgr& bmgr = *bmgr_copy.get();
            InputBasketProvider& bprovider = *bprovider_copy.get();
            auto tallyptr = tally_so.get();

            std::function<void(const UniversalBasket&)> result_fct//fixme: typedef TallyFct
              = [tallyptr,cbmgrptr,
                 &thread_ntalliedStats] (const UniversalBasket& b)
              {
                {
                  nc_assert( b.valid()&& b.neutrons != nullptr );
                  const std::size_t n = b.size();
                  double w = 0.0;
                  const double * it = b.neutrons->fields.w.data;
                  const double * itE = it + n;
                  for ( ; it!=itE; ++it )
                    w += *it;
                  thread_ntalliedStats.count += n;
                  thread_ntalliedStats.weight += w;
                }

                tallyptr->registerResultsUB(b);//fixme remove UB post migration
                if ( cbmgrptr )
                  cbmgrptr->registerData( BasketView_UniversalBasket(&b) );
              };

            RNG * rawrng = rng.get();
            nc_assert(rawrng!=nullptr);
            thread_missStats = doWork( *rawrng, theengine,
                                       bmgr,bprovider,
                                       result_fct, common,
                                       do_ignoreMiss );

            NCRYSTAL_DEBUGMMCMSG("Thread provides results.");
            tallymgr_copy->addResult( std::move(tally_so) );
            NCRYSTAL_DEBUGMMCMSG("Thread ends.");
          }));
        }
        //join:
        for(auto& t : m_workers)
          t.join();
        m_workers.clear();

        //Flush callback buffers:
        if ( m_cbmgr.has_value() )
          m_cbmgr.value().flush();

        //combine stats in return value:
        LaunchSimReturnVal rv;
        for ( auto& e : missStats ) {
          rv.miss.count += e.count;
          rv.miss.weight += e.weight;
        };
        for ( auto& e : ntalliedStats ) {
          rv.tallied.count += e.count;
          rv.tallied.weight += e.weight;
        }
        return rv;
      }
#else
      //Fixme: test non-MT builds in CI.
      //Fixme: Less code duplication between MT and non-MT paths!!!!
      LaunchSimReturnVal launchSimulationsImpl( std::uint64_t seed )
      {
        const bool do_ignoreMiss = ( m_engineOpts.ignoreMiss
                                     == EngineOpts::IgnoreMiss::YES );
        //Single threaded:
        auto rng = createBuiltinRNG( seed );
        TallyPtr tally_so = m_tallymgr->getIndependentTallyPtr();
        ParticleCountSum ntalliedStats;
        CB::CBMgr* cbmgrptr = nullptr;
        if ( m_cbmgr.has_value() )
          cbmgrptr = &m_cbmgr.value();
        auto tallyptr = tally_so.get();


        std::function<void(const UniversalBasket&)>//fixme TallyFct typedef
          result_fct = [tallyptr,
                        &ntalliedStats,
                        cbmgrptr](const UniversalBasket& b)
          {
            {
              nc_assert( b.valid()&& b.neutrons != nullptr );
              const std::size_t n = b.size();
              double w = 0.0;
              const double * it = b.neutrons->fields.w.data;
              const double * itE = it + n;
              for ( ; it!=itE; ++it )
                w += *it;
              ntalliedStats.count += n;
              ntalliedStats.weight += w;
            }
            tallyptr->registerResultsUB(b);//fixme remove UB post migration
            if ( cbmgrptr )
              cbmgrptr->registerData( BasketView_UniversalBasket(&b) );
          };
        std::function<void(const UniversalBasket&)> result_fct_srcmiss(nullptr);
        if ( !do_ignoreMiss )
          result_fct_srcmiss = result_fct;

        ParticleCountSum missStats;
        while ( true ) {
          auto basket = m_bprovider->getInputBasket( rng,
                                                     result_fct_srcmiss,
                                                     missStats );
          if ( !basket.valid() ) {
            //Nothing more to process apparently!
            break;
          }
          nc_assert(basket.size()>0);

          //Got something, and we are marked as active. Process it:
          m_engine->step( std::move(basket), rng, m_bmgr, result_fct );
        }
        m_tallymgr->addResult( std::move(tally_so) );

        //Flush callback buffers:
        if ( m_cbmgr.has_value() )
          m_cbmgr.value().flush();

        LaunchSimReturnVal rv;
        rv.miss = missStats;
        rv.tallied = ntalliedStats;
        return rv;
      }
#endif
    };

  }
}
#endif

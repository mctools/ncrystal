#ifndef NCrystal_MMC_SimMgrMT_hh
#define NCrystal_MMC_SimMgrMT_hh

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

#include "NCrystal/internal/minimc/NCMMC_Geom.hh"
#include "NCrystal/internal/minimc/NCMMC_Basket.hh"
#include "NCrystal/internal/minimc/NCMMC_BasketSrcFiller.hh"
#include "NCrystal/internal/minimc/NCMMC_BasketMgr.hh"
#include "NCrystal/internal/minimc/NCMMC_Source.hh"
#include "NCrystal/internal/minimc/NCMMC_Tally.hh"
#ifndef NCRYSTAL_DISABLE_THREADS
#  include <thread>
#  include <condition_variable>
#endif

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {


    //SimMgr for multi-threaded simulations:
    template<class TSim>
    class SimMgrMT : NoCopyMove {
    public:
      using basket_t = typename TSim::basket_t;
      using basketmgr_t = BasketMgr<basket_t>;
      using basket_holder_t = typename basketmgr_t::basket_holder_t;
      using basket_srcfiller_t = BasketSrcFiller<basket_t>;
      using tally_t = Tally<basket_t>;

      SimMgrMT( GeometryPtr geom,
                SourcePtr src,
                shared_obj<TSim> sim,
                shared_obj<TallyMgr> tallymgr,
                Optional<shared_obj<basketmgr_t>> bm = NullOpt )
        : m_geom(geom),
          m_basketmgr( bm.has_value()
                       ? std::move(bm.value())
                       : makeSO<basketmgr_t>() ),
          m_srcfiller( makeSO<basket_srcfiller_t>(std::move(geom),
                                                  std::move(src),
                                                  m_basketmgr,
                                                  ThreadedUsage::Multi) ),
          m_sim( std::move(sim) ),
          m_tallymgr( std::move(tallymgr) )
      {
      }

      basketmgr_t& basketMgr() { return *m_basketmgr; }
      const basketmgr_t& basketMgr() const { return *m_basketmgr; }
      shared_obj<basketmgr_t> basketMgrPtr() { return m_basketmgr; }
      shared_obj<const basketmgr_t> basketMgrPtr() const { return m_basketmgr; }

      void launchSimulations( ThreadCount nthreads, bool await_end = true )
      {
#ifndef NCRYSTAL_DISABLE_THREADS
        launchSimulationsImpl(nthreads,await_end);
#else
        (void)await_end;
        if ( nthreads.get() > 1 )
          NCRYSTAL_WARN("NCrystal installation does not support threads."
                        " Running simulation in single thread and not "
                        "the requested "<<nthreads<<" threads");
        launchSimulationsImpl();
#endif
      }
      void waitForEnd()
      {
#ifndef NCRYSTAL_DISABLE_THREADS
        waitForEndImpl();
#endif
      }

    private:
      GeometryPtr m_geom;
      shared_obj<basketmgr_t> m_basketmgr;
      shared_obj<basket_srcfiller_t> m_srcfiller;
      shared_obj<TSim> m_sim;
      shared_obj<TallyMgr> m_tallymgr;
#ifndef NCRYSTAL_DISABLE_THREADS
      SmallVector<std::thread,64> m_workers;
      struct CommonThreadWaitingInfo {
        ThreadCount nthreads = ThreadCount{0};
        std::atomic<ThreadCount::value_type> nthreads_inactive = {0};
        std::mutex mutex;
        std::condition_variable condvar;
      };

      static void doWork( RNG& rng,
                          TSim& sim,
                          basket_srcfiller_t& srcfiller,
                          std::function<void(const basket_t&)>& resultfct,
                          CommonThreadWaitingInfo& common )
      {
        basket_holder_t basket{no_init};
        while( true ) {
          //Get pending basket. If there are no more pending, it might
          //indicate that we are done. But before concluding that, we do first
          //retry after a while, since other threads might have something to
          //provide.
          auto tryFillBasket =  [&basket,&srcfiller,&rng,&resultfct,&common]
          {
            if ( !basket.valid() )
              basket = srcfiller.getPendingBasket( common.nthreads, rng, resultfct );
            nc_assert(!basket.valid()||basket.basket().size()>0);
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
                return;
              }
            }
          }

          while( !basket.valid() ) {
            nc_assert_always( we_are_marked_inactive );
            //Ok, there was no pending baskets available and we have marked
            //ourselves as inactive, but there might be some other active
            //threads. Let us wait for them to potentially give us something (via the common conditions variable):
            std::unique_lock<std::mutex> lock(common.mutex);
            common.condvar.wait( lock, [&tryFillBasket,&allThreadsInactive]()
            {
              return tryFillBasket() || allThreadsInactive();
            });
            if ( !basket.valid()
                 && allThreadsInactive()
                 && !tryFillBasket() ) {
              common.condvar.notify_all();//make sure other waiting threads might notice
              return;//Truly no more src particles!
            }
            if ( basket.valid() && we_are_marked_inactive ) {
              common.nthreads_inactive.fetch_sub(1);
              we_are_marked_inactive = false;
            }
          }

          nc_assert_always(!we_are_marked_inactive);
          nc_assert_always(basket.valid());
          nc_assert_always(basket.basket().size()>0);
          //Got something, and we are marked as active. Process it:
          sim.advanceSimulation( rng,
                                 srcfiller.geometry(),
                                 std::move(basket),
                                 srcfiller.basketMgr(),
                                 resultfct );
          nc_assert_always(!basket.valid());//currently, the engine MUST consume
                                            //the basket.

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

      void launchSimulationsImpl( ThreadCount nthreads, bool await_end )
      {
        if ( nthreads.indicatesAutoDetect() )
          nthreads = ThreadCount{ std::thread::hardware_concurrency() };
        if ( nthreads.get() < 1 )
          nthreads = ThreadCount{1};
        nc_assert_always( nthreads.get() >= 1 && nthreads.get() < 9999 );
        auto sf_copy = m_srcfiller.getsp();
        auto sim_copy = m_sim.getsp();
        auto tallymgr_copy = m_tallymgr.getsp();
        CommonThreadWaitingInfo common;
        common.nthreads = nthreads;
        for ( auto i : ncrange(nthreads.get()) ) {
          (void)i;
          m_workers.push_back(std::thread([sf_copy,sim_copy,tallymgr_copy,&common]()
          {
            auto rng = getIndependentRNG();
            NCRYSTAL_DEBUGMMCMSG( "In thread "<<std::this_thread::get_id()
                                  <<" RNG @ "<<(void*)rng.get()
                                  <<" (first val gen: "<<rng->generate()<<")" );
            auto thesim_so = sim_copy->clone_so();
            TSim& thesim = *thesim_so.get();
            TallyPtr thetallyptr = tallymgr_copy->getIndependentTallyPtr();
            basket_srcfiller_t& thesf = *sf_copy.get();

            //We must downcast the tally ptr to tally_t, which depends on the
            //basket type:
            tally_t * downcast_tallyptr = dynamic_cast<tally_t*>(thetallyptr.get());
            nc_assert_always(downcast_tallyptr!=nullptr);
            std::function<void(const basket_t&)> result_collect_function
              = [downcast_tallyptr](const basket_t& b) { return downcast_tallyptr->registerResults(b); };
            doWork(rng,thesim,thesf,result_collect_function,common);
            NCRYSTAL_DEBUGMMCMSG("Thread provides results.");
            tallymgr_copy->addResult( std::move(thetallyptr) );
            NCRYSTAL_DEBUGMMCMSG("Thread ends.");
          }));
        }
        if ( await_end )
          waitForEnd();
      }

      void waitForEndImpl()
      {
        for(auto& t : m_workers)
          t.join();
        m_workers.clear();
      }
#else
      void launchSimulationsImpl()
      {
        //Single threaded:
        auto rng = getRNG();
        auto tallyptr = m_tallymgr->getIndependentTallyPtr();
        //We must downcast the tally ptr to tally_t, which depends on the
        //basket type:
        tally_t * downcast_tallyptr = dynamic_cast<tally_t*>(tallyptr.get());
        nc_assert_always(downcast_tallyptr!=nullptr);
        std::function<void(const basket_t&)> result_collect_fct
          = [downcast_tallyptr](const basket_t& b) { return downcast_tallyptr->registerResults(b); };
        while ( true ) {
          auto basket = m_srcfiller->getPendingBasket( ThreadCount{1}, rng, result_collect_fct );
          if ( !basket.valid() ) {
            //Nothing more to process apparently!
            break;
          }
          m_sim->advanceSimulation( rng,
                                    m_geom,
                                    std::move(basket),
                                    m_basketmgr,
                                    result_collect_fct );
        }
        m_tallymgr->addResult( std::move(tallyptr) );
      }
#endif
    };

  }
}
#endif

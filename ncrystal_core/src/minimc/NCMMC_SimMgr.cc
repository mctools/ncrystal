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

#include "NCrystal/internal/minimc/NCMMC_SimMgr.hh"
#ifndef NCRYSTAL_DISABLE_THREADS
#  include <thread>
#  include <condition_variable>
#endif

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

class NCMMC::SimMgr::Impl
{
  EngineOpts m_engineOpts;
  shared_obj<BasketMgr> m_bmgr;
  shared_obj<SimEngine> m_engine;
  shared_obj<TallyMgr> m_tallymgr;
  Optional<CB::CBMgr> m_cbmgr;
public:

  Impl( const EngineOpts& eopts,
        shared_obj<SimEngine> engine,
        shared_obj<BasketMgr> bmgr,
        shared_obj<TallyMgr> tallymgr,
        Optional<CB::CBMgrInput> callback )
    : m_engineOpts(eopts),
      m_bmgr( std::move( bmgr ) ),
      m_engine( std::move(engine) ),
      m_tallymgr( std::move(tallymgr) )
  {
    if ( callback.has_value() ) {
      std::weak_ptr<BasketMgr> weak_bpv = m_bmgr.getsp();
      std::function<void()> haltSrc = [weak_bpv]()
      {
        //halt the source if ptr is available (obviously no need otherwise)
        auto bpv = weak_bpv.lock();
        if ( bpv != nullptr )
          bpv->haltSource();
      };
      m_cbmgr.emplace( std::move(callback.value()),
                       std::move(haltSrc) );
    }
  }

  NCMMC::SimMgr::LaunchSimReturnVal
  launchSim(ThreadCount, std::uint64_t seed);
};

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace {
      ThreadCount determineNThreads( ThreadCount n )
      {
#ifdef NCRYSTAL_DISABLE_THREADS
        if ( n.get() > 1 ) {
          NCRYSTAL_WARN("NCrystal installation does not support threads."
                        " Running simulation in single thread and not "
                        "the requested "<<n<<" threads");
        }
        n = ThreadCount{ 1 };
#endif
        if ( n.indicatesAutoDetect() ) {
#ifdef NCRYSTAL_DISABLE_THREADS
          n = ThreadCount{ 1 };
#else
          n = ThreadCount{ std::thread::hardware_concurrency() };
#endif
        }
        if ( n.get() < 1 )
          n = ThreadCount{1};
        nc_assert_always( n.get() >= 1 && n.get() < 9999 );
        return n;
      }

      ParticleCountSum doSimWork( RNG& rng,
                                  SimEngine& engine,
                                  shared_obj<BasketMgr> bmgr,
                                  const TallyFct& result_fct,
                                  bool do_ignoreMiss )
      {
        NCRYSTAL_DEBUGMMCMSG("Thread "<<std::this_thread::get_id()
                             <<" doSimWork begins");
        TallyFct result_fct_srcmiss;
        if ( !do_ignoreMiss )
          result_fct_srcmiss = result_fct;
        ParticleCountSum missStats;
        WorkerToken token(bmgr);
        while( true ) {
          NCRYSTAL_DEBUGMMCMSG("Thread "<<std::this_thread::get_id()
                               <<" requestInputBasket");

          Basket basket = bmgr->getInputBasket( token,
                                                rng,
                                                result_fct_srcmiss,
                                                missStats );
          if ( !basket.valid() ) {
            NCRYSTAL_DEBUGMMCMSG("Thread "<<std::this_thread::get_id()
                                 <<" doSimWork returns");
            token.dispose();//best to not just rely on destructor (cleaner in
                            //case of exceptions).
            return missStats;//done
          }
          NCRYSTAL_DEBUGMMCMSG("Thread "<<std::this_thread::get_id()
                               <<" start step");
          engine.step( std::move(basket), rng, result_fct );
          NCRYSTAL_DEBUGMMCMSG("Thread "<<std::this_thread::get_id()
                               <<" step done");

        }
      }
    }
  }
}

NCMMC::SimMgr::SimMgr( const EngineOpts& eopts,
                       shared_obj<SimEngine> engine,
                       shared_obj<BasketMgr> bmgr,
                       shared_obj<TallyMgr> tallymgr,
                       Optional<CB::CBMgrInput> callback )
  : m_impl(new Impl(eopts,
                    std::move(engine),
                    std::move(bmgr),
                    std::move(tallymgr),
                    std::move(callback)))
{
}

NCMMC::SimMgr::~SimMgr()
{
  delete m_impl;
}

NCMMC::SimMgr::LaunchSimReturnVal
NCMMC::SimMgr::launchSimulations( ThreadCount nthreads, std::uint64_t seed )
{
  return m_impl->launchSim( nthreads, seed );
}

NCMMC::SimMgr::LaunchSimReturnVal
NCMMC::SimMgr::Impl::launchSim( ThreadCount nthreads, std::uint64_t seed)
{
  nthreads = determineNThreads(nthreads);
  auto nthrval = nthreads.get();
  nc_assert_always( nthrval >= 1 && nthrval < 9999 );

  auto bmgr_copy = m_bmgr.getsp();
  auto engine_copy = m_engine.getsp();
  auto tallymgr_copy = m_tallymgr.getsp();
  const bool do_ignoreMiss = ( m_engineOpts.ignoreMiss
                               == EngineOpts::IgnoreMiss::YES );
  std::vector<ParticleCountSum> missStats;
  std::vector<ParticleCountSum> ntalliedStats;
  std::vector<std::function<void()>> worker_functions;
  missStats.resize( nthrval );
  ntalliedStats.resize( nthrval );
  worker_functions.resize(nthrval);

  auto rng_next = createBuiltinRNG( seed );

  for ( auto i : ncrange(nthrval) ) {
    if ( i )
      rng_next = rng_next->createJumped();
    std::shared_ptr<RNGStream> rng = rng_next;
    ParticleCountSum& thread_missStats = missStats.at(i);
    ParticleCountSum& thread_ntalliedStats = ntalliedStats.at(i);
    CB::CBMgr* cbmgrptr = nullptr;
    if ( m_cbmgr.has_value() )
      cbmgrptr = &m_cbmgr.value();
    worker_functions.at(i) = [ rng,
                               bmgr_copy,
                               engine_copy,
                               tallymgr_copy,
                               do_ignoreMiss,
                               &thread_missStats,
                               &thread_ntalliedStats,
                               cbmgrptr ]()
    {
      NCRYSTAL_DEBUGMMCMSG( "In thread "<<std::this_thread::get_id()
                            <<" RNG @ "<<(void*)rng.get()
                            <<" (first val gen: "<<rng->generate()<<")" );
      auto theengine_so = engine_copy->clone();
      SimEngine& theengine = *theengine_so.get();
      TallyPtr tally_so = tallymgr_copy->getIndependentTallyPtr();
      auto tallyptr = tally_so.get();
      TallyFct result_fct
        = [tallyptr,cbmgrptr,&thread_ntalliedStats] (const Basket& b)
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
          tallyptr->registerResults(b);
          if ( cbmgrptr )
            cbmgrptr->registerData(b);
        };

      RNG * rawrng = rng.get();
      nc_assert(rawrng!=nullptr);
      thread_missStats = doSimWork( *rawrng, theengine,
                                    bmgr_copy,
                                    result_fct,
                                    do_ignoreMiss );

      NCRYSTAL_DEBUGMMCMSG("Thread provides results.");
      tallymgr_copy->addResult( std::move(tally_so) );
      NCRYSTAL_DEBUGMMCMSG("Thread "<<std::this_thread::get_id()
                           <<" processing ends.");
    };
  }

  if ( worker_functions.size() == 1 ) {
    NCRYSTAL_DEBUGMMCMSG( "Single threaded invocation of worker begin." );
    worker_functions.at(0)();
    NCRYSTAL_DEBUGMMCMSG( "Single threaded invocation of worker done." );
  } else {
#ifdef NCRYSTAL_DISABLE_THREADS
    nc_assert_always(false);
#else
    //Must launch subthreads. We must not let an exception escape from a
    //thread, so instead use std::exception_ptr.
    nc_assert_always(worker_functions.size()>1);
    std::exception_ptr thread_exception;
    std::mutex thread_exception_mtx;
    std::vector<std::thread> workers_threads;
    workers_threads.reserve( worker_functions.size() );
    NCRYSTAL_DEBUGMMCMSG( "In main thread "<<std::this_thread::get_id()
                          <<" about to launch worker threads." );
    for ( auto& workfct : worker_functions ) {
      workers_threads.push_back(std::thread([&workfct,
                                             bmgr_copy,
                                             &thread_exception,
                                             &thread_exception_mtx]()
      {
        std::exception_ptr err_ptr;
        try {
          NCRYSTAL_DEBUGMMCMSG( "Work in thread "<<std::this_thread::get_id()
                                <<" begins");
          workfct();
          NCRYSTAL_DEBUGMMCMSG( "Work in thread "<<std::this_thread::get_id()
                                <<" ends");
        } catch (...) {
          err_ptr = std::current_exception();
        }
        if ( err_ptr ) {
          //First of all, halt everything to wind down other threads:
          NCRYSTAL_DEBUGMMCMSG( "Exception detected in thread "
                                <<std::this_thread::get_id()
                                <<". Halting source.");
          bmgr_copy->haltError();
          //Now move exception to shared location (in case multiple threads
          //store something, we will just rethrow the first of them)
          NCRYSTAL_LOCK_GUARD(thread_exception_mtx);
          if (!thread_exception)
            thread_exception = std::move(err_ptr);
          NCRYSTAL_DEBUGMMCMSG( "Work in thread "<<std::this_thread::get_id()
                                <<" ends after exception");

        }
      }));
    }

    //join threads:
    NCRYSTAL_DEBUGMMCMSG("Thread joining.");
    for(auto& t : workers_threads)
      t.join();
    NCRYSTAL_DEBUGMMCMSG("Thread joining done.");
    worker_functions.clear();
    workers_threads.clear();

    //Rethrow any exception:
    if ( thread_exception ) {
      NCRYSTAL_DEBUGMMCMSG("Rethrowing exception from worker thread in"
                           " main thread.");
      std::rethrow_exception(thread_exception);
      thread_exception = nullptr;
    }
#endif
  }

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

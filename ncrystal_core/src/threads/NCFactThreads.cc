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

#include "NCrystal/threads/NCFactThreads.hh"
namespace NC = NCrystal;

#ifdef NCRYSTAL_DISABLE_THREADS

#include "NCrystal/internal/utils/NCString.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace FactoryThreadPool {
    namespace detail {
      namespace {
        //Threads are never used, but a user configuration is still
        //tracked, for consistency with the normal case:
        std::atomic<bool>& userConfiguredFlag()
        {
          static std::atomic<bool> b( ncgetenv_int64( "FACTORY_THREADS",
                                                      -1 ) >= 0 );
          return b;
        }
      }
      void beginTemporaryThreads( ThreadCount ) {}
      void endTemporaryThreads() {}
    }
  }
}

void NC::FactoryThreadPool::enable( ThreadCount )
{
  detail::userConfiguredFlag().store( true );
}
void NC::FactoryThreadPool::queue( voidfct_t job ) { job(); }
NC::FactoryThreadPool::detail::FactoryJobsHandler
NC::FactoryThreadPool::detail::getFactoryJobsHandler() { return {}; }
NC::ThreadCount NC::FactoryThreadPool::currentThreadCount()
{
  return ThreadCount{ 1 };
}
bool NC::FactoryThreadPool::userConfigured()
{
  return detail::userConfiguredFlag().load();
}
bool NC::FactoryThreadPool::threadsAvailable() { return false; }

#else

#include "NCThreadPool.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace FactoryThreadPool {
    namespace detail {
      namespace {
        struct FJH {
          std::mutex mtx;
          detail::FactoryJobsHandler fjh;
        };

        FJH& getFJH()
        {
          static FJH fjh;
          return fjh;
        }

        void setFJH( detail::FactoryJobsHandler&& fjh )
        {
          auto& db = getFJH();
          NCRYSTAL_LOCK_GUARD(db.mtx);
          db.fjh = std::move(fjh);
        }

        ThreadPool::ThreadPool& getTP() {
          static ThreadPool::ThreadPool tp;
          return tp;
        }

        voidfct_t detail_get_pending_job()
        {
          return getTP().getPendingJob();
        }

        //State of the thread-pool configuration:
        struct State {
          std::mutex mtx;
          std::atomic<unsigned> nthreads_total{1};
          std::atomic<bool> user_configured{false};
          //Set once the env var no longer needs to be checked:
          std::atomic<bool> envvar_done{false};
          unsigned ntemp_users = 0;//protected by mtx
        };

        State& getState()
        {
          static State s;
          return s;
        }

        //Change number of threads, without any bookkeeping of who
        //requested it (must be called with State::mtx locked):
        void setThreadsUnlocked( ThreadCount nthreads )
        {
          if ( nthreads.indicatesAutoDetect() )
            nthreads
              = ThreadCount{ std::thread::hardware_concurrency() };
          const unsigned nt = nthreads.get();
          const unsigned n_extra = nt >= 2 ? nt - 1 : 0;
          setFJH( detail::FactoryJobsHandler{ nullptr, nullptr } );
          getTP().changeNumberOfThreads( n_extra );
          if ( n_extra > 0 )
            setFJH( detail::FactoryJobsHandler{
                ::NC::FactoryThreadPool::queue,
                ::NC::FactoryThreadPool::detail::detail_get_pending_job
              } );
          getState().nthreads_total.store( n_extra + 1 );
        }

        //Apply the NCRYSTAL_FACTORY_THREADS env var if set, unless done
        //already, or if enable(..) was called first (which takes
        //precedence):
        void ensureEnvVarProcessed()
        {
          auto& st = getState();
          if ( st.envvar_done.load() )
            return;
          NCRYSTAL_LOCK_GUARD( st.mtx );
          if ( st.envvar_done.load() )
            return;
          st.envvar_done.store( true );
          std::int64_t n = ncgetenv_int64( "FACTORY_THREADS", -1 );
          if ( n < 0 )
            return;
          st.user_configured.store( true );
          const auto nclamped
            = static_cast<unsigned>( n > 9999 ? 9999 : n );
          setThreadsUnlocked( ThreadCount{ nclamped } );
        }
      }
    }
  }
}

void NC::FactoryThreadPool::enable( ThreadCount nthreads )
{
  auto& st = detail::getState();
  NCRYSTAL_LOCK_GUARD( st.mtx );
  st.envvar_done.store( true );//explicit calls take precedence
  st.user_configured.store( true );
  detail::setThreadsUnlocked( nthreads );
}

void NC::FactoryThreadPool::queue( voidfct_t job )
{
  detail::ensureEnvVarProcessed();
  detail::getTP().queue( std::move(job) );
}

NC::FactoryThreadPool::detail::FactoryJobsHandler
NC::FactoryThreadPool::detail::getFactoryJobsHandler()
{
  ensureEnvVarProcessed();
  auto& db = getFJH();
  NCRYSTAL_LOCK_GUARD(db.mtx);
  FactoryJobsHandler fjh = db.fjh;
  return fjh;
}

NC::ThreadCount NC::FactoryThreadPool::currentThreadCount()
{
  detail::ensureEnvVarProcessed();
  return ThreadCount{ detail::getState().nthreads_total.load() };
}

bool NC::FactoryThreadPool::userConfigured()
{
  detail::ensureEnvVarProcessed();
  return detail::getState().user_configured.load();
}

bool NC::FactoryThreadPool::threadsAvailable() { return true; }

void NC::FactoryThreadPool::detail::beginTemporaryThreads(
  ThreadCount n )
{
  ensureEnvVarProcessed();
  auto& st = getState();
  NCRYSTAL_LOCK_GUARD( st.mtx );
  if ( st.ntemp_users++ == 0 && !st.user_configured.load() )
    setThreadsUnlocked( n );
}

void NC::FactoryThreadPool::detail::endTemporaryThreads()
{
  auto& st = getState();
  NCRYSTAL_LOCK_GUARD( st.mtx );
  nc_assert_always( st.ntemp_users > 0 );
  if ( --st.ntemp_users == 0 && !st.user_configured.load() )
    setThreadsUnlocked( ThreadCount{ 1 } );
}

#endif

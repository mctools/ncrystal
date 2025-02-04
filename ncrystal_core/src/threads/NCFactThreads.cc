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

#include "NCrystal/threads/NCFactThreads.hh"
namespace NC = NCrystal;

#ifdef NCRYSTAL_DISABLE_THREADS

void NC::FactoryThreadPool::enable( ThreadCount ) {}
void NC::FactoryThreadPool::queue( voidfct_t job ) { job(); }
NC::FactoryThreadPool::detail::FactoryJobsHandler
NC::FactoryThreadPool::detail::getFactoryJobsHandler() { return {}; }

#else

#include "NCThreadPool.hh"

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

        std::atomic<bool>& getFactThreadsCalledAB()
        {
          static std::atomic<bool> called(false);
          return called;
        }

        ThreadPool::ThreadPool& getTP() {
          static ThreadPool::ThreadPool tp;
          return tp;
        }

        voidfct_t detail_get_pending_job()
        {
          return getTP().getPendingJob();
        }
      }//end anon namespace

    }//end detail namespace
  }

  namespace detail {
    //NOTICE: This thread-safe function will be called by NCFactImpl.cc and is
    //fwd declared there:
    bool factThreadsEnableCalledExplicitly()
    {
      return FactoryThreadPool::detail::getFactThreadsCalledAB().load();
    }
  }
}

void NC::FactoryThreadPool::enable( ThreadCount nthreads )
{
  if ( nthreads.indicatesAutoDetect() )
    nthreads = ThreadCount{ std::thread::hardware_concurrency() };

  detail::getFactThreadsCalledAB().store(true);
  unsigned n_extra_threads = nthreads.get() >= 2 ? nthreads.get() - 1 : 0;
  {
    detail::setFJH( detail::FactoryJobsHandler{ nullptr, nullptr } );
    //detail::setFactoryJobsHandler(nullptr,nullptr);
    detail::getTP().changeNumberOfThreads( n_extra_threads );
    if ( n_extra_threads > 0 )
      detail::setFJH( detail::FactoryJobsHandler{
          ::NC::FactoryThreadPool::queue,
          ::NC::FactoryThreadPool::detail::detail_get_pending_job
        } );
  }
}

void NC::FactoryThreadPool::queue( voidfct_t job )
{
  detail::getTP().queue( std::move(job) );
}

NC::FactoryThreadPool::detail::FactoryJobsHandler
NC::FactoryThreadPool::detail::getFactoryJobsHandler()
{
  auto& db = getFJH();
  NCRYSTAL_LOCK_GUARD(db.mtx);
  FactoryJobsHandler fjh = db.fjh;
  return fjh;
}

#endif

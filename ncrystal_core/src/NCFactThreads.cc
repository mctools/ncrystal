////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

#include "NCrystal/NCFactThreads.hh"
namespace NC = NCrystal;

#ifdef NCRYSTAL_DISABLE_THREADS

void NC::FactoryThreadPool::enable( ThreadCount ) {}
void NC::FactoryThreadPool::queue( std::function<void()> job ) { job(); }

#else

#include "NCThreadPool.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace detail {
    //NOTICE: This is an internal fwd declared function from NCFactoryJobs.cc:
    using voidfct_t = std::function<void()>;
    void setFactoryJobsHandler( std::function<void(voidfct_t)> job_queuefct,
                                std::function<voidfct_t()> get_pending_job_fct );


    std::atomic<bool>& getFactThreadsCalledAB()
    {
      static std::atomic<bool> called(false);
      return called;
    }

    //NOTICE: This thread-safe function will be called by NCFactImpl.cc and is
    //fwd declared there:
    bool factThreadsEnableCalledExplicitly()
    {
      return getFactThreadsCalledAB().load();
    }

  }

  namespace FactoryThreadPool {
    namespace {
      ThreadPool::ThreadPool& getTP() {
        static ThreadPool::ThreadPool tp;
        return tp;
      }
      using voidfct_t = std::function<void()>;
      voidfct_t detail_get_pending_job()
      {
        return getTP().getPendingJob();
      }
    }
  }
}

void NC::FactoryThreadPool::enable( ThreadCount nthreads )
{
  if ( nthreads.indicatesAutoDetect() )
    nthreads = ThreadCount{ std::thread::hardware_concurrency() };

  detail::getFactThreadsCalledAB().store(false);
  unsigned n_extra_threads = nthreads.get() >= 2 ? nthreads.get() - 1 : 0;
  {
    detail::setFactoryJobsHandler(nullptr,nullptr);
    getTP().changeNumberOfThreads( n_extra_threads );
    if ( n_extra_threads > 0 )
      detail::setFactoryJobsHandler(::NC::FactoryThreadPool::queue,
                                    ::NC::FactoryThreadPool::detail_get_pending_job);
  }
}

void NC::FactoryThreadPool::queue( std::function<void()> job )
{
  getTP().queue( std::move(job) );
}

#endif

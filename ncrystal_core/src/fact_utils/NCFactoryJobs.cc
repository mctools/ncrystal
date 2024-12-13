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

#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#ifndef NCRYSTAL_DISABLE_THREADS
#  include <condition_variable>
#  include <chrono>
#endif

namespace NC = NCrystal;

#ifdef NCRYSTAL_DISABLE_THREADS

NC::FactoryJobs::~FactoryJobs() = default;
NC::FactoryJobs::FactoryJobs() = default;
void NC::FactoryJobs::queueMT( voidfct_t ) {}
void NC::FactoryJobs::waitAllMT() {}
NC::FactoryJobs::voidfct_t NC::FactoryJobs::getGloballyPendingJob()
{
  return {};
}

#else

namespace NCRYSTAL_NAMESPACE {

  struct FactoryJobs::MTImpl {
    unsigned m_unfinishedjobs = 0;
    std::function<void(voidfct_t)> m_job_queuefct;
    std::function<voidfct_t()> m_get_pending_job_fct;
    std::mutex m_mutex;
    std::condition_variable m_condvar;
  };

  namespace detail {
    using voidfct_t = std::function<void()>;
    namespace {
      struct FactJobsQueueAccess final : NoCopyMove {
        std::mutex mtx;
        std::function<void(voidfct_t)> job_queuefct;
        std::function<voidfct_t()> get_pending_job_fct;
      };
      FactJobsQueueAccess& getThreadQueueInfo() {
        static FactJobsQueueAccess db;
        return db;
      }
    }

    //NOTE: This next function is fwd declared in NCFactThreads.cc, and is used
    //to avoid a requirement to link NCFactoryJobs.cc with libpthreads (or
    //whatever), which is useful in the NCrystal development environment:
    void setFactoryJobsHandler( std::function<void(voidfct_t)> job_queuefct,
                                std::function<voidfct_t()> get_pending_job_fct )
    {
      auto& tqi = detail::getThreadQueueInfo();
      NCRYSTAL_LOCK_GUARD(tqi.mtx);
      tqi.job_queuefct = std::move(job_queuefct);
      tqi.get_pending_job_fct = std::move(get_pending_job_fct);
    }
  }
}

NC::FactoryJobs::~FactoryJobs()
{
  delete m_mt;
}

NC::FactoryJobs::FactoryJobs()
{
  std::function<void(voidfct_t)> job_queuefct;
  std::function<voidfct_t()> get_pending_job_fct;
  {
    auto& tqi = detail::getThreadQueueInfo();
    NCRYSTAL_LOCK_GUARD(tqi.mtx);
    job_queuefct = tqi.job_queuefct;
    get_pending_job_fct = tqi.get_pending_job_fct;
  }
  nc_assert( bool(job_queuefct) == bool(get_pending_job_fct) );
  if ( job_queuefct ) {
    m_mt = new FactoryJobs::MTImpl;
    m_mt->m_job_queuefct = std::move(job_queuefct);
    m_mt->m_get_pending_job_fct = std::move(get_pending_job_fct);
  }
}

void NC::FactoryJobs::queueMT( voidfct_t job )
{
  nc_assert( m_mt != nullptr );
  {
    std::unique_lock<std::mutex> lock(m_mt->m_mutex);
    ++m_mt->m_unfinishedjobs;
  }
  MTImpl * mt = m_mt;
  m_mt->m_job_queuefct( [mt,job]()
  {
    job();
    std::unique_lock<std::mutex> lock(mt->m_mutex);
    --(mt->m_unfinishedjobs);
    mt->m_condvar.notify_one();
  });
}

void NC::FactoryJobs::waitAllMT()
{
  nc_assert( m_mt != nullptr );
  while ( true ) {
    {
      std::unique_lock<std::mutex> lock(m_mt->m_mutex);
      if ( m_mt->m_unfinishedjobs == 0 )
        break;
    }

    //We still have some of our associated jobs running in the thread
    //pool. Rather than risking a dead-lock (in case we are ourselves running
    //inside the thread pool and thus consuming a thread already), let us
    //process a job from the thread pool, and then check again.
    nc_assert(m_mt->m_get_pending_job_fct);
    std::function<void()> job = m_mt->m_get_pending_job_fct();
    if (job) {
      job();
      continue;
    }

    //No pending jobs in the global queue, so all of our own associated jobs
    //must have at least started to run. Let us simply wait for them to finish
    //(however, we recheck after a small amount of time has elapsed, in case
    //something else showed up in the global job queue that we could help with
    //rather than simply waiting). We are also aware that we might get spurious
    //wakeups:
    std::unique_lock<std::mutex> lock(m_mt->m_mutex);
    m_mt->m_condvar.wait_for(lock, std::chrono::milliseconds(10));
  }
}

NC::FactoryJobs::voidfct_t NC::FactoryJobs::getGloballyPendingJob()
{
  voidfct_t job;
  auto& tqi = detail::getThreadQueueInfo();
  {
    NCRYSTAL_LOCK_GUARD(tqi.mtx);
    if ( tqi.get_pending_job_fct )
      job = tqi.get_pending_job_fct();
  }
  return job;
}

#endif

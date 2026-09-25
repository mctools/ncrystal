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

#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#ifndef NCRYSTAL_DISABLE_THREADS
#  include "NCrystal/threads/NCFactThreads.hh"
#  include <condition_variable>
#  include <chrono>
#  include <exception>
#endif

namespace NC = NCrystal;

#ifdef NCRYSTAL_DISABLE_THREADS

NC::FactoryJobs::~FactoryJobs() = default;
NC::FactoryJobs::FactoryJobs() = default;
void NC::FactoryJobs::queueMT( voidfct_t ) {}
void NC::FactoryJobs::waitAllMT() {}
NC::voidfct_t NC::FactoryJobs::getGloballyPendingJob()
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
    std::exception_ptr m_exception;//first exception from a job
    void waitAllNoThrow();
  };
}

NC::FactoryJobs::~FactoryJobs()
{
  if ( m_mt ) {
    //Never leave jobs running which might reference objects on the
    //stack of the caller (e.g. when an exception is unwinding the stack
    //before waitAll() was called). Any exception is discarded:
    m_mt->waitAllNoThrow();
    delete m_mt;
  }
}

NC::FactoryJobs::FactoryJobs()
{
  auto fjh = FactoryThreadPool::detail::getFactoryJobsHandler();
  nc_assert( bool(fjh.jobQueueFct) == bool(fjh.getPendingJobFct) );
  if ( fjh.jobQueueFct ) {
    m_mt = new FactoryJobs::MTImpl;
    m_mt->m_job_queuefct = std::move(fjh.jobQueueFct);
    m_mt->m_get_pending_job_fct = std::move(fjh.getPendingJobFct);
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
  try {
    m_mt->m_job_queuefct( [mt,job]()
    {
      //Exceptions must not escape (the job might run in a worker thread),
      //so the first is kept and rethrown by waitAll() in the calling thread:
      std::exception_ptr err;
      try {
        job();
      } catch (...) {
        err = std::current_exception();
      }
      std::unique_lock<std::mutex> lock(mt->m_mutex);
      if ( err && !mt->m_exception )
        mt->m_exception = std::move(err);
      --(mt->m_unfinishedjobs);
      mt->m_condvar.notify_one();
    });
  } catch (...) {
    //The job could not be queued (e.g. a memory allocation failed), so it
    //will never run. Undo the book-keeping, to not wait for it forever:
    std::unique_lock<std::mutex> lock(mt->m_mutex);
    --(mt->m_unfinishedjobs);
    throw;
  }
}

void NC::FactoryJobs::waitAllMT()
{
  nc_assert( m_mt != nullptr );
  m_mt->waitAllNoThrow();
  std::exception_ptr err;
  {
    std::unique_lock<std::mutex> lock(m_mt->m_mutex);
    std::swap( err, m_mt->m_exception );
  }
  if ( err )
    std::rethrow_exception( err );
}

void NC::FactoryJobs::MTImpl::waitAllNoThrow()
{
  while ( true ) {
    {
      std::unique_lock<std::mutex> lock(m_mutex);
      if ( m_unfinishedjobs == 0 )
        break;
    }
    //We still have some of our associated jobs running in the thread
    //pool. Rather than risking a dead-lock (in case we are ourselves
    //running inside the thread pool and thus consuming a thread
    //already), let us process a job from the thread pool, and then check
    //again. NB: Such jobs must never throw (jobs queued via FactoryJobs
    //catch exceptions, and other jobs must not throw, as documented in
    //NCFactThreads.hh).
    nc_assert(m_get_pending_job_fct);
    voidfct_t job = m_get_pending_job_fct();
    if (job) {
      FactoryThreadPool::detail::runJobNoThrow( job );
      continue;
    }
    //No pending jobs in the global queue, so all of our own associated
    //jobs must have at least started to run. Let us simply wait for them
    //to finish (however, we recheck after a small amount of time has
    //elapsed, in case something else showed up in the global job queue
    //that we could help with rather than simply waiting). We are also
    //aware that we might get spurious wakeups:
    std::unique_lock<std::mutex> lock(m_mutex);
    m_condvar.wait_for(lock, std::chrono::milliseconds(10));
  }
}

NC::voidfct_t NC::FactoryJobs::getGloballyPendingJob()
{
  voidfct_t job;
  auto fjh = FactoryThreadPool::detail::getFactoryJobsHandler();
  if ( fjh.getPendingJobFct )
    job = fjh.getPendingJobFct();
  return job;
}

#endif

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

#include "NCThreadPool.hh"
namespace NC = NCrystal;

NC::ThreadPool::ThreadPool::ThreadPool() = default;

#ifdef NCRYSTAL_DISABLE_THREADS

void NC::ThreadPool::ThreadPool::changeNumberOfThreads(unsigned) {}
void NC::ThreadPool::ThreadPool::queue( voidfct_t job) { job(); }
NC::voidfct_t NC::ThreadPool::ThreadPool::getPendingJob() { return {}; }
NC::ThreadPool::ThreadPool::~ThreadPool() {}

#else

NC::ThreadPool::ThreadPool::~ThreadPool()
{
  endAllThreads();
}

void NC::ThreadPool::ThreadPool::changeNumberOfThreads( unsigned nthreads )
{
  //Todo: check that this means that we can change number of threads dynamically
  //(multiple times) as we wish

  //  if ( nthreads == 0 )
  //    nthreads = std::thread::hardware_concurrency();//auto

  std::unique_lock<std::mutex> lock(m_mutex);
  if ( nthreads == m_threads.size() ) {
    //no change
  } else if ( nthreads > m_threads.size() ) {
    m_threads_should_end = false;
    m_threads.reserve(nthreads);
    while ( (unsigned)m_threads.size() < nthreads )
      m_threads.emplace_back(std::thread(&ThreadPool::threadWorkFct,this));
  } else {
    nc_assert( nthreads < m_threads.size() );
    //For simplicity, go to a complete halt, then restart (this is anyway most
    //likely to happen when there are no running jobs):
    lock.unlock();
    this->endAllThreads();
    this->changeNumberOfThreads( nthreads );
  }
}

NC::voidfct_t NC::ThreadPool::ThreadPool::getPendingJob()
{
  std::unique_lock<std::mutex> lock(m_mutex);
  if ( m_jobqueue.empty() )
    return nullptr;
  voidfct_t job = m_jobqueue.front();
  m_jobqueue.pop();
  //NB: Not notifying anyone via m_condvar since we *removed* a job.
  return job;
}

// #include <cxxabi.h>

// namespace {

//   // std::string demangle( const char * symbol )
//   // {
//   //   int status;
//   //   std::size_t nbuf = 0;//ok if too short
//   //   //char * buf = static_cast<char*>(std::malloc(nbuf));//MUST use std::malloc
//   //   char * rawres = abi::__cxa_demangle(symbol, nullptr, &nbuf, &status);
//   //   nc_assert_always(status==0);
//   //   std::string res{ rawres };
//   //   std::free(rawres);
//   //   return res;
//   // }
//   //demangle(job.target_type().name()));

// }

void NC::ThreadPool::ThreadPool::threadWorkFct()
{
  while (true) {
    //Wait for more jobs to be available in the queue OR m_threads_should_end to
    //be set. Note that we want our queue to always run all jobs, even if we are
    //trying to end all threads.
    voidfct_t job;
    std::unique_lock<std::mutex> lock(m_mutex);
    m_condvar.wait( lock, [this]
    {
      return !m_jobqueue.empty() || m_threads_should_end;
    });
    if ( !m_jobqueue.empty() ) {
      //there is a job to run, so run it:
      job = std::move(m_jobqueue.front());
      m_jobqueue.pop();
      lock.unlock();
      job();
    } else {
      nc_assert_always( m_threads_should_end );
      return;//end thread
    }
  }
}

void NC::ThreadPool::ThreadPool::queue(voidfct_t job)
{
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    if ( m_threads_should_end ) {
      //TP not active or is winding down:
      lock.unlock();
      job();
      return;
    }
    m_jobqueue.push(std::move(job));
  }
  m_condvar.notify_one();
}

void NC::ThreadPool::ThreadPool::endAllThreads()
{
  {
    std::unique_lock<std::mutex> lock(m_mutex);
    m_threads_should_end = true;
  }
  m_condvar.notify_all();
  std::unique_lock<std::mutex> lock(m_mutex);
  while ( !m_threads.empty() ) {
    {
      std::thread t = std::move( m_threads.back() );
      m_threads.pop_back();
      lock.unlock();
      t.join();
    }
    lock.lock();
  }
}

#endif


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

#include "NCThreadPool.hh"
#include "NCrystal/threads/NCFactThreads.hh"
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

void NC::ThreadPool::ThreadPool::threadWorkFctTrampoline( void* arg )
{
  static_cast<ThreadPool*>( arg )->threadWorkFct();
}

#ifdef _WIN32

NC::ThreadPool::WorkerThread::WorkerThread( void (*fct)( void* ),
                                            void* arg )
  : m_t( fct, arg )
{
}
void NC::ThreadPool::WorkerThread::join() { m_t.join(); }
NC::ThreadPool::WorkerThread::~WorkerThread() = default;
NC::ThreadPool::WorkerThread::WorkerThread( WorkerThread&& ) noexcept
  = default;
NC::ThreadPool::WorkerThread&
NC::ThreadPool::WorkerThread::operator=( WorkerThread&& ) noexcept
  = default;

#else

namespace NCRYSTAL_NAMESPACE {
  namespace ThreadPool {
    namespace {
      struct WTLaunchData {
        void (*fct)( void* );
        void* arg;
      };
      void* wt_trampoline( void* raw )
      {
        WTLaunchData d = *static_cast<WTLaunchData*>( raw );
        delete static_cast<WTLaunchData*>( raw );
        d.fct( d.arg );
        return nullptr;
      }
    }
  }
}

NC::ThreadPool::WorkerThread::WorkerThread( void (*fct)( void* ),
                                            void* arg )
{
  pthread_attr_t attr;
  if ( pthread_attr_init( &attr ) != 0 )
    throw std::runtime_error("pthread_attr_init failed");
  std::size_t stacksize = 8 * 1024 * 1024;
#  ifdef PTHREAD_STACK_MIN
  //NB: PTHREAD_STACK_MIN might expand to a sysconf call (glibc>=2.34):
  const std::size_t psm
    = static_cast<std::size_t>( PTHREAD_STACK_MIN );
  if ( stacksize < psm )
    stacksize = psm;
#  endif
  pthread_attr_setstacksize( &attr, stacksize );
  auto d = new WTLaunchData{ fct, arg };
  const int ec = pthread_create( &m_t, &attr, &wt_trampoline, d );
  pthread_attr_destroy( &attr );
  if ( ec != 0 ) {
    delete d;
    throw std::runtime_error("pthread_create failed");
  }
  m_joinable = true;
}

void NC::ThreadPool::WorkerThread::join()
{
  nc_assert_always( m_joinable );
  pthread_join( m_t, nullptr );
  m_joinable = false;
}

NC::ThreadPool::WorkerThread::~WorkerThread()
{
  //Mirror std::thread semantics loosely, but never terminate(): a
  //still-joinable thread here would be a logic error caught in debug
  //builds:
  nc_assert( !m_joinable );
}

NC::ThreadPool::WorkerThread::WorkerThread( WorkerThread&& o ) noexcept
  : m_t( o.m_t ), m_joinable( o.m_joinable )
{
  o.m_joinable = false;
}

NC::ThreadPool::WorkerThread&
NC::ThreadPool::WorkerThread::operator=( WorkerThread&& o ) noexcept
{
  nc_assert( !m_joinable );
  m_t = o.m_t;
  m_joinable = o.m_joinable;
  o.m_joinable = false;
  return *this;
}

#endif

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
      m_threads.emplace_back( &ThreadPool::threadWorkFctTrampoline, this );
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
      FactoryThreadPool::detail::runJobNoThrow( job );
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
      WorkerThread t = std::move( m_threads.back() );
      m_threads.pop_back();
      lock.unlock();
      t.join();
    }
    lock.lock();
  }
}

#endif


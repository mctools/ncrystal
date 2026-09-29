#ifndef NCrystal_ThreadPool_hh
#define NCrystal_ThreadPool_hh

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

#include "NCrystal/core/NCDefs.hh"

#ifndef NCRYSTAL_DISABLE_THREADS
#  include <thread>
#  include <condition_variable>
#  include <queue>
#  ifndef _WIN32
#    include <pthread.h>
#  endif
#endif

namespace NCRYSTAL_NAMESPACE {

  namespace ThreadPool {

#ifndef NCRYSTAL_DISABLE_THREADS
    //Worker thread with an explicit 8MB stack on POSIX: std::thread
    //offers no stack-size control, and e.g. macOS gives secondary
    //threads just 512kB by default (vs. glibc's usual 8MB), observed
    //(as CI SIGBUS crashes) to be too small for factory jobs loading
    //complex materials. On Windows std::thread is kept (1MB default,
    //no issues observed).
    class WorkerThread final : private MoveOnly {
    public:
      WorkerThread( void (*fct)( void* ), void* arg );
      void join();
      ~WorkerThread();
      WorkerThread( WorkerThread&& ) noexcept;
      WorkerThread& operator=( WorkerThread&& ) noexcept;
    private:
#  ifdef _WIN32
      std::thread m_t;
#  else
      pthread_t m_t;
      bool m_joinable = false;
#  endif
    };
#endif

    class ThreadPool final : NoCopyMove {
    public:

      ThreadPool();
      ~ThreadPool();

      void changeNumberOfThreads( unsigned nthreads = 0 );
      void queue( voidfct_t );

      //To avoid deadlocks in case of recursive usage of the thread pool (jobs
      //queueing new jobs and waiting for them to finish), a thread might wish
      //to "help out" by running any pending jobs (returns an empty function if
      //there are no pending jobs). This requires very careful usage.
      voidfct_t getPendingJob();

#ifndef NCRYSTAL_DISABLE_THREADS
    private:
      std::vector<WorkerThread> m_threads;
      std::queue<voidfct_t> m_jobqueue;
      std::mutex m_mutex;
      std::condition_variable m_condvar;
      bool m_threads_should_end = true;
      void endAllThreads();
      void threadWorkFct();
      static void threadWorkFctTrampoline( void* );
#endif
    };
  }
}
#endif

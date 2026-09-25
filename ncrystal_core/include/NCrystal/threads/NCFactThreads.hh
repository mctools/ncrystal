#ifndef NCrystal_FactThreads_hh
#define NCrystal_FactThreads_hh

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

#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  // Factories of Info objects or physics processes might optionally utilise
  // multi-threading to perform some of their work. This will be disabled by
  // default, unless enabled by a call to the enable(..) function below, or the
  // environment variable NCRYSTAL_FACTORY_THREADS is set(*). An explicit
  // enable(..) call always takes precedence over the environment variable.
  // Values of 0 or 1 in either mean that threads are disabled. If NCrystal is
  // built with the NCRYSTAL_DISABLE_THREADS setting, threads can not be
  // enabled and any attempt to enable them will be silently ignored.
  //
  // Note for NCrystal plugin developers: there is no generic mechanism to signal
  // when a job has finished running in the thread pool, so callers must devise
  // their own appropriate mechanism for signaling this (most likely based on
  // captured condition or atomic variables): For this reason, it is strongly
  // recommended that any plugin wishing to take advantage of multi-threading,
  // uses the utility classes in NCProcCompBldr.hh or NCFactoryJobs.hh.
  //
  //
  // (*): The NCRYSTAL_FACTORY_THREADS environment variable is only queried
  // once, the first time the thread-pool is used or queried (e.g. when
  // materials are first created), and only if enable(..) was not already
  // called. Later changes to that variable will not have any effect.

  namespace FactoryThreadPool {

    //Enable threading during object initialisation phase. Auto detection
    //implies using a number of threads appropriate for the system (using the
    //std::thread::hardware_concurrency() C++ function).
    //
    //Assuming that user code runs in single thread (at least while initialising
    //materials), this requested value is the TOTAL number of threads utilised
    //INCLUDING that user thread. Thus, a value of 0 or 1 number will disable
    //this thread pool, while for instance calling FactoryThreadPool::enable(8)
    //will result in 7 secondary worker threads being allocated.
    NCRYSTAL_API void enable( ThreadCount = ThreadCount::auto_detect() );

    //Schedule a job to be run. If a thread-pool was not enabled, the job will
    //simply be run immediately in the current thread.
    NCRYSTAL_API void queue( voidfct_t );

    //Current total number of threads used (including the user
    //thread), so 1 means that the thread-pool is disabled:
    NCRYSTAL_API ThreadCount currentThreadCount();

    //Whether the user configured the thread-pool, by calling enable(..)
    //or via the NCRYSTAL_FACTORY_THREADS env var (with any value,
    //including 0 or 1 which explicitly disables threads):
    NCRYSTAL_API bool userConfigured();

    //Whether NCrystal was built with thread support (if not, enable(..)
    //has no effect and currentThreadCount() always returns 1):
    NCRYSTAL_API bool threadsAvailable();

    namespace detail {
      //For internal usage (e.g. by queries loading many materials): use
      //the given number of threads until the matching
      //endTemporaryThreads() call, but only if the user did not
      //configure the thread-pool (cf. userConfigured()). This does not
      //count as a user configuration. Calls can be nested or concurrent
      //(only the outermost begin/end calls change the thread-pool).
      NCRYSTAL_API void beginTemporaryThreads( ThreadCount );
      NCRYSTAL_API void endTemporaryThreads();
    }

  }
}

////////////////////////////
// Inline implementations //
////////////////////////////
namespace NCRYSTAL_NAMESPACE {
  namespace FactoryThreadPool {
    namespace detail {
      struct NCRYSTAL_API FactoryJobsHandler {
        std::function<void(voidfct_t)> jobQueueFct;
        std::function<voidfct_t()> getPendingJobFct;
      };
      NCRYSTAL_API FactoryJobsHandler getFactoryJobsHandler();
    }
  }
}
#endif

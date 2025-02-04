#ifndef NCrystal_FactThreads_hh
#define NCrystal_FactThreads_hh

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

#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  // Factories of Info objects or physics processes might optionally utilise
  // multi-threading to perform some of their work. This will be disabled by
  // default, unless enabled by a call to the enable(..) function below, or the
  // environment variable NCRYSTAL_FACTORY_THREADS is set(*). If NCrystal is built
  // with the NCRYSTAL_DISABLE_THREADS setting, threads can not be enabled and
  // any attempt to enable them will be silently ignored.
  //
  // Note for NCrystal plugin developers: there is no generic mechanism to signal
  // when a job has finished running in the thread pool, so callers must devise
  // their own appropriate mechanism for signaling this (most likely based on
  // captured condition or atomic variables): For this reason, it is strongly
  // recommended that any plugin wishing to take advantage of multi-threading,
  // uses the utility classes in NCProcCompBldr.hh or NCFactoryJobs.hh.
  //
  //
  // (*): The NCRYSTAL_FACTORY_THREADS environment variable is only queried the
  // first time one of the standard factory methods (createInfo, createScatter,
  // createAbsorption) is invoked, and later changes to that variable will not
  // have any effect.

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

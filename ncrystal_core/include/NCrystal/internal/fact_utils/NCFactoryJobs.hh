#ifndef NCrystal_FactoryJobs_hh
#define NCrystal_FactoryJobs_hh

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

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  ///////////////////////////////////////////////////////////////////////////
  //                                                                       //
  // Code wishing to potentially use MT to perform work (initialisation    //
  // only!) dispatch worker functions with the following FactoryJobs       //
  // utility class. This class is cheap to instantiate, and its usage      //
  // should in non-MT mode not incur any overhead beyond a few function    //
  // calls.                                                                //
  //                                                                       //
  // If factory threading is enabled (at the time of the FactoryJobs       //
  // constructor), the jobs might be run via NCrystal's dedicated          //
  // factory thread-pool.                                                  //
  //                                                                       //
  // Note that for scatter/absorption factories, the utility class in      //
  // NCProcCompBldr.hh might be more appropriate to use (it uses a         //
  // FactorJobs internally, but has ComponentList building capabilities).  //
  //                                                                       //
  ///////////////////////////////////////////////////////////////////////////

  class FactoryJobs final : NoCopyMove {
  public:

    FactoryJobs();
    ~FactoryJobs();

    //Either run job fct immediately, or dispatch to an MT queue, as
    //appropriate:
    void queue( voidfct_t job )
    {
      if ( isMT() )
        queueMT( std::move(job) );
      else
        job();
    }

    //Wait for all queued jobs to finish (note that in some cases of recursive
    //usage of JobGroups, this might actually entail helping to run pending jobs
    //from the global thread-pool, thus avoiding a potential deadlock).
    void waitAll()
    {
      if ( isMT() )
        waitAllMT();
    }

    //Advanced functions, for example for thread-pool integration into factory
    //infrastructure:

    bool isMT() const { return m_mt != nullptr; }
    static voidfct_t getGloballyPendingJob();

  private:
    void queueMT( voidfct_t );
    void waitAllMT();
    struct MTImpl;
    MTImpl * m_mt = nullptr;
  };

}

#endif

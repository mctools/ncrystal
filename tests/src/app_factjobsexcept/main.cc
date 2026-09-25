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


// Test exception handling in FactoryJobs (cf. GitHub issue #339):
// exceptions from jobs (in worker or calling threads) must be
// rethrown in the calling thread by waitAll(), and FactoryJobs must
// wait for all jobs when destroyed during stack unwinding. Without
// thread support, jobs run immediately in queue(..), so exceptions
// propagate from there.

#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#include "NCrystal/threads/NCFactThreads.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include <atomic>
#include <iostream>
namespace NC = NCrystal;

namespace {
  //Some busy work, so jobs are still running when others throw:
  double busyWork( unsigned n )
  {
    double s = 0.0;
    for ( auto i : NC::ncrange( n ) )
      s += std::sqrt( static_cast<double>( i ) );
    return s;
  }

  bool isMT()
  {
    NC::FactoryJobs jobs;
    return jobs.isMT();
  }

  void testThrowingJobs()
  {
    std::atomic<unsigned> nstarted( 0 ), ndone( 0 );
    bool caught = false;
    try {
      NC::FactoryJobs jobs;
      for ( auto i : NC::ncrange( 50u ) ) {
        jobs.queue( [i,&nstarted,&ndone]()
        {
          ++nstarted;
          if ( i % 10 == 3 )
            NCRYSTAL_THROW2( CalcError, "failure in job " << i );
          (void)busyWork( 20000 );
          ++ndone;
        } );
      }
      jobs.waitAll();
    } catch ( NC::Error::CalcError& e ) {
      caught = true;
      const bool ok
        = NC::StrView( e.what() ).startswith("failure in job");
      std::cout << "  caught CalcError from a job: "
                << ( ok ? "OK" : "BAD" ) << std::endl;
      nc_assert_always( ok );
    }
    nc_assert_always( caught );
    //In MT mode all jobs run, otherwise it stops at the first failure:
    if ( isMT() ) {
      nc_assert_always( nstarted.load() == 50 );
      nc_assert_always( ndone.load() == 45 );
    } else {
      nc_assert_always( nstarted.load() == 4 );
      nc_assert_always( ndone.load() == 3 );
    }
  }

  void testNonStdException()
  {
    bool caught = false;
    try {
      NC::FactoryJobs jobs;
      for ( auto i : NC::ncrange( 8u ) )
        jobs.queue( [i]() { if ( i == 5 ) throw 17; } );
      jobs.waitAll();
    } catch ( int v ) {
      caught = ( v == 17 );
    }
    nc_assert_always( caught );
    std::cout << "  caught non-std exception from a job: OK"
              << std::endl;
  }

  void testUnwindingWhileJobsRun()
  {
    std::atomic<unsigned> ndone( 0 );
    const unsigned njobs = 20;
    bool caught = false;
    try {
      NC::FactoryJobs jobs;
      for ( auto i : NC::ncrange( njobs ) ) {
        (void)i;
        jobs.queue( [&ndone]()
        {
          (void)busyWork( 200000 );
          ++ndone;
        } );
      }
      NCRYSTAL_THROW( BadInput, "exception in calling thread" );
    } catch ( NC::Error::BadInput& ) {
      caught = true;
    }
    nc_assert_always( caught );
    //The FactoryJobs destructor must have waited for all jobs:
    nc_assert_always( ndone.load() == njobs );
    std::cout << "  all jobs finished when unwinding: OK" << std::endl;
  }

  void runAll()
  {
    testThrowingJobs();
    testNonStdException();
    testUnwindingWhileJobsRun();
    //FactoryJobs still works normally afterwards:
    std::atomic<unsigned> n( 0 );
    NC::FactoryJobs jobs;
    for ( auto i : NC::ncrange( 10u ) ) {
      (void)i;
      jobs.queue( [&n]() { ++n; } );
    }
    jobs.waitAll();
    nc_assert_always( n.load() == 10 );
    std::cout << "  normal usage afterwards: OK" << std::endl;
  }
}

int main()
{
  std::cout << "Without thread-pool:" << std::endl;
  NC::FactoryThreadPool::enable( NC::ThreadCount{ 1 } );
  runAll();
  std::cout << "With thread-pool (if thread support available):"
            << std::endl;
  NC::FactoryThreadPool::enable( NC::ThreadCount{ 4 } );
  runAll();
  NC::FactoryThreadPool::enable( NC::ThreadCount{ 1 } );
  return 0;
}

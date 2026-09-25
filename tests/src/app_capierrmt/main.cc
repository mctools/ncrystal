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


// Test that the error state of the C API is kept per thread, so errors
// in one thread are never reported in another. Threads are only used if
// available, via NCrystal's factory thread pool (so
// NCRYSTAL_DISABLE_THREADS is respected).

#include "NCrystal/ncrystal.h"
#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/threads/NCFactThreads.hh"
#include <cstring>
#include <iostream>
namespace NC = NCrystal;

namespace {
  //Returns number of problems seen. Even jobs trigger errors (and check
  //that they see exactly their own error), odd jobs make successful
  //calls (and check that they never see an error):
  unsigned runJob( unsigned ijob, unsigned nrepeat )
  {
    unsigned nproblems = 0;
    const std::string fn
      = "nonexistent_" + std::to_string(ijob) + ".ncmat";
    for ( auto i : NC::ncrange( nrepeat ) ) {
      (void)i;
      if ( ijob % 2 == 0 ) {
        ncrystal_info_t info = ncrystal_create_info( fn.c_str() );
        const char * msg = ncrystal_lasterror();
        if ( !ncrystal_error() || info.internal || !msg
             || !std::strstr( msg, fn.c_str() ) )
          ++nproblems;
        ncrystal_clearerror();
        if ( ncrystal_error() )
          ++nproblems;
      } else {
        double wl = ncrystal_ekin2wl( 0.025 );
        if ( ncrystal_error() || !( wl > 1.0 ) )
          ++nproblems;
      }
    }
    return nproblems;
  }
}

int main()
{
  ncrystal_sethaltonerror( 0 );
  ncrystal_setquietonerror( 1 );
  const unsigned njobs = 16;
  const unsigned nrepeat = 200;
  NC::FactoryThreadPool::enable( NC::ThreadCount{ 8 } );
  std::vector<unsigned> nproblems( njobs, 0 );
  {
    NC::FactoryJobs jobs;
    for ( auto ijob : NC::ncrange( njobs ) )
      jobs.queue( [ijob,&nproblems]()
      {
        NC::vectAt( nproblems, ijob ) = runJob( ijob, nrepeat );
      } );
    jobs.waitAll();
  }
  NC::FactoryThreadPool::enable( NC::ThreadCount{ 0 } );
  unsigned ntot = 0;
  for ( auto e : nproblems )
    ntot += e;
  std::cout << "Ran " << njobs << " jobs with " << nrepeat
            << " C API calls each (half failing), concurrently if"
            << " possible. Problems: " << ntot << std::endl;
  nc_assert_always( ntot == 0 );
  nc_assert_always( !ncrystal_error() );
  return 0;
}

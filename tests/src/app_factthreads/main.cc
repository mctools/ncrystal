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


// Test the factory thread-pool configuration logic, in particular that
// temporary threads (for internal usage) never override a user
// configuration and leave no trace. With NCRYSTAL_DISABLE_THREADS, the
// thread count is always 1.

#include "NCrystal/threads/NCFactThreads.hh"
#include <iostream>
namespace NC = NCrystal;
namespace FTP = NC::FactoryThreadPool;

namespace {
  void check( unsigned expected_nthreads, bool expected_userconf )
  {
    //NB: Print the expected values (the same for all builds):
    const unsigned nexp
      = FTP::threadsAvailable() ? expected_nthreads : 1;
    std::cout << "  nthreads=" << expected_nthreads
              << " (1 without thread support) userConfigured="
              << expected_userconf << std::endl;
    nc_assert_always( FTP::currentThreadCount().get() == nexp );
    nc_assert_always( FTP::userConfigured() == expected_userconf );
  }
}

int main()
{
  if ( FTP::userConfigured() ) {
    //NCRYSTAL_FACTORY_THREADS set in environment, can not run test:
    std::cout << "Skipping test (user configured threads)" << std::endl;
    return 0;
  }
  std::cout << "Initial state:" << std::endl;
  check( 1, false );
  std::cout << "Nested temporary threads:" << std::endl;
  FTP::detail::beginTemporaryThreads( NC::ThreadCount{ 4 } );
  check( 4, false );
  FTP::detail::beginTemporaryThreads( NC::ThreadCount{ 8 } );
  check( 4, false );
  FTP::detail::endTemporaryThreads();
  check( 4, false );
  FTP::detail::endTemporaryThreads();
  check( 1, false );
  std::cout << "Explicitly disabled:" << std::endl;
  FTP::enable( NC::ThreadCount{ 1 } );
  check( 1, true );
  FTP::detail::beginTemporaryThreads( NC::ThreadCount{ 4 } );
  check( 1, true );
  FTP::detail::endTemporaryThreads();
  check( 1, true );
  std::cout << "Explicitly enabled, also during temporary threads:"
            << std::endl;
  FTP::detail::beginTemporaryThreads( NC::ThreadCount{ 4 } );
  FTP::enable( NC::ThreadCount{ 3 } );
  check( 3, true );
  FTP::detail::endTemporaryThreads();
  check( 3, true );
  FTP::enable( NC::ThreadCount{ 0 } );
  check( 1, true );
  return 0;
}

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
// NCRYSTAL_DISABLE_THREADS is respected). Also test that exceptions which do
// not derive from std::exception (here thrown by a custom random generator
// function) do not escape the C API, but results in the usual error state.

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

namespace {
  double badRandGen()
  {
    throw 17;
    return 0.5;
  }

  void testOutputsOnError()
  {
    //All output parameters must be set to well defined values on errors:
    ncrystal_atomdata_t ad;
    ad.internal = nullptr;//invalid handle
    const char * lbl = nullptr;
    const char * descr = nullptr;
    double mass(1.0), incxs(1.0), cohsl(1.0), absxs(1.0);
    unsigned ncomp(17), zval(18), aval(19);
    ncrystal_atomdata_getfields( ad, &lbl, &descr, &mass, &incxs, &cohsl,
                                 &absxs, &ncomp, &zval, &aval );
    nc_assert_always( ncrystal_error() );
    ncrystal_clearerror();
    nc_assert_always( lbl && descr && !lbl[0] && !descr[0] );
    nc_assert_always( mass < 0.0 && incxs < 0.0 && cohsl < 0.0
                      && absxs < 0.0 );
    nc_assert_always( ncomp == 0 && zval == 0 && aval == 0 );
    double fraction(17.0);
    ncrystal_atomdata_t sub
      = ncrystal_create_atomdata_subcomp( ad, 0, &fraction );
    nc_assert_always( ncrystal_error() );
    ncrystal_clearerror();
    nc_assert_always( sub.internal == nullptr && fraction == -1.0 );
    std::cout << "Output parameters well defined after error" << std::endl;
  }

  void testNonStdException()
  {
    ncrystal_setrandgen( badRandGen );
    nc_assert_always( !ncrystal_error() );
    ncrystal_scatter_t sc
      = ncrystal_create_scatter( "stdlib::Al_sg225.ncmat;comp=inelas" );
    nc_assert_always( !ncrystal_error() );
    double ekin_final(-1.0), mu(-1.0);
    bool escaped = false;
    try {
      ncrystal_samplescatterisotropic( sc, 0.025, &ekin_final, &mu );
    } catch ( ... ) {
      escaped = true;
    }
    std::cout << "Non-standard exception escaped C API: "
              << ( escaped ? "yes" : "no" ) << std::endl;
    nc_assert_always( !escaped );
    nc_assert_always( ncrystal_error() );
    std::cout << "Error type: " << ncrystal_lasterrortype() << std::endl;
    std::cout << "Error message: " << ncrystal_lasterror() << std::endl;
    ncrystal_clearerror();
    ncrystal_unref( &sc );
    ncrystal_setbuiltinrandgen();
    nc_assert_always( !ncrystal_error() );
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
  testNonStdException();
  testOutputsOnError();
  return 0;
}

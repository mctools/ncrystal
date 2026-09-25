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


// Test that readEntireFileToString can be used concurrently (it once
// used a shared static read buffer). Threads are only used if
// available, via NCrystal's factory thread pool (so
// NCRYSTAL_DISABLE_THREADS is respected).

#include "NCrystal/internal/utils/NCFileUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#include "NCrystal/threads/NCFactThreads.hh"
#include <cstdio>
#include <fstream>
#include <iostream>
namespace NC = NCrystal;

namespace {
  //Content spanning many 4096 byte read blocks, distinct for each file:
  std::string fileContent( unsigned ifile )
  {
    std::string s;
    for ( auto i : NC::ncrange( 20000u ) )
      s += std::to_string( ifile * 1000003u + i ) + '\n';
    return s;
  }
}

int main()
{
  const unsigned nfiles = 8;
  const unsigned nrepeat = 20;
  NC::VectS paths, contents;
  for ( auto i : NC::ncrange( nfiles ) ) {
    paths.push_back( "mtfileread_" + std::to_string(i) + ".txt" );
    contents.push_back( fileContent(i) );
    std::ofstream( paths.back(), std::ios::binary ) << contents.back();
  }

  NC::FactoryThreadPool::enable( NC::ThreadCount{ 8 } );
  std::vector<unsigned> nbad( nfiles * nrepeat, 0 );
  {
    NC::FactoryJobs jobs;
    for ( auto i : NC::ncrange( nfiles * nrepeat ) ) {
      jobs.queue( [i,&paths,&contents,&nbad]()
      {
        const unsigned ifile = i % nfiles;
        auto res
          = NC::readEntireFileToString( NC::vectAt( paths, ifile ) );
        if ( !res.has_value()
             || res.value() != NC::vectAt( contents, ifile ) )
          NC::vectAt( nbad, i ) = 1;
      } );
    }
    jobs.waitAll();
  }
  NC::FactoryThreadPool::enable( NC::ThreadCount{ 0 } );

  unsigned ntotbad = 0;
  for ( auto e : nbad )
    ntotbad += e;
  std::cout << "Read " << nfiles << " files (each "
            << contents.front().size() << " bytes) " << nrepeat
            << " times each, concurrently if possible. Bad reads: "
            << ntotbad << std::endl;
  nc_assert_always( ntotbad == 0 );
  for ( auto& p : paths )
    std::remove( p.c_str() );
  return 0;
}

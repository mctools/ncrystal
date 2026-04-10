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

#include "NCrystal/internal/minimc/NCMMC_Query.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/threads/NCFactThreads.hh"
#include <iostream>
#include <fstream>

namespace NC = NCrystal;

int main( int argc, char** argv )
{
  if ( argc!=6 && argc!=7 ) {
    std::cout
      <<"Please provide args: [cfgstr] [geomcfg] [srccfg] [enginecfg] [outfile] [factory_nthreads]"
      <<std::endl;
    std::cout
      <<"(outfile should be a json file, \"none\" or \"stdout\")."
      <<std::endl;
    std::cout
      <<"(factory_threads should be \"auto\", \"off\" or an integer value - default is \"auto\")."
      <<std::endl;
    return 1;
  }
  const std::string dest(argv[5]);

  std::string fact_threads("auto");
  if ( argc==7 )
    fact_threads = argv[6];
  if ( fact_threads == "off" ) {
    //do nothing
  } else if ( fact_threads == "auto" ) {
    NC::FactoryThreadPool::enable();
  } else {
    auto nt = NC::str2int32( fact_threads,
                             "factory threads should be auto/off/<integer>");
    nc_assert_always(nt>=0&&nt<=9999);
    NC::FactoryThreadPool
      ::enable(NC::ThreadCount(static_cast<std::uint32_t>(nt)));
  }

  std::string res;
  {
    //Invoke MiniMC:
    NC::MiniMC::Query::Query query;
    query.emplace_back("mmc");
    query.emplace_back("run");
    query.emplace_back(argv[1]);
    query.emplace_back(argv[2]);
    query.emplace_back(argv[3]);
    query.emplace_back(argv[4]);
    std::ostringstream tmp;
    NC::MiniMC::Query::JSONQuery( tmp, query );
    res = std::move(tmp).str();
  }

  if ( dest == "none " ) {
    //nothing
  } else if ( dest == "stdout" ) {
    std::cout<<res<<std::endl;
  } else {
    {
      std::ofstream fh(dest);
      fh << res;
      fh.close();
    }
    std::cout<<"Wrote: "<<dest<<std::endl;
  }
  return 0;
}

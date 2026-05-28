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

#include "NCrystal/internal/query/NCQuery.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <iostream>
#include <fstream>

namespace NC = NCrystal;

int main( int argc, char** argv )
{
  if ( argc<3 ) {
    std::cout
      <<"Please provide args: [outfile] [query]"
      <<std::endl;
    std::cout
      <<"(outfile should be a json file, \"none\" or \"stdout\")."
      <<std::endl;
    return 1;
  }
  const std::string dest(argv[1]);

  std::string res;
  {
    NC::Query query;
    for ( int i = 2; i < argc; ++ i)
      query.emplace_back(argv[i]);
    std::ostringstream tmp;
    NC::JSONQuery( tmp, query );
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

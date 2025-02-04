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

#include "NCrystal/internal/utils/NCPointwiseDist.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>

//this test creates a uniform distribution in a fancy way
int main(int , char**)
{
  const double xmax = 4.0;
  NCrystal::VectD xv = {0, 1 , 2 , 3 , xmax};
  NCrystal::VectD w2v = {4, 4 , 4, 4 , 4  };
  NCrystal::PointwiseDist dist2 (xv, w2v);
  namespace NC = NCrystal;
  for ( auto r : NC::linspace(0.0,1.0,117) ) {
    const double expected_x_at_percentile_r  = r * xmax;
    const double p = dist2.percentile(r);
    nc_assert_always( NC::floateq(expected_x_at_percentile_r, p, 1e-12, 1e-12) );
    std::cout<<"percentile( "<<NC::fmtg(r)<<" ) = "<<NC::fmtg(p)<<std::endl;
  }

  return 0;
}

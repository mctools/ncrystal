////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

#include "NCrystal/internal/NCPointwiseDist.hh"
#include <vector>

//this test creates a uniform distribution in a fancy way
int main(int , char**)
{
  double x[5]  = {0, 1 , 2 , 3 , 4};
  //  double w1[5] = {0, 1 , 2, 3 , 4 };
  double w2[5] = {4, 4 , 4, 4 , 4  };
  NCrystal::VectD xv(x,x+5);
  // NCrystal::VectD w1v(w1,w1+5);
  NCrystal::VectD w2v(w2,w2+5);

  // NCrystal::PointwiseDist dist1 (xv, w1v, 0.5);
  NCrystal::PointwiseDist dist2 (xv, w2v);

  // dist1 += dist2;
  for(unsigned i=0;i<10;i++)
  {
    double r=i*0.1;
    //the accuracy of percentile(double r) is 1e-12 (fixme: actually now it is machine precision...)
    printf("%0.8e\n", dist2.percentile(r)*0.25-r);
  }

  return 0;
}

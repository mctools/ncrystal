#ifndef NCrystal_PointwiseDist_hh
#define NCrystal_PointwiseDist_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2019 NCrystal developers                                   //
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

#include "NCrystal/NCDefs.hh"
#include <vector>
#include <stdexcept>
#include <cstdio>

namespace NCrystal {
  //Class which provides random sampling of a 1D piece-wise linear
  //distribution function. The function is defined by its non-negative
  //values on a given set of points, which must be specified in
  //ascending order. The integral_weight parameter is only used when
  //merging two distributions.
  class PointwiseDist {
  public:
    PointwiseDist(const std::vector<double> &x, const std::vector<double> &y, double integral_weight=1. );
    ~PointwiseDist();
    double percentile( double percentile_value ) const;//p in [0,1]
    double sample(RandomBase* rng) const { return percentile(rng->generate()); }
    PointwiseDist& operator+=(const PointwiseDist&);
    PointwiseDist& operator*=(double frac);
    void setIntegralWeight(double);
    void print() const;
  private:
    std::vector<double> m_cdf;
    std::vector<double> m_x;
    std::vector<double> m_y;
    double m_iweight;
  };
}

#endif

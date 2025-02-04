#ifndef NCrystal_PointwiseDist_hh
#define NCrystal_PointwiseDist_hh

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

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  // Utility class which provides integration or sampling of a 1D piece-wise
  // linear distribution function. The function is defined by its non-negative
  // values on a given set of points, which must form a proper grid of
  // increasing non-identical values.

  class PointwiseDist {
  public:
    PointwiseDist( const VectD& x, const VectD& y );
    PointwiseDist( VectD&& x, VectD&& y );

    //Percentile (argument must be in [0,1]):
    double percentile( double percentile_value ) const { return percentileWithIndex(percentile_value).first; }

    //Sample:
    double sample(RNG& rng) const { return percentileWithIndex(rng()).first; }

    const VectD& getXVals() const { return m_x; }
    const VectD& getYVals() const { return m_y; }

    //Convenience constructor (would not be needed if we had C++17's std::make_from_tuple):
    PointwiseDist(const std::pair<VectD,VectD>& xy ) : PointwiseDist(xy.first,xy.second) {}
    PointwiseDist(std::pair<VectD,VectD>&& xy ) : PointwiseDist(std::move(xy.first),std::move(xy.second)) {}

    //versions which also returns index of bin in which returned value resides
    //(i.e returns (value,idx) where value will lie in interval
    //[getXVals().at(idx),getXVals().at(idx+1)]):
    std::pair<double,unsigned> percentileWithIndex( double percentile_value ) const;
    std::pair<double,unsigned> sampleWithIndex( RNG& rng ) const { return percentileWithIndex(rng()); }

    //Sample distribution, truncated at some value (throws BadInput exception if
    //xtrunc is less than first x-value in distribution):
    double sampleBelow( RNG& rng, double xtrunc ) const;

    //Access CDF (normalised, so last value is 1.0):
    const VectD& getCDF() const { return m_cdf; }

    //Access the integral from -inf to x (again, of the normalised function, so
    //returns 1.0 if x >= last value, and 0.0 if x <= first value)):
    double commulIntegral( double x ) const;

  private:
    //todo: We have both m_cdf and m_y, although they essentially contain the
    //same info. Could we implement more light-weight version? Could we
    //implement as a non-owning view, i.e. which keeps m_x in span (but likely
    //needs to be possible to be owning still). Or using shared ptrs?
    VectD m_cdf;
    VectD m_x;
    VectD m_y;
  };
}

#endif

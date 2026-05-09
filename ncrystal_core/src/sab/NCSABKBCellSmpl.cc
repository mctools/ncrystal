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

#include "NCrystal/internal/sab/NCSABKBCellSmpl.hh"
#include "NCrystal/internal/utils/NCString.hh"
namespace NC = NCrystal;
namespace NCSU = NCrystal::SABUtils;

NCSU::ParabolicBandBoxSampler::ParabolicBandBoxSampler( double x0, double y0,
                                                        double x1, double y1 )
{
  //m_data encoding:
  // (-1.0,-1.0,0,0)  [if can not sample]
  // (x0,y0,x1,y1)    [if box fully inside band]
  // (x0,y0,-x1,-y1)  [reject_parallelogram]
  // (x0,y0,-x1,y1)   [reject_box]

  //First we want to trim the box as tightly as possible to the parabolic band,
  //while still keeping the box axis-aligned.
  //
  //We have y^\pm=(1\pm\sqrt(x))^2 and only need to consider the quadrant x>=0
  //and y>=0. The trimming is more conveniently done in the variables
  //(u,v)=(sqrt(x),sqrt(y)). Here the y+ edge simply corresponds to v=1+u, and
  //y- corresponds to v=|1-u|, which is 1-u when u<=1 and u-1 when u>=1.

  nc_assert( std::isfinite(x0) );
  nc_assert( std::isfinite(y0) );
  nc_assert( std::isfinite(x1) );
  nc_assert( std::isfinite(y1) );
  nc_assert( x1 >= x0 );
  nc_assert( y1 >= y0 );

  double u0 = std::sqrt(x0);
  double u1 = std::sqrt(x1);
  double v0 = std::sqrt(y0);
  double v1 = std::sqrt(y1);

  v1 = ncmin( v1, 1.0 + u1 );//trim down to v=1+u line
  u0 = ncmax( u0, v0 - 1.0 );//trim rightwards to v=1+u line
  v0 = ncmax( v0, u0 - 1.0 );//trim upwards to v=u-1 line
  u1 = ncmin( u1, v1 + 1.0 );//trim leftwards to v=u-1 line
  v0 = ncmax( v0, 1.0 - u1 );//trim upwards to v=1-u line
  u0 = ncmax( u0, 1.0 - v1 );//trim rightwards to v=1-u line

  if ( !( (v1>v0) && (u1>u0) ) ) {
    //not a valid box -> it must have started empty or must not have overlapped
    //with the parabolic band.
    m_data[0] = -1.0;
    m_data[1] = -1.0;
    m_data[2] = 0.0;
    m_data[3] = 0.0;
  } else {
    //Let us check if the box is now fully inside, still in (u,v) space.
    const bool fully_inside = ( v1 <= 1.0+u0 && v0 >= 1.0-u0 && v0 >= u1-1.0 );
    const double res_x0(u0*u0), res_y0(v0*v0), res_x1(u1*u1), res_y1(v1*v1);
    m_data[0] = res_x0;
    m_data[1] = res_y0;
    m_data[2] = res_x1;
    m_data[3] = res_y1;
    if ( !fully_inside ) {
      m_data[2] = -m_data[2];
      //Detect "parallelogram" scenario where an axis-aligned box overlay has
      //low acceptance rate. The required parallelogram will end up with an area
      //of 4sqrt(x1)*(x1-x0), whereas the simple box has an area of
      //(y1-y0)*(x1-x0). So we simply pick the parallelogram method whenever it
      //would have a smaller area than the box. However, as the parallelogram
      //method also requires a few extra calculations (most notable a sqrt
      //evaluation), we artificially increase its area for the purposes of this
      //decision.
      constexpr double parallelogram_penalty = 1.2;//tunable! (fixme)
      constexpr double kkk = ncsquare(4.0*parallelogram_penalty);
      if ( kkk * res_x1 < ncsquare(res_y1-res_y0) )
        m_data[3] = -m_data[3];
    }
  }
}

void NCSU::ParabolicBandBoxSampler::toJSON( std::ostream& os ) const
{
  const char * samplemode = ( canSample()
                              ? ( fullyInside()
                                  ? "fullbox"
                                  : ( m_data[3]<0.0
                                      ? "parallelogram_reject"
                                      : "box_reject" )
                                  )
                              : "unable" );
  os << "{\"samplemode\":";
  streamJSON(os,samplemode);
  os <<",\"overlay_region_curve\":";
  std::vector<PairDD> v;
  if (canSample()) {
    const double x0 = m_data[0];
    const double y0 = m_data[1];
    const double x1 = ncabs(-m_data[2]);
    const double y1 = ncabs(-m_data[3]);
    if ( m_data[3]<0.0 ) {
      //parallelogram
      const double sqrtx1 = std::sqrt( x1 );
      const double oneminustwosqrtx1 = 1.0 - 2.0*sqrtx1;
      const double foursqrtx1 = 4.0*sqrtx1;
      v.emplace_back(x0,oneminustwosqrtx1+x0);
      v.emplace_back(x1,oneminustwosqrtx1+x1);
      v.emplace_back(x1,oneminustwosqrtx1+foursqrtx1+x1);
      v.emplace_back(x0,oneminustwosqrtx1+foursqrtx1+x0);
    } else {
      //box
      v.emplace_back(x0,y0);
      v.emplace_back(x1,y0);
      v.emplace_back(x1,y1);
      v.emplace_back(x0,y1);
    }
  }
  if (!v.empty())
    v.emplace_back(v.front());//close loop
  streamJSON(os,v);
  os <<'}';
}

NC::PairDD
NCSU::ParabolicBandBoxSampler::sampleCrossingParallelogram( RNG& rng ) const
{
  //When x1 >> x0, usage of an axis-aligned overlay box has very poor acceptance
  //rate. So in this case, we instead overlay with a parallelogram which covers
  //the parabolic band and is still easy to sample uniformly.
  nc_assert(canSample());
  const double x0 =  m_data[0];
  const double y0 =  m_data[1];
  const double x1 = -m_data[2];
  const double y1 =  -m_data[3];
  nc_assert( x1 > 0.0 );
  nc_assert( y1 > 0.0 );
  const double sqrtx1 = std::sqrt( x1 );
  const double oneminustwosqrtx1 = 1.0 - 2.0*sqrtx1;
  const double foursqrtx1 = 4.0*sqrtx1;
  while( true ) {
    const double x = ncmin(x1,x0 + rng.generate()*(x1-x0));
    const double y = oneminustwosqrtx1 + rng.generate() * foursqrtx1 + x;
    if ( valueInInterval(y0,y1,y) && ncsquare(y-(1.0+x)) < 4.0 * x )
      return { x, y };
  }
}

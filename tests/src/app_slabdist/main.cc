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

#include "NCrystal/internal/minimc/NCMMC_Utils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;
namespace NCMMCU = NCrystal::MiniMC::Utils;

void test( double x, double ux, double slab_halfthickness = 1.0 )
{
  //Using very robust and simple (can be expensive) math to verify the output of
  //the more optimised algs.
  nc_assert_always(slab_halfthickness>0.0);
  nc_assert_always( ux >= -1.0 && ux <= 1.0 );
  const double dx = slab_halfthickness;
  const bool is_on_edge = ( std::fabs(x) == dx );
  const bool is_completely_inside = ( std::fabs(x) < dx );
  const bool is_completely_outside = ( std::fabs(x) > dx );
  std::cout<<"\nTesting x="<<NC::fmt(x)<<", ux="
           <<NC::fmt(ux)<<" dx="<<NC::fmt(dx)<<std::endl;

  double dist_outside;
  NCMMCU::distToSlabEntry( &x, &ux, &dist_outside, 1, dx );
  std::cout<<"  -> distFromOutside = "<<NC::fmt(dist_outside)<<std::endl;

  NC::Optional<double> dist_inside;
  if ( !is_completely_outside ) {
    double tmp;
    NCMMCU::distToSlabExit( &x, &ux, &tmp, 1, dx );
    dist_inside = tmp;
    std::cout<<"  -> distFromInside = "<<NC::fmt(dist_inside.value())
             <<std::endl;
  }

  //Carefully figure out what the correct values should be:
  double ref_dist_outside;
  double ref_dist_inside;
  //For simplicity we mirror around x=0 if x<0.
  if ( x<0.0 ) {
    ux = -ux;
    x = -x;
  }
  if ( is_on_edge ) {
    //ON EDGE
    nc_assert_always(x>=0.0);
    if ( ux == 0.0 ) {
      //already skirting along the edge, should
      //attenuate/scatter.
      ref_dist_outside = 0.0;
      ref_dist_inside = NC::kInfinity;
    } else if ( ux > 0.0 ) {
      //moving immediately out of the volume
      ref_dist_outside = -1.0;
      ref_dist_inside = 0.0;
    } else {
      //heading inside, moving to the other edge of the volume.
      nc_assert_always( ux < 0.0 );
      ref_dist_inside = 2*dx/std::fabs(ux);
      ref_dist_outside = 0.0;
    }
  } else if ( is_completely_outside ) {
    //OUTSIDE
    nc_assert_always(!dist_inside.has_value());
    ref_dist_inside = 0.1234567;//silence compiler warning
    nc_assert_always( x>=0.0 );
    if ( ux >= 0.0 ) {
      ref_dist_outside = -1.0;//missing
    } else {
      //moving to surface at a speed of 1/ux.
      nc_assert_always( ux < 0.0 );
      ref_dist_outside = (x-dx)/std::fabs(ux);
    }
  } else {
    //INSIDE
    nc_assert_always( x>=0.0 );
    nc_assert_always( is_completely_inside );
    nc_assert_always( dist_inside.has_value() );
    ref_dist_outside = 0.0;//already inside
    if ( ux == 0.0 ) {
      ref_dist_inside = NC::kInfinity;//infinite path inside
    } else if ( ux > 0.0 ) {
      nc_assert(x<dx);
      ref_dist_inside = (dx-x)/ux;
    } else {
      nc_assert_always( ux < 0.0 );
      ref_dist_inside = (x+dx)/(-ux);
    }
  }

  std::cout<<"  -> distFromOutside should be = "<<NC::fmt(ref_dist_outside)
           <<std::endl;
  if ( !is_completely_outside ) {
    std::cout<<"  -> distFromInside should be = "<<NC::fmt(ref_dist_inside)
             <<std::endl;
  }

  if ( !NC::floateq( dist_outside, ref_dist_outside ) )
    nc_assert_always(false);
  if ( !is_completely_outside ) {
    if ( !NC::floateq( dist_inside.value(), ref_dist_inside ) )
      nc_assert_always(false);
  }


}

int main(int,char**) {

  const std::vector<double> xvals =
    { -1e99, -10.0, -1.0+1e-15, -1.0, -1.0+1e-15, -0.5, -1e-15, 0.0, 1e-15,
      1.0-1e-15, 0.5, 1.0, 1.0+1e-15, 10.0, 1e99 };//fixme: if we add inf here
                                                   //we trigger asserts. We
                                                   //should test infinity slab
                                                   //gun with dir (1,0,0)
  const std::vector<double> uxvals =
    { -1.0, -1.0+1e-15, -0.5, -1e-15, 0.0, 1e-15, 0.5, 1.0-1e-15, 1.0 };
  for ( auto& slab_halfthickness : xvals ) {
    if ( slab_halfthickness <= 0.0 )
      continue;
    for ( auto& x : xvals ) {
      for ( auto& ux : uxvals ) {
        test( x, ux, slab_halfthickness );
      }
    }
  }
  std::cout<<"All tests passed!"<<std::endl;
  return 0;
}

#ifndef NCrystal_ExtnUtils_hh
#define NCrystal_ExtnUtils_hh

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

  // Various utilities for extinction models. Several functions concerns the
  // Sabine model, described in:

  // Sabine, T. M. "The flow of radiation in a real crystal." (2006): 609-616.
  // doi: 10.1107/97809553602060000603

  namespace Extn {

    double calcSabineEl( double x, double y );//Sabine 6.4.5.3-4
    double calcSabineEb( double x, double y );//Sabine 6.4.5.5

    //Same but y=0 or Eb with pre-calculated A+B factors:
    double calcSabineEl_y0( double x );
    double calcSabineEb_y0( double x );
    double calcSabineA( double y );//Sabine 6.4.5.6
    double calcSabineB( double y );//Sabine 6.4.5.7
    double calcSabineEb_cachedAB( double x, double A, double B );
  }

}

////////////////////////////
// Inline implementations //
////////////////////////////

inline double NCrystal::Extn::calcSabineEl( double x, double y )
{
  nc_assert( y>=0.0 );
  nc_assert( std::isfinite(y) );
  return std::exp(-y)*calcSabineEl_y0(x);
}

inline double NCrystal::Extn::calcSabineEb_y0( double x )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );
  return 1.0 / std::sqrt( 1.0 + x );
}

inline double NCrystal::Extn::calcSabineEb_cachedAB( double x,
                                                     double A,
                                                     double B )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );
  return A / std::sqrt( 1.0 + B * x );
}

#endif

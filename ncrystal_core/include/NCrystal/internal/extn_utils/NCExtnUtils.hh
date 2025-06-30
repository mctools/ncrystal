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

  //////////////////////////////////////////////////////////////////////////////
  //
  // Various utilities for extinction models. Several functions concerns the
  // Sabine model, described in:
  //
  // Sabine, T. M. "The flow of radiation in a real crystal." (2006): 609-616.
  // doi: 10.1107/97809553602060000603
  //
  // Note that the recipe for the Eb factor provided in Sabine's paper in
  // eqs. 6.4.5.3-4 is not used directly in NCrystal by default, since the
  // expansions were not deemed precise enough. Instead an improved formula is
  // provided, as described in the following, and the El value given by the
  // original recipe can instead be accessed via the calcSabineElOriginal
  // functions if needed.
  //
  // The problem was that the expansions for El presented by Sabine in
  // eqs. 6.4.5.3-4 have very few terms, reaching an imprecision above 1% near
  // x=1, where a visible kink of ~2% is introduced in the curves due to this
  // imprecision. Such a kink unfortunately could give the impression of a new
  // unphysical reflection in a Bragg edge plot.
  //
  // Therefore for NCrystal T. Kittelmann redid the calculations leading to
  // those equations and found that the exact closed-form expressions for El
  // should be (for any x>=0):
  //
  // El = ( I(0,x) + I(1,x) ) * exp( -x ) * exp( -y )
  //
  // Where I(nu,x) is the modified Bessel function of the first kind.
  //
  // Based on that, we implemented a higher-precision improved evaluation which
  // provides >6 significant digits of precision everywhere. This was achieved
  // by:
  //
  //   * Adding more terms to the expansions shown in Sabine's paper used for
  //     low and large x.
  //   * Adding another Taylor expansion around x=4, to be used for
  //     intermediate x-values in [1.563,6.396].
  //
  // For the exact coefficients and terms, please refer to the implementation.
  //
  //////////////////////////////////////////////////////////////////////////////

  namespace Extn {

    double calcSabineEb( double x, double y );//Sabine 6.4.5.5
    double calcSabineEl( double x, double y );//Improved version of Sabine
                                              //eqs. 6.4.5.3-4 as discussed
                                              //above. See "Original" functions
                                              //below for reference.

    //Same but y=0 or Eb with pre-calculated A(y) and B(y) factors:
    double calcSabineEb_y0( double x );
    double calcSabineEl_y0( double x );
    double calcSabineA( double y );//Sabine 6.4.5.6
    double calcSabineB( double y );//Sabine 6.4.5.7
    double calcSabineEb_cachedAB( double x, double A, double B );

    //Original recipe for El as it appeared in Sabine's paper:
    double calcSabineElOriginal( double x, double y );
    double calcSabineElOriginal_y0( double x );
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

inline double NCrystal::Extn::calcSabineElOriginal( double x, double y )
{
  nc_assert( y>=0.0 );
  nc_assert( std::isfinite(y) );
  return std::exp(-y)*calcSabineElOriginal_y0(x);
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

inline double NCrystal::Extn::calcSabineEb( double x, double y )
{
  nc_assert( y>=0.0 );
  nc_assert( std::isfinite(y) );
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );
  return calcSabineA( y ) / std::sqrt( 1.0 + calcSabineB( y ) * x );
}



#endif

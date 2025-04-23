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

#include "NCrystal/internal/extinction/NCExtinctionUtils.hh"

namespace NC = NCrystal;
namespace NCE = NCrystal::Extinction;

double NCE::calcSabineA( double y )
{
  //Evaluate exp(-y)*sinh(y)/y for y >= 0
  nc_assert( y>=0.0 );
  nc_assert( std::isfinite(y) );

  if ( y < 0.2 ) {
    //Taylor expand

    // The function exp(-y)*sinh(y)/y = (1-exp(-2y))/(2y) must be evaluated with
    // a Taylor expansion near y=0.
    // For the record we simply got the Taylor coefficients with sagemath:
    // > sage: f=(1-exp(-2*y))/(2*y)
    // First investigating number of orders needed for u=0.2 with command:
    // > sage: u=0.2;(f-f.taylor(y,0,14))(y=u).n()
    // Then generate the coefficients with:
    // > sage: print( '\n'.join(('constexpr double c%i = %s;'%(c[1],str(c[0]))).replace('/','./') for c in (f.taylor(y,0,14)).coefficients()))

    constexpr double c0 = 1;
    constexpr double c1 = -1;
    constexpr double c2 = 2./3;
    constexpr double c3 = -1./3;
    constexpr double c4 = 2./15;
    constexpr double c5 = -2./45;
    constexpr double c6 = 4./315;
    constexpr double c7 = -1./315;
    constexpr double c8 = 2./2835;
    constexpr double c9 = -2./14175;
    constexpr double c10 = 4./155925;
    constexpr double c11 = -2./467775;
    constexpr double c12 = 4./6081075;
    constexpr double c13 = -4./42567525;
    constexpr double c14 = 8./638512875;
    return c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*(c8+y*(c9+y*(c10+y*(c11+y*(c12+y*(c13+y*c14)))))))))))));
  } else {
    double m2y = -2.0*y;
    return std::expm1( m2y ) / m2y;
  }

}

double NCE::calcSabineB( double y )
{
  //Evaluate 1/y - exp(-y)/sinh(y) for y>=0

  nc_assert( y>=0.0 );
  nc_assert( std::isfinite(y) );

  if ( y < 0.3 ) {
    //Taylor expand

    // The function exp(-y)*sinh(y)/y = (1-exp(-2y))/(2y) must be evaluated with
    // a Taylor expansion near y=0.
    // For the record we simply got the Taylor coefficients with sagemath:
    // > sage: f=1/y-exp(-y)/sinh(y)
    // First investigating number of orders needed for u=0.3 with command:
    // > sage: u=0.3;(f-f.taylor(y,0,15))(y=u).n()
    // Then generate the coefficients with:
    // > sage: print( '\n'.join(('constexpr double c%i = %s;'%(c[1],str(c[0]))).replace('/','./') for c in (f.taylor(y,0,15)).coefficients()))
    constexpr double c0 = 1;
    constexpr double c1 = -1./3;
    constexpr double c3 = 1./45;
    constexpr double c5 = -2./945;
    constexpr double c7 = 1./4725;
    constexpr double c9 = -2./93555;
    constexpr double c11 = 1382./638512875;
    constexpr double c13 = -4./18243225;
    constexpr double c15 = 3617./162820783125;
    double s = y * y;
    return c0 + y*(c1+s*(c3+s*(c5+s*(c7+s*(c9+s*(c11+s*(c13+s*c15)))))));
  } else {
    return 1.0/y + 2.0/( 1.0 - std::exp( 2.0*y ) );
  }

}

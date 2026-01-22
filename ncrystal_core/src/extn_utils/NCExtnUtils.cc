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

#include "NCrystal/internal/extn_utils/NCExtnUtils.hh"

namespace NC = NCrystal;
namespace NCE = NCrystal::Extn;

double NCE::calcSabineA( double y )
{
  //Evaluate exp(-y)*sinh(y)/y for y >= 0
  nc_assert( y>=0.0 );
  nc_assert( std::isfinite(y) );

  if ( y < 0.2 ) {
    // The function exp(-y)*sinh(y)/y = (1-exp(-2y))/(2y) must be evaluated with
    // a Taylor expansion near y=0.
    // For the record we simply got the Taylor coefficients with sagemath:
    // > sage: f=(1-exp(-2*y))/(2*y)
    // First investigating number of orders needed for u=0.2 with command:
    // > sage: u=0.2;(f-f.taylor(y,0,14))(y=u).n()
    // Then generate the coefficients with:
    // > sage: print( '\n'.join(('constexpr double c%i = %s;'%(c[1],str(c[0])))
    //          .replace('/','./') for c in (f.taylor(y,0,14)).coefficients()))

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
    return c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*(c7+y*(c8
           +y*(c9+y*(c10+y*(c11+y*(c12+y*(c13+y*c14)))))))))))));
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

    // The function 1/y-exp(-y)*sinh(y) must be evaluated with
    // a Taylor expansion near y=0.
    // For the record we simply got the Taylor coefficients with sagemath:
    // > sage: f=1/y-exp(-y)/sinh(y)
    // First investigating number of orders needed for u=0.3 with command:
    // > sage: u=0.3;(f-f.taylor(y,0,15))(y=u).n()
    // Then generate the coefficients with:
    // > sage: print( '\n'.join(('constexpr double c%i = %s;'%(c[1],str(c[0])))
    //         .replace('/','./') for c in (f.taylor(y,0,15)).coefficients()))
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

double NCE::calcSabineElOriginal_y0( double x )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );

  //Note that the following is not our own Taylor expansions, but the formula as
  //presented in Sabine 6.4.5.3-5.
  if ( x <= 1.0 ) {
    constexpr double c1 = -1./2.;
    constexpr double c2 = 1./4.;
    constexpr double c3 = -5./48.;
    constexpr double c4 = 7./192.;
    return ( 1.0 + x*(c1+x*(c2+x*(c3+x*c4))) );
  } else {
    constexpr double c1 = -1./8.;
    constexpr double c2 = -3./128.;
    constexpr double c3 = -15./1024.;
    constexpr double k = 2.0 * kInvSqrt2Pi; // = sqrt(2/pi)
    const double invx = 1.0 / x;
    return ( k*std::sqrt(invx) ) * ( 1.0+invx*(c1+invx*(c2+invx*c3)) );
  }
}

double NCE::calcSabineEl_y0( double x )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );

  if ( x <= 1.563 ) {
    //Taylor expansion around x=0

    //For the record, this was generated in sagemath with:
    //sage: f
    //(bessel_I(1, x) + bessel_I(0, x))*e^(-x)
    //sage: print( '\n'.join(('constexpr double c%i = %.19g;'%(c[1],c[0]
    //      .n(digits=20))) for c in (f.taylor(x,0,14).coefficients())))
    constexpr double c0 = 1;
    constexpr double c1 = -0.5;
    constexpr double c2 = 0.25;
    constexpr double c3 = -0.1041666666666666713;
    constexpr double c4 = 0.03645833333333333565;
    constexpr double c5 = -0.01093749999999999931;
    constexpr double c6 = 0.002864583333333333131;
    constexpr double c7 = -0.0006649925595238095465;
    constexpr double c8 = 0.0001385401165674603132;
    constexpr double c9 = -2.616868868496472552e-05;
    constexpr double c10 = 4.520046227402997952e-06;
    constexpr double c11 = -7.190982634504769757e-07;
    constexpr double c12 = 1.060208978164164819e-07;
    constexpr double c13 = -1.456331013961764904e-08;
    constexpr double c14 = 1.872425589379411794e-09;
    return (c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*(c8
           +x*(c9+x*(c10+x*(c11+x*(c12+x*(c13+x*c14))))))))))))));
  }

  if ( x <= 6.396 ) {
    //Taylor expansion around x=4

    //For the record, this was generated in sagemath with:
    //sage: f
    //(bessel_I(1, x) + bessel_I(0, x))*e^(-x)
    //sage: print( '\n'.join(('constexpr double c%i = %.19g;'%(c[1],c[0]
    //      .n(digits=20))) for c in (f.taylor(x,4,14)(x=u+4).coefficients())))
    constexpr double c0 = 0.385752760726422006;
    constexpr double c1 = -0.04468770987560883451;
    constexpr double c2 = 0.007640542253708286317;
    constexpr double c3 = -0.001417866664781492056;
    constexpr double c4 = 0.0002675101804915058475;
    constexpr double c5 = -4.981828415109825685e-05;
    constexpr double c6 = 8.997402217944376013e-06;
    constexpr double c7 = -1.558125440195910587e-06;
    constexpr double c8 = 2.568920556636835166e-07;
    constexpr double c9 = -4.015777071986909606e-08;
    constexpr double c10 = 5.93998308426714981e-09;
    constexpr double c11 = -8.309071416664671749e-10;
    constexpr double c12 = 1.09961981812523942e-10;
    constexpr double c13 = -1.378133839491552616e-11;
    constexpr double c14 = 1.637930268726293705e-12;
    const double u = x - 4.0;
    return (c0+u*(c1+u*(c2+u*(c3+u*(c4+u*(c5+u*(c6+u*(c7+u*(c8
           +u*(c9+u*(c10+u*(c11+u*(c12+u*(c13+u*c14))))))))))))));
  }

  {
    //Large x approximation.
    //
    //Based on asymptotic behaviour functions for the Bessel functions, we find
    //the asymptotic behaviour like in Sabine's paper, but with more terms:
    //
    constexpr double c1 = 1./8;
    constexpr double c2 = 3./128;
    constexpr double c3 = 15./1024;
    constexpr double c4 = 525./32768;
    constexpr double c5 = 6615./262144;
    constexpr double c6 = 218295./4194304;
    constexpr double c7 = 4459455./33554432;
    constexpr double c8 = 869593725./2147483648;
    constexpr double k = 2.0 * kInvSqrt2Pi; // = sqrt(2/pi)
    const double u = 1.0 / x;
    return ( k*std::sqrt(u)
             * (1.0-u*(c1+u*(c2+u*(c3+u*(c4+u*(c5+u*(c6+u*(c7+u*c8)))))))) );
  }
}

double NCE::calcSabineEb_ScndTriang_y0( double x )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );

  //Limit of Sabine 6.4.9.5 as y->0 and thus A->1 and B->1 is:
  //
  // Eb = 2*(x-ln(1+x))/x^2
  //
  if ( x < 0.1 ) {
    // The function (x-ln(1+x))/(2*x^2) must be evaluated with
    // a Taylor expansion near x=0.
    // For the record we simply got the Taylor coefficients with sagemath:
    // > sage: f=2*(x-ln(1+x))/x^2
    // First investigating number of orders needed for u=0.1 with command:
    // > sage: u=0.1;(f-f.taylor(x,0,13))(x=u).n()
    // Then generate the coefficients with:
    // > sage: print( '\n'.join(('constexpr double c%i = %s;'%(c[1],str(c[0])))
    //          .replace('/','./') for c in (f.taylor(x,0,13)).coefficients()))
    constexpr double c0 = 1;
    constexpr double c1 = -2./3;
    constexpr double c2 = 1./2;
    constexpr double c3 = -2./5;
    constexpr double c4 = 1./3;
    constexpr double c5 = -2./7;
    constexpr double c6 = 1./4;
    constexpr double c7 = -2./9;
    constexpr double c8 = 1./5;
    constexpr double c9 = -2./11;
    constexpr double c10 = 1./6;
    constexpr double c11 = -2./13;
    constexpr double c12 = 1./7;
    constexpr double c13 = -2./15;
    return (c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*(c8
           +x*(c9+x*(c10+x*(c11+x*(c12+x*c13)))))))))))));
  } else {
    return 2.0*( x - std::log1p(x) ) / (x*x);
  }

}
double NCE::calcSabineEl_ScndRect_y0( double x )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );

  // Evaluating Sabine 6.4.9.2 [y=0], in other words:
  //
  // El = ( 1-exp(-2x) ) / (2x)
  //
  if ( x < 0.2 ) {
    // The function must be evaluated with
    // a Taylor expansion near x=0.
    // For the record we simply got the Taylor coefficients with sagemath:
    // > sage: f = ( 1 - exp(-2*x) ) / (2*x)
    // First investigating number of orders needed for u=0.1 with command:
    // > sage: u=0.2;(f-f.taylor(x,0,13))(x=u).n()
    // Then generate the coefficients with:
    // > sage: print( '\n'.join(('constexpr double c%i = %s;'%(c[1],str(c[0]))).replace('/','./') for c in (f.taylor(x,0,13)).coefficients()))
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
    return (c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*(c8
            +x*(c9+x*(c10+x*(c11+x*(c12+x*c13)))))))))))));
  } else {
    return ( -0.5 / x ) *  std::expm1( -2.0 * x );
  }
}

double NCE::calcSabineEl_ScndTriang_y0( double x )
{
  nc_assert( x>=0.0 );
  nc_assert( std::isfinite(x) );

  // Evaluating Sabine 6.4.9.4 [y=0], in other words:
  //
  // El = ( 1 - (1-exp(-2x))/2x ) / x
  //
  if ( x < 0.2 ) {
    // The function must be evaluated with
    // a Taylor expansion near x=0.
    // For the record we simply got the Taylor coefficients with sagemath:
    // > sage: f=( 1 - (1-exp(-2*x))/(2*x) ) / x
    // First investigating number of orders needed for u=0.1 with command:
    // > sage: u=0.2;(f-f.taylor(x,0,12))(x=u).n()
    // Then generate the coefficients with:
    // > sage: print( '\n'.join(('constexpr double c%i = %s;'%(c[1],str(c[0])))
    //          .replace('/','./') for c in (f.taylor(x,0,12)).coefficients()))
    constexpr double c0 = 1;
    constexpr double c1 = -2./3;
    constexpr double c2 = 1./3;
    constexpr double c3 = -2./15;
    constexpr double c4 = 2./45;
    constexpr double c5 = -4./315;
    constexpr double c6 = 1./315;
    constexpr double c7 = -2./2835;
    constexpr double c8 = 2./14175;
    constexpr double c9 = -4./155925;
    constexpr double c10 = 2./467775;
    constexpr double c11 = -4./6081075;
    constexpr double c12 = 4./42567525;
    return (c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*(c8
           +x*(c9+x*(c10+x*(c11+x*c12))))))))))));
  } else {
    const double invx = 1.0 / x;
    return invx * ( 1.0 + 0.5 * invx * std::expm1( -2.0 * x ) );
  }
}

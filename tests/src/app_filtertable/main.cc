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

//Tests of the filter table creation with synthetic cross section curves: the
//tolerance must be met for complicated shapes, and problems which could lead
//to wrong tables must result in exceptions.

#include "NCrystal/internal/filter/NCFilterTable.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>

namespace NC = NCrystal;
namespace NCF = NCrystal::Filter;

namespace {

  using Fct = std::function<double(double)>;

  //Worst relative error of the table vs. the function, at log-spaced points
  //and at points close to the given discontinuities:
  double worstRelErr( const NCF::Table& t, const Fct& f, const NC::VectD& discont,
                      double wlmax )
  {
    NC::VectD wls = NC::geomspace( 1e-7, wlmax, 200000 );
    for ( auto e : discont )
      for ( auto r : { 1-1e-8, 1+1e-8, 1-1e-5, 1+1e-5 } )
        wls.push_back( e * r );
    //The tolerance is relative to max(xs,1e-12):
    const double abs_floor = 1e-12;
    double worst = 0.0;
    for ( auto wl : wls ) {
      const double exact = f( wl );
      const double tab = NCF::evalTable( t.wl.data(), t.xs.data(), t.wl.size(), wl );
      const double yref = std::max( exact, abs_floor );
      const double err = ( yref > 0.0
                           ? NC::ncabs( tab - exact ) / yref
                           : NC::ncabs( tab - exact ) );
      worst = std::max( worst, err );
    }
    return worst;
  }

  void testOK( const char* descr, const Fct& f, NC::VectD discont = {},
               double tol = 1e-3, double wlmax = 500.0 )
  {
    NCF::TableParams params;
    params.tol = tol;
    params.wlmax = wlmax;
    auto t = NCF::createTable( f, discont, params );
    const double worst = worstRelErr( t, f, discont, wlmax );
    nc_assert_always( worst <= tol );
    std::cout << descr << ": npts=" << t.wl.size()
              << " xs(0)=" << std::round( t.xs.front() * 1e6 ) / 1e6
              << " worst/tol<1: OK" << std::endl;
  }

  template <class TError>
  void testError( const char* descr, const Fct& f, NC::VectD discont = {},
                  NCF::TableParams params = {} )
  {
    try {
      NCF::createTable( f, discont, params );
    } catch ( TError& e ) {
      std::cout << descr << ": got expected exception: " << e.what() << std::endl;
      return;
    }
    std::cout << descr << ": did not get the expected exception!" << std::endl;
    nc_assert_always( false );
  }

  double step( double wl, double e, double below, double above )
  {
    return wl < e ? below : above;
  }

}

int main()
{
  //Shapes which must be handled:
  testOK( "constant", []( double ) { return 5.0; } );
  testOK( "zero", []( double ) { return 0.0; } );
  testOK( "a+b*wl", []( double wl ) { return 3.0 + 0.5 * wl; } );
  testOK( "a/wl^0.1 at long wavelengths", []( double wl ) { return 1.0 + wl / std::pow( 1.0 + wl, 0.1 ); } );
  testOK( "vanishing at wl=0", []( double wl ) { return wl * wl / ( 1.0 + wl ); } );
  testOK( "vanishing at wl=0 with large curvature", []( double wl ) { return 10.0 * wl * wl / ( 1.0 + wl * wl ); } );
  testOK( "oscillating", []( double wl ) { return 2.0 + std::sin( 20.0 * wl ); }, {}, 1e-3, 20.0 );
  //Resonance-like peaks (Lorentzians in energy), as long as they are wider
  //than the dense sampling:
  testOK( "peaks", []( double wl )
  {
    const double e = NC::wl2ekin( wl );
    double v = 1.0 + 0.3 * wl;
    for ( double e0 : { 0.5, 3.0, 20.0 } ) {
      const double g = 0.01 * e0;//width relative to the energy: 1%
      v += 50.0 * g * g / ( ( e - e0 ) * ( e - e0 ) + g * g );
    }
    return v;
  } );
  //Declared discontinuities (up and down), also with kinks:
  testOK( "steps", []( double wl )
  {
    return ( 1.0 + 0.1 * wl + step( wl, 2.0, 0.0, 3.0 )
             + step( wl, 4.0, 5.0, 0.0 ) + step( wl, 7.0, 0.0, 0.2 * ( wl - 7.0 ) ) );
  }, { 2.0, 4.0, 7.0 } );
  testOK( "drop to zero", []( double wl ) { return step( wl, 5.0, 2.0 + wl, 0.0 ); }, { 5.0 } );

  //Problems which must result in exceptions:
  testError<NC::Error::CalcError>( "unexpected step", []( double wl )
  {
    return 1.0 + step( wl, 2.0, 0.0, 3.0 );
  } );
  testError<NC::Error::CalcError>( "no limit for wl->0", []( double wl ) { return 1.0 / wl; } );
  testError<NC::Error::CalcError>( "negative", []( double wl ) { return 1.0 - wl; } );
  testError<NC::Error::CalcError>( "NaN", []( double wl ) { return wl > 3.0 ? std::numeric_limits<double>::quiet_NaN() : 1.0; } );
  testError<NC::Error::CalcError>( "infinite", []( double wl ) { return wl > 3.0 ? NC::kInfinity : 1.0; } );
  {
    NCF::TableParams p;
    p.tol = 0.0;
    testError<NC::Error::BadInput>( "tol=0", []( double ) { return 1.0; }, {}, p );
    p = NCF::TableParams();
    p.wlmax = 1e-9;
    testError<NC::Error::BadInput>( "tiny wlmax", []( double ) { return 1.0; }, {}, p );
  }

  //evalTable: pairs at discontinuities, and clamped extrapolation:
  {
    const NC::VectD wl = { 0.0, 1.0, 2.0, 2.0, 3.0 };
    const NC::VectD xs = { 1.0, 2.0, 3.0, 1.0, 0.5 };
    auto ev = [&]( double w ) { return NCF::evalTable( wl.data(), xs.data(), wl.size(), w ); };
    nc_assert_always( ev( -1.0 ) == 1.0 );
    nc_assert_always( ev( 0.5 ) == 1.5 );
    nc_assert_always( ev( 1.999 ) > 2.99 );
    nc_assert_always( ev( 2.0 ) == 1.0 );
    nc_assert_always( ev( 2.5 ) == 0.75 );
    nc_assert_always( ev( 4.0 ) == 0.0 );//extrapolated: -0.0, clamped
    nc_assert_always( ev( 3.5 ) == 0.25 );
    std::cout << "evalTable: OK" << std::endl;
  }
  return 0;
}

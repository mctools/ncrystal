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

#include "NCrystal/internal/utils/NCSpline.hh"

#define REQUIRE(x) nc_assert_always(x)

namespace NC = NCrystal;

namespace {
  using PWLFct = NC::PiecewiseLinearFct1D;

  template <class TFunc>
  PWLFct createPWL( const NC::VectD& xvals, const TFunc& func, PWLFct::OutOfBoundsYValues ofvals = {} )
  {
    return PWLFct{ xvals, NC::vectorTrf( xvals, func ), ofvals };
  }

  template <class TestCodeFct>
  void ensureLogicError( const TestCodeFct& fct )
  {
    bool expectedError = false;
    try {
      fct();
    } catch ( NC::Error::LogicError& e ) {
      expectedError = true;
    }
    REQUIRE(expectedError);
  }

  template <class TestCodeFct>
  void ensureCalcError( const TestCodeFct& fct )
  {
    bool expectedError = false;
    try {
      fct();
    } catch ( NC::Error::CalcError& e ) {
      expectedError = true;
    }
    REQUIRE(expectedError);
  }
}

int main()
{
  PWLFct{ { 1.0, 2.0, 3.0 }, { 0.0, 0.0, 0.0 } };
#ifndef NDEBUG
  ensureLogicError([](){ PWLFct{ { 1.0, 2.0, NC::kInfinity }, { 0.0, 0.0, 0.0 } }; });
  ensureLogicError([](){ PWLFct{ { 1.0, 2.0, 3.0 }, { 0.0, 0.0, NC::kInfinity } }; });
#endif

  auto linfct = [](double u) { return 1.0-2.0*u; };
  {
    //first test with a function which is exactly reproduced with the piecewise
    //linear approximation:
    const auto xgrid = NC::linspace(-10.0,10.0,50);
    const auto fct_exact = createPWL( xgrid, linfct);
    REQUIRE( fct_exact.xValues() == xgrid );
    for ( auto& x : NC::linspace(xgrid.front(),xgrid.back(),1000) )
      REQUIRE( NC::floateq(linfct(x),fct_exact.eval(x),1e-14,1e-14) );
    for ( auto& x : xgrid )
      REQUIRE( linfct(x) == fct_exact.eval(x) );
    //overflow/underflow were not allowed:
    ensureCalcError([&fct_exact,&xgrid](){ fct_exact( xgrid.front() - 0.00001 ); });
    ensureCalcError([&fct_exact](){ fct_exact( -NC::kInfinity ); });
    ensureCalcError([&fct_exact,&xgrid](){ fct_exact( xgrid.back() + 0.00001 ); });
    ensureCalcError([&fct_exact](){ fct_exact( NC::kInfinity ); });
  }

  {
    //try with overflow but not underflow:
    const auto xgrid = NC::linspace(-10.0,10.0,50);
    const auto fct_exact = createPWL( xgrid, linfct, { NC::NullOpt, 2.6 } );
    for ( auto& x : NC::linspace(xgrid.front(),xgrid.back(),1000) )
      REQUIRE( NC::floateq(linfct(x),fct_exact.eval(x),1e-14,1e-14) );
    for ( auto& x : xgrid )
      REQUIRE( linfct(x) == fct_exact.eval(x) );
    ensureCalcError([&fct_exact,&xgrid](){ fct_exact( xgrid.front() - 0.00001 ); });
    ensureCalcError([&fct_exact](){ fct_exact( -NC::kInfinity ); });
    REQUIRE( 2.6 == fct_exact.eval(xgrid.back()+1e-13 ) );
    REQUIRE( 2.6 == fct_exact.eval(NC::kInfinity) );
  }

  {
    //try with underflow but not overflow:
    const auto xgrid = NC::linspace(-10.0,10.0,50);
    const auto fct_exact = createPWL( xgrid, linfct, { -3.4, NC::NullOpt } );
    for ( auto& x : NC::linspace(xgrid.front(),xgrid.back(),1000) )
      REQUIRE( NC::floateq(linfct(x),fct_exact.eval(x),1e-14,1e-14) );
    for ( auto& x : xgrid )
      REQUIRE( linfct(x) == fct_exact.eval(x) );
    ensureCalcError([&fct_exact,&xgrid](){ fct_exact( xgrid.back() + 0.00001 ); });
    ensureCalcError([&fct_exact](){ fct_exact( NC::kInfinity ); });
    REQUIRE( -3.4 == fct_exact.eval(xgrid.front()-1e-13 ) );
    REQUIRE( -3.4 == fct_exact.eval(-NC::kInfinity) );
  }


  {
    //try with both underflow and overflow:
    const auto xgrid = NC::linspace(-10.0,10.0,50);
    const auto fct_exact = createPWL( xgrid, linfct, { -3.4, 2.6 } );
    for ( auto& x : NC::linspace(xgrid.front(),xgrid.back(),1000) )
      REQUIRE( NC::floateq(linfct(x),fct_exact.eval(x),1e-14,1e-14) );
    for ( auto& x : xgrid )
      REQUIRE( linfct(x) == fct_exact.eval(x) );
    REQUIRE( -3.4 == fct_exact.eval(xgrid.front()-1e-13 ) );
    REQUIRE( -3.4 == fct_exact.eval(-NC::kInfinity) );
    REQUIRE( 2.6 == fct_exact.eval(xgrid.back()+1e-13 ) );
    REQUIRE( 2.6 == fct_exact.eval(NC::kInfinity) );
  }

  return 0;
}

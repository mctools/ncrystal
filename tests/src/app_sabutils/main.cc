////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

#include "NCrystal/internal/NCSABUtils.hh"
#include "NCrystal/NCFmt.hh"
#include <iostream>
namespace NC = NCrystal;
namespace NS = NC::SABUtils;

void printrange(NC::Span<double> v,bool eol = true) {
  //fixme: we should think about numpy-like adapters to print span's of values in ostreams. Perhaps even overloading adapters:
  //std::cout << NC::stream(v) <<std::endl;
  //std::cout << NC::stream_shortform(v) <<std::endl;
  std::cout<<"[";
  if (!v.empty()) {
    for (auto e: NC::Span<const double>(v).first(v.size()-1))
      std::cout<<e<<", ";
    std::cout<<v.back();
  }
  std::cout<<"]";
  if (eol)
    std::cout<<std::endl;
}

void test_loglin( double a, double fa, double b, double fb, double x, double calculated_value ) {
  double expected;
  std::ostringstream ss;
  ss << "loglin(a="<<a<<", fa="<<fa<<", b="<<b<<", fb="<<fb<<", x="<<x<<")";
  bool linlin(false);
  if (fa==fb) {
    expected = fa;//trivial
  } else if ( fa!=0 && fb!=0 ) {
    expected = std::exp(std::log(fa)+(std::log(fb/fa))*(x-a)/(b-a));//loglin definition
  } else {
    linlin = true;
    expected = fa+(fb-fa)*(x-a)/(b-a);//linlin
  }
  if (!NC::floateq(expected, calculated_value,1e-12,1e-12))
    NCRYSTAL_THROW2(CalcError,ss.str()<<" gave "<<calculated_value<<" not expected "<<expected<<" (absdiff "<<
                    calculated_value-expected<<")"<<(linlin?"linlin":"loglin")<<")");
  //std::cout<<"Testing "<<ss.str()<<" [expected "<<expected<<", got "<<calculated_value<<"] ("<<(linlin?"linlin":"loglin")<<")"<<std::endl;

}

void test_alpha_integrals() {
  //Reference values created with:
  //
  //#!/usr/bin/env python3
  //import numpy as np
  //import mpmath
  //mpmath.mp.dps = 1000
  //mpf=mpmath.mpf
  //def exact_eval_mpf(y):
  //    if y==0.0:
  //        return mpf(0.5)
  //    s1=(mpf(1)-mpf(y))/mpf(2)
  //    s2=(mpf(1)+mpf(y))/mpf(2)
  //    return (s2-s1)/mpmath.mp.log(s2/s1)
  //yvals = np.asarray([1e-200,1e-20,1e-10,1e-4,0.005,0.007,0.09,0.11,0.2,0.5,0.9,0.99,1.0-1e-14])
  //yvals = np.hstack((-np.flip(yvals),np.asarray([0.0]),yvals))
  //for y in yvals:
  //    print('   ref.emplace_back(',y,',',mpmath.nstr(exact_eval_mpf(y),30),');')
  //
  std::vector<NC::PairDD> ref;
  ref.emplace_back( -0.99999999999999 , 0.0303673187635422435892942080719 );
  ref.emplace_back( -0.99 , 0.18702871509983897278285159318 );
  ref.emplace_back( -0.9 , 0.305660944705597718728739460966 );
  ref.emplace_back( -0.5 , 0.455119613313418696807120082868 );
  ref.emplace_back( -0.2 , 0.493260692475286336445129433371 );
  ref.emplace_back( -0.11 , 0.497976784653061267607378554863 );
  ref.emplace_back( -0.09 , 0.49864707156245853926222023281 );
  ref.emplace_back( -0.007 , 0.499991833226619483216329399507 );
  ref.emplace_back( -0.005 , 0.499995833305555191792870951536 );
  ref.emplace_back( -0.0001 , 0.499999998333333328888888705869 );
  ref.emplace_back( -1e-10 , 0.499999999999999999998333333333 );
  ref.emplace_back( -1e-20 , 0.5 );
  ref.emplace_back( -1e-200 , 0.5 );
  ref.emplace_back( 0.0 , 0.5 );
  ref.emplace_back( 1e-200 , 0.5 );
  ref.emplace_back( 1e-20 , 0.5 );
  ref.emplace_back( 1e-10 , 0.499999999999999999998333333333 );
  ref.emplace_back( 0.0001 , 0.499999998333333328888888705869 );
  ref.emplace_back( 0.005 , 0.499995833305555191792870951536 );
  ref.emplace_back( 0.007 , 0.499991833226619483216329399507 );
  ref.emplace_back( 0.09 , 0.49864707156245853926222023281 );
  ref.emplace_back( 0.11 , 0.497976784653061267607378554863 );
  ref.emplace_back( 0.2 , 0.493260692475286336445129433371 );
  ref.emplace_back( 0.5 , 0.455119613313418696807120082868 );
  ref.emplace_back( 0.9 , 0.305660944705597718728739460966 );
  ref.emplace_back( 0.99 , 0.18702871509983897278285159318 );
  ref.emplace_back( 0.99999999999999 , 0.0303673187635422435892942080719 );

  //for (auto& [ y, refval ] : ref ) {
  for (auto& e : ref ) {
    // auto& y = e.first;
    // auto& refval = e.second;
    const double y = e.first;
    const double refval = e.second;
    const double s1 = (1.0-y)/2.0;
    const double s2 = (1.0+y)/2.0;
    const double res = NS::integrateAlphaInterval(0.0, s1, 1.0, s2);//claims precision better than 1e-15
    const double resfast = NS::integrateAlphaInterval_fast(0.0, s1, 1.0, s2, std::log(s1), std::log(s2) );//claims precision better than 1e-14
    if (!NC::floateq(res,refval,1e-15,0.0))
      NCRYSTAL_THROW2(LogicError,"problem with integrateAlphaInterval for y="
                      <<NC::fmt(y)<<" : "<<res<<" differs too much from expected "<<refval);
    if (!NC::floateq(resfast,refval,1e-14,0.0))
      NCRYSTAL_THROW2(LogicError,"problem with integrateAlphaInterval_fast for y="
                      <<NC::fmt(y)<<" : "<<resfast<<" differs too much from expected "<<refval);
  }
}

int main()
{
  {
    //Unit test alpha limits taylor expansion (really for a dedicated KinUtils test)
    const double ekin_div_kT = 6.345703125e-08;//probably doesn't matter (using value where we once saw a failure)
    const double betalim=0.01*ekin_div_kT;//NB: 0.01 should be same factor as used in alphaMinusNeedsTaylor!
    const double b0 = betalim*(1+1e-15);
    const double b1 = betalim*(1-1e-15);
    nc_assert_always( NC::detail::alphaMinusNeedsTaylor(ekin_div_kT,b0) != NC::detail::alphaMinusNeedsTaylor(ekin_div_kT,b1) );
    nc_assert_always( NC::detail::alphaMinusNeedsTaylor(ekin_div_kT,-b0) != NC::detail::alphaMinusNeedsTaylor(ekin_div_kT,-b1) );
    nc_assert_always( NC::floateq( NC::getAlphaLimits(ekin_div_kT,b0).first, NC::getAlphaLimits(ekin_div_kT,b1).first, 1e-11, 1e-299 ) );
    nc_assert_always( NC::floateq( NC::getAlphaLimits(ekin_div_kT,-b0).first, NC::getAlphaLimits(ekin_div_kT,-b1).first, 1e-11, 1e-299 ) );
  }

  ////////////////////////

  using V = NCrystal::VectD;

  ////////////////////////

  V outb, outsab;
  auto insab = V{  0.111, 0.333, 0.777,
                   1.111, 1.333, 1.777,
                   4.111, 4.333, 4.777 };
  NS::expandBetaAndSABToAllBetas( NC::Span<const double>(V{0,1,4}),
                                  V{ 11.1, 33.3, 77.7 },
                                  insab,
                                  outb,outsab);

  nc_assert_always(outb==V({-4, -1, 0, 1, 4}));
  nc_assert_always(outsab==V({4.111, 4.333, 4.777, 1.111, 1.333, 1.777, 0.111, 0.333, 0.777, 1.111, 1.333, 1.777, 4.111, 4.333, 4.777}));
  nc_assert_always(V(outsab.begin()+6,outsab.end()) == insab);
  printrange(outb);
  printrange(outsab);

  //////////////////////////

  for (auto& a :  { -1e5, -1e-5, 0.0, 0.5, 10.0 ,1e5 } ) {
    for (auto& brel :  { 0.001, 0.9, 1000.0 } ) {
      double b = a + std::abs(a?a:1.0)*brel;
      for (auto& fa :  { 0.0, 1e-5, 1.0, 1e5 } ) {
        for (auto& fb :  { 0.0, 1e-5, 1.0, 1e5 } ) {
          for (auto& xrel : { 0.0, 0.1, 0.5, 0.8, 1.0 } ) {
            double x = (xrel == 1.0 ? b : a + xrel * (b-a) );
            test_loglin( a,fa,b,fb,x, NS::interpolate_loglin_fallbacklinlin( a,fa,b,fb,x ) );
            double logfa(fa>0.0?std::log(fa):-1000.0);
            double logfb(fb>0.0?std::log(fb):-1000.0);
            test_loglin( a,fa,b,fb,x, NS::interpolate_loglin_fallbacklinlin_fast( a,fa,b,fb,x, logfa, logfb ) );
          }
        }
      }
    }
  }

  //////////////////////////
  test_alpha_integrals();

}

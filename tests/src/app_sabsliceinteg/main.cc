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

//Test createTailedBreakdown function by integrating simple function,
//SAB=exp(alpha), which is interpolated exactly by the loglin-interpolation
//assumed for SAB tables (at least internally in NCSABUtils).

#include "NCrystal/internal/NCSABUtils.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCIter.hh"
#include <iostream>
namespace NC = NCrystal;

namespace {

  struct FakeSABSlice {
    const NC::VectD alphaGrid, sab, logsab, alphaIntegrals_cumul;
  };

  template<class Func_F, class Func_DefIntegralOfF>
  FakeSABSlice createFakeSABSlice( const NC::VectD& alphaGrid_,
                                   const Func_F& f,
                                   const Func_DefIntegralOfF&  f_defintegral)
  {
    NC::VectD sab, logsab, aic;
    for ( auto&& alpha : NC::enumerate(alphaGrid_) ) {
      sab.push_back(f(alpha.val));
      nc_assert_always(sab.back());
      logsab.push_back(sab.back()>0.0?std::log(sab.back()):-NC::kInfinity);
    }
    aic.push_back(0.0);
    for (unsigned i=1; i < alphaGrid_.size();++i)
      aic.push_back(aic.back()+f_defintegral(alphaGrid_.at(i-1),alphaGrid_.at(i)));
    return FakeSABSlice{alphaGrid_,std::move(sab),std::move(logsab),std::move(aic)};
  }

  NC::SABUtils::TailedBreakdown
  makeTailedBreakdown(const FakeSABSlice& fakesabslice,
                      const double alpha_low,
                      const double alpha_upp)
  {
    const auto& ag = fakesabslice.alphaGrid;
    auto it = std::lower_bound(ag.begin(),ag.end(),alpha_low);
    if (it > ag.begin())
      it = std::prev(it);
    while (std::next(it)!=ag.end() && *std::next(it)<=alpha_low )
      it = std::next(it);
    unsigned aidx_low = std::distance(ag.begin(),it);
    it = std::lower_bound(ag.begin(),std::prev(ag.end()),alpha_upp);
    while (it>ag.begin() && *std::prev(it)>=alpha_upp )
      it = std::prev(it);
    unsigned aidx_upp = std::distance(ag.begin(),it);
    nc_assert_always(aidx_low<=aidx_upp && aidx_upp < ag.size());
    std::cout<<"makeTailedBreakdown DBG: aidx_low="
             <<aidx_low<<", aidx_upp="<<aidx_upp<<std::endl;
    return NC::SABUtils::createTailedBreakdown( fakesabslice.alphaGrid,
                                                fakesabslice.sab,
                                                fakesabslice.logsab,
                                                fakesabslice.alphaIntegrals_cumul,
                                                alpha_low, alpha_upp,
                                                aidx_low, aidx_upp );
  }

  template<class Func_DefIntegralOfF>
  void checkXS( const FakeSABSlice& fakesabslice,
                double a0,
                double a1,
                const Func_DefIntegralOfF& func_f_defintegral_)
  {
    auto tb = makeTailedBreakdown(fakesabslice,a0,a1);
    double xs_sum(tb.xs_front+tb.xs_middle+tb.xs_back), xs_expect(func_f_defintegral_(a0,a1));
    static int i = 0;
    double reldiff = (xs_expect==xs_sum?0.0:(xs_sum/xs_expect-1));
    std::cout<<" checkXS("<<++i<<"): xs_front = "<<tb.xs_front
             <<" xs_middle = "<<tb.xs_middle
             <<" xs_back = "<<tb.xs_back
             <<" xs_sum = "<<xs_sum
             <<" expected = "<<xs_expect
           //<<" reldiff = "<<100.0*reldiff<<"%"
             <<" reldiff<1e-15? = "<<(reldiff<1e-15?"yes":"no")<<"%"
             <<std::endl;
  }

  template<class Func_F, class Func_DefIntegralOfF>
  void check_integrateAlphaInterval (double a0, double a1,
                                     const Func_F& func_f_,
                                     const Func_DefIntegralOfF&  func_f_defintegral_)
  {
    double f0(func_f_(a0)),f1(func_f_(a1));
    double calc = NC::SABUtils::integrateAlphaInterval(a0,f0, a1, f1 );
    double calc_fast = NC::SABUtils::integrateAlphaInterval_fast(a0,f0, a1, f1, std::log(f0), std::log(f1) );
    double expect = func_f_defintegral_(a0,a1);
    constexpr double tol = 3e-14;
    auto rd_ok = [expect](double v)
    {
      return ( ( v==expect || NC::ncabs(v/expect-1.0) < tol )
               ? "yes" : "no" );
    };

    std::cout<<" integrateAlphaInterval : "<<calc<<" vs. "<<calc_fast<<" vs. "<<expect<<" reldiff < "<<tol<<": "
             <<rd_ok(calc)<<" and "<<rd_ok(calc_fast)<<std::endl;
  }
}

int main () {

#if 0
  auto alphaGrid = NC::linspace(0.1,2.0,20);
#elif 0
  //This emulates issue seen on osx:
  const NC::VectD alphaGrid = { 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
                                1.1,1.2,1.3,1.4,1.5,1.6,1.7,
                                std::nextafter(1.8, -NC::kInfinity),
                                1.9,2.0 };
#else
  //Hardcode input for now, to test sabutils and not NC::linspace:
  const NC::VectD alphaGrid = { 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
                                1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0 };
#endif

  auto func_f = [](double alpha) { return std::exp(alpha); };
  auto func_f_defintegral = [&alphaGrid](double a, double b)
                            {
                              //Implement cutoff on alphaGrid range:
                              nc_assert_always(b>=a);
                              if (a>=alphaGrid.back()||b<=alphaGrid.front())
                                return 0.0;
                              b = std::min(b,alphaGrid.back());
                              a = std::max(a,alphaGrid.front());
                              nc_assert_always(b>=a);
                              //Evaluate:
                              //return std::exp(b)-std::exp(a);//numerically unstable
                              return std::exp(a)*std::expm1(b-a);//numerically stable
                            };
  auto sabslice = createFakeSABSlice(alphaGrid,func_f,func_f_defintegral);
  checkXS(sabslice, 0.1,2.0,func_f_defintegral);//everything
  checkXS(sabslice, 0.2,1.8,func_f_defintegral);//a range of whole bins
  checkXS(sabslice, 0.1123,1.8,func_f_defintegral);//+front
  checkXS(sabslice, 0.112,1.8123,func_f_defintegral);//+back
  checkXS(sabslice, 0.123,1.8123,func_f_defintegral);//+front and back
  checkXS(sabslice, 0.15,1.85,func_f_defintegral);
  checkXS(sabslice, 0.1234,0.1567,func_f_defintegral);
  checkXS(sabslice, 0.20000001,0.2567,func_f_defintegral);
  checkXS(sabslice, 0.2,0.2567,func_f_defintegral);
  checkXS(sabslice, 0.2,0.5567,func_f_defintegral);
  checkXS(sabslice, 0.2,1.4,func_f_defintegral);

  checkXS(sabslice, 100.0,101.0,func_f_defintegral);//outside
  checkXS(sabslice, 1e-10,1e-9,func_f_defintegral);//outside

  checkXS(sabslice, 0.05, 1.123,func_f_defintegral);//partly outside
  checkXS(sabslice, 0.05, 2.3,func_f_defintegral);//partly outside
  checkXS(sabslice, 1.05, 2.3,func_f_defintegral);//partly outside

  //inside:
  check_integrateAlphaInterval(0.11,1.22, func_f, func_f_defintegral);
  check_integrateAlphaInterval(0.1,2.0, func_f, func_f_defintegral);

  // //Extending outside:
  // check_integrateAlphaInterval(0.01,1.2, func_f, func_f_defintegral);
  // check_integrateAlphaInterval(0.0,1.2, func_f, func_f_defintegral);

}

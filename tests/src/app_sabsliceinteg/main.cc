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
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCIter.hh"
#include <iostream>
namespace NC = NCrystal;

// struct TailedBreakdown {
//       double xs_front=0, xs_middle=0, xs_back=0;
//       unsigned imiddle_low=0, imiddle_upp=0;
//       struct TailPoint {
//         double alpha=0, sval=0, logsval=0;
//       } front, back;
//       bool narrow = false;
//     };
//     TailedBreakdown createTailedBreakdown( const span<const double>& alphaGrid,
//                                            const span<const double>& sab,
//                                            const span<const double>& logsab,
//                                            const span<const double>& alphaIntegrals_cumul,
//                                            const double alpha_low, const double alpha_upp,
//                                            const unsigned aidx_low, const unsigned aidx_upp );

// namespace NCTest {
//   struct FakeSABSlice {
//     const vector<const double> alphaGrid, sab, logsab, alphaIntegrals_cumul;
//   };
//}

int main () {

  //Test createTailedBreakdown function by integrating simple function,
  //SAB=exp(alpha), which is interpolated exactly by the loglin-interpolation
  //assumed for SAB tables (at least internally in NCSABUtils).

  struct FakeSABSlice {
    const NC::VectD alphaGrid, sab, logsab, alphaIntegrals_cumul;
  };
  auto createFakeSABSlice = [](auto alphaGrid_, auto f, auto f_defintegral)
                            {
                              NC::VectD sab, logsab, aic;
                              for ( auto&& alpha : NC::enumerate(alphaGrid_) ) {
                                sab.push_back(f(alpha.val));
                                nc_assert(sab.back());
                                logsab.push_back(sab.back()>0.0?std::log(sab.back()):-NC::kInfinity);
                              }
                              aic.push_back(0.0);
                              for (unsigned i=1; i < alphaGrid_.size();++i)
                                aic.push_back(aic.back()+f_defintegral(alphaGrid_.at(i-1),alphaGrid_.at(i)));
                              return FakeSABSlice{alphaGrid_,std::move(sab),std::move(logsab),std::move(aic)};
                            };

  auto makeTailedBreakdown = [](const auto& fakesabslice, const double alpha_low, const double alpha_upp)
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
                               nc_assert(aidx_low<=aidx_upp && aidx_upp < ag.size());
                               return NC::SABUtils::createTailedBreakdown( fakesabslice.alphaGrid,
                                                                           fakesabslice.sab,
                                                                           fakesabslice.logsab,
                                                                           fakesabslice.alphaIntegrals_cumul,
                                                                           alpha_low, alpha_upp,
                                                                           aidx_low, aidx_upp );
                             };

  auto checkXS = [&makeTailedBreakdown](const auto& fakesabslice, double a0, double a1, auto func_f_defintegral_)
                 {
                   auto tb = makeTailedBreakdown(fakesabslice,a0,a1);
                   double xs_sum(tb.xs_front+tb.xs_middle+tb.xs_back), xs_expect(func_f_defintegral_(a0,a1));
                   static int i = 0;
                   std::cout<<" checkXS("<<++i<<"): xs_front = "<<tb.xs_front
                            <<" xs_middle = "<<tb.xs_middle
                            <<" xs_back = "<<tb.xs_back
                            <<" xs_sum = "<<xs_sum
                            <<" expected = "<<xs_expect
                            <<" reldiff = "<<100.0*(xs_expect==xs_sum?0.0:(xs_sum/xs_expect-1))<<"%"
                            <<std::endl;
                 };

  //
  auto alphaGrid = NC::linspace(0.1,2.0,20);
  auto func_f = [](double alpha) { return std::exp(alpha); };
  auto func_f_defintegral = [&alphaGrid](double a, double b)
                            {
                              //Implement cutoff on alphaGrid range:
                              nc_assert(b>=a);
                              if (a>=alphaGrid.back()||b<=alphaGrid.front())
                                return 0.0;
                              b = std::min(b,alphaGrid.back());
                              a = std::max(a,alphaGrid.front());
                              nc_assert(b>=a);
                              //Evaluate:
                              return std::exp(b)-std::exp(a);

                            };
  auto sabslice = createFakeSABSlice(NC::linspace(0.1,2.0,20),func_f,func_f_defintegral);
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


  auto check_integrateAlphaInterval = [](double a0, double a1, auto func_f_, auto func_f_defintegral_)
                                      {
                                        double f0(func_f_(a0)),f1(func_f_(a1));
                                        double calc = NC::SABUtils::integrateAlphaInterval(a0,f0, a1, f1 );
                                        double calc_fast = NC::SABUtils::integrateAlphaInterval_fast(a0,f0, a1, f1, std::log(f0), std::log(f1) );
                                        double expect = func_f_defintegral_(a0,a1);
                                        std::cout<<" integrateAlphaInterval : "<<calc<<" vs. "<<calc_fast<<" vs. "<<expect<<" reldiff: "
                                                 <<100.0*(calc/expect-1.0)<<"% and "<<100.0*(calc_fast/expect-1.0)<<"%"<<std::endl;
                                      };
  //inside:
  check_integrateAlphaInterval(0.11,1.22, func_f, func_f_defintegral);
  check_integrateAlphaInterval(0.1,2.0, func_f, func_f_defintegral);

  // //Extending outside:
  // check_integrateAlphaInterval(0.01,1.2, func_f, func_f_defintegral);
  // check_integrateAlphaInterval(0.0,1.2, func_f, func_f_defintegral);

}

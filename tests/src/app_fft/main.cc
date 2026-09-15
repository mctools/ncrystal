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

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCFastConvolve.hh"
#include <iostream>
#include "refvals.hh"

namespace NC=NCrystal;

#define REQUIRE(x) nc_assert_always(x)
#define REQUIREFLTEQ(x, y) nc_assert_always(::NC::floateq((x), (y)))

namespace {

  // Tests linear convolution of two real VectD inputs.
  // The output has size a1.size()+a2.size()-1.
  // dt scales every output element.
  // Negative values must remain negative; this detects abs() misuse.
  // Inputs whose FFT size is not a power of two are also tested.
  // Empty-input behaviour is intentionally not tested here.

  using NC::VectD;
  using NC::vectAt;
  using NC::fmt;
  using NC::ncrange;

  void check(const VectD& got, const VectD& want)
  {
    auto prtvect = [](const VectD&v)
    {
      nc_assert_always(!v.empty());
      std::cout<<"[";
      for ( auto i : ncrange(v.size()) ) {
        double val = v.at(i);
        if (val==0.0)
          val = 0.0;//no -0
        if ( i )
          std::cout<<", ";
        std::cout<<fmt(val,"%.12g");
      }
      std::cout<<"]"<<std::endl;
    };
    prtvect(got);
    prtvect(want);
#ifndef NDEBUG
    REQUIRE(got.size() == want.size());
#endif
    auto i = std::size_t(0);
    for ( auto f : want ) {
      REQUIREFLTEQ(vectAt(got,i), f);
      ++i;
    }
  }

  void testFastConvolve()
  {
    NC::FastConvolve c;
    VectD y;

    c.convolve(VectD{1.0}, VectD{-2.0}, y, 1.0);
    check(y, VectD{-2.0});
    c.convolve(VectD{1.0,-2.0}, VectD{3.0,4.0}, y, 1.0);
    check(y, VectD{3.0,-2.0,-8.0});
    c.convolve(VectD{1.0,2.0,3.0},
               VectD{4.0,5.0}, y, 1.0);
    check(y, VectD{4.0,13.0,22.0,15.0});
    c.convolve(VectD{1.0,2.0},
               VectD{3.0,4.0}, y, 0.25);
    check(y, VectD{0.75,2.5,2.0});
    c.convolve(VectD{1.0,2.0,3.0,4.0},
               VectD{2.0,-1.0,0.5}, y, 1.0);
    check(y, VectD{2.0,3.0,4.5,6.0,-2.5,2.0});
    c.convolve(VectD{0.0,0.0,0.0},
               VectD{2.0,-3.0}, y, 1.0);
    check(y, VectD{0.0,0.0,0.0,0.0});
  }

  void testCalcPhase()
  {
    using NC::ncabs;
    using NC::ncmax;
    using NC::FastConvolve;
    using NC::fmt;
    using NC::fmtg;
    auto rd = [](double val, double stdval)
    {
      if ( val == 0.0 && ncabs(stdval)<1e-15 )
        return 0.0;//std library does not always return strict 0 like we do
      return ncabs(val-stdval)/(ncmax(1e-30,0.5*(ncabs(val)+ncabs(stdval))));
    };
    double statn_worstrd;
    unsigned long statn_currentn;
    unsigned long statn_kvals;
    auto reset_statn = [&statn_worstrd,&statn_kvals,&statn_currentn](){
      statn_worstrd = -1.0;
      statn_kvals = 0;
      statn_currentn = std::numeric_limits<unsigned long>::max();
    };
    auto add_statn = [&statn_worstrd,&statn_kvals](double rdval)
    {
      ++statn_kvals;
      statn_worstrd = ncmax(statn_worstrd,rdval);
    };
    auto allowed_rd = [](unsigned long n)
    {
      if ( n<=3 )
        return 0.0;
      if ( n<=8 )
        return 9e-16;
      return 9e-15;
    };
    auto flush_statn = [&statn_worstrd,&statn_kvals,
                        &statn_currentn,&reset_statn,&allowed_rd](){
      const auto n = statn_currentn;
      if ( n == std::numeric_limits<unsigned long>::max() )
        return;
      std::cout<<"Number of k values tested for n="<<n
               <<": "<<statn_kvals;
      if ( statn_kvals == (1ul<<n) )
        std::cout<<" [ALL]";
      const double eps = allowed_rd(n);
      const bool ok = statn_worstrd<=eps;
      if ( ok )
        std::cout<<" (worst reldiff < "<<fmtg(eps)<<": yes)"<<std::endl;
      else
        std::cout<<" (worst reldiff above "<<fmtg(eps)<<": "
                 <<fmt(statn_worstrd)<<')'<<std::endl;
      nc_assert_always(ok);
      reset_statn();
    };

    reset_statn();

    double maxrd(0.0);
    for ( auto itest : NC::ncrange(calcPhaseTests::ntests) ) {
      const auto n = calcPhaseTests::nvals[itest];
      if ( statn_currentn != n ) {
        flush_statn();
        statn_currentn = n;
        maxrd = allowed_rd(n);
      }
      const auto k = calcPhaseTests::kvals[itest];
      const double ref_cosval = calcPhaseTests::cosvals[itest];
      const double ref_sinval = calcPhaseTests::sinvals[itest];
      auto res = FastConvolve::calcPhase(k,n);
      const double rd_c = rd(res.first,ref_cosval);
      const double rd_s = rd(res.second,ref_sinval);
      add_statn(ncmax(rd_c,rd_s));
      const bool ok = ( ncmax(rd_c,rd_s) <= maxrd
                        && ncabs(res.first-ref_cosval)<1e-14
                        && ncabs(res.second-ref_sinval)<1e-14 );
      if (!ok) {
        std::cout<<"Real[Phase(k="<<k<<", n="<<n<<")] = "<<fmt(res.first)
                 <<" (ref calc gives "<<fmt(ref_cosval)
                 <<", reldiff="<<fmtg(rd_c)<<")."<<std::endl;
        std::cout<<"Imag[Phase(k="<<k<<", n="<<n<<")] = "<<fmt(res.second)
                 <<" (ref calc gives "<<fmt(ref_sinval)
                 <<", reldiff="<<fmtg(rd_s)<<")."<<std::endl;
      }
      REQUIRE(ok);
    }
    flush_statn();
  }
}

int main() {
  std::cout<<"testFastConvolve start..."<<std::endl;
  testFastConvolve();
  std::cout<<"testFastConvolve done."<<std::endl;
  std::cout<<"testCalcPhase start..."<<std::endl;
  testCalcPhase();
  std::cout<<"testCalcPhase done."<<std::endl;
  return 0;
}

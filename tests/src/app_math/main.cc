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

#include "NCrystal/internal/NCMath.hh"
#include <iostream>
namespace NC=NCrystal;

template <class Func>
void testIntegration(const Func& f, double a, double b, double expected_integral, std::string descr) {
  std::cout<<"\nTest integral of f(x)="<<descr<<" over ["<<a<<", "<<b<<"]:"<<std::endl;
  std::cout<<"  Exact                : "<<expected_integral<<std::endl;
  double val;
  val = NC::integrateTrapezoidal(f,a,b,5);
  std::cout<<"  Trapezoidal(n=5)     : "<<val<<" (relative deviation: "<<(val/expected_integral-1.0)<<")"<<std::endl;
  val = NC::integrateTrapezoidal(f,a,b,10000);
  std::cout<<"  Trapezoidal(n=10000) : "<<val<<" (relative deviation: "<<(val/expected_integral-1.0)<<")"<<std::endl;
  val = NC::integrateSimpsons(f,a,b,2);
  std::cout<<"  Simpsons(n=2)        : "<<val<<" (relative deviation: "<<(val/expected_integral-1.0)<<")"<<std::endl;
  val = NC::integrateSimpsons(f,a,b,10000);
  std::cout<<"  Simpsons(n=10000)    : "<<val<<" (relative deviation: "<<(val/expected_integral-1.0)<<")"<<std::endl;

  val = NC::integrateTrapezoidal(f,a,b,17);
  std::cout<<"  Trapezoidal(n=17)     : "<<val<<" (relative deviation: "<<(val/expected_integral-1.0)<<")"<<std::endl;
  val = NC::integrateSimpsons(f,a,b,17);
  std::cout<<"  Simpsons(n=17)        : "<<val<<" (relative deviation: "<<(val/expected_integral-1.0)<<")"<<std::endl;
  val = NC::integrateRomberg17(f,a,b);
  std::cout<<"  Romberg(n=17)        : "<<val<<" (relative deviation: "<<(val/expected_integral-1.0)<<")"<<std::endl;


  val = NC::integrateTrapezoidal(f,a,b,33);
  std::cout<<"  Trapezoidal(n=33)     : "<<val<<" (relative deviation: "<<(val/expected_integral-1.0)<<")"<<std::endl;
  val = NC::integrateSimpsons(f,a,b,33);
  std::cout<<"  Simpsons(n=33)        : "<<val<<" (relative deviation: "<<(val/expected_integral-1.0)<<")"<<std::endl;
  val = NC::integrateRomberg33(f,a,b);
  std::cout<<"  Romberg(n=33)        : "<<val<<" (relative deviation: "<<(val/expected_integral-1.0)<<")"<<std::endl;
}

int main() {

  testIntegration([](double x){return x*x+0.1;}, -1,2, 3.3, "x^2+1/10");
  testIntegration([](double x){return std::exp(-x*x)*x*x*x;}, 0,2,0.5-2.5*std::exp(-4.0), "exp(-x^2)*x^3 ");

  for( auto val : {0.5, 1.0, 2.0, 0.0, NC::kInfinity} ) {
    std::cout<<"NC::ekin2wl("<<val<<") = "<<NC::ekin2wl(val)<<std::endl;
    std::cout<<"NC::ekin2wlsq("<<val<<") = "<<NC::ekin2wlsq(val)<<std::endl;
    std::cout<<"NC::ekin2wlsqinv("<<val<<") = "<<NC::ekin2wlsqinv(val)<<std::endl;
    std::cout<<"NC::constexpr_ekin2wl("<<val<<") = "<<NC::constexpr_ekin2wl(val)<<std::endl;
    std::cout<<"NC::wl2ekin("<<val<<") = "<<NC::wl2ekin(val)<<std::endl;
    std::cout<<"NC::ekin2ksq("<<val<<") = "<<NC::ekin2ksq(val)<<std::endl;
    std::cout<<"NC::ekin2k("<<val<<") = "<<NC::ekin2k(val)<<std::endl;
    std::cout<<"NC::constexpr_ekin2k("<<val<<") = "<<NC::constexpr_ekin2k(val)<<std::endl;
    std::cout<<"NC::wl2k("<<val<<") = "<<NC::wl2k(val)<<std::endl;
    std::cout<<"NC::wl2ksq("<<val<<") = "<<NC::wl2ksq(val)<<std::endl;

    std::cout<<"NC::wlsq2ekin("<<val<<") = "<<NC::wlsq2ekin(val)<<std::endl;
    std::cout<<"NC::ksq2ekin("<<val<<") = "<<NC::ksq2ekin(val)<<std::endl;
    std::cout<<"NC::k2ekin("<<val<<") = "<<NC::k2ekin(val)<<std::endl;
    std::cout<<"NC::k2wl("<<val<<") = "<<NC::k2wl(val)<<std::endl;

  }

  std::cout<<"linspace test: ";
  for ( auto val : NC::linspace(0.1,2.0,20) )
    std::cout<<", "<<NC::fmt(val,"%.15g");
  std::cout<<std::endl;

}

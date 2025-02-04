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
#include "NCrystal/internal/utils/NCSpan.hh"
#include <iostream>
#include <array>

namespace NC = NCrystal;

namespace {
  void testfct1(NC::Span<double> sp ) { for (auto e : sp ) std::cout << e << std::endl; }
  void testfct2(const NC::Span<double>& sp ) { for (auto e : sp ) std::cout << e << std::endl; }
  void testfct3(NC::Span<const double> sp ) { for (auto e : sp ) std::cout << e << std::endl; }
  void testfct4(const NC::Span<const double>& sp ) { for (auto e : sp ) std::cout << e << std::endl; }

}

int main() {
  NC::VectD v = { 1.0, 1.0, 2.0, 3.0, 5.0, 8.0, 11.0 };
  std::array<double,4> a = { 0.1, 0.2, 0.3, 0.4 };
  double ca[3] = { 0.01, 0.02, 0.03 };

  NC::Span<double> span(v);
  std::cout<<"---> vector v:"<<std::endl;;
  for (auto e : v )
    std::cout<<e<<std::endl;

  std::cout<<"---> span(v):"<<std::endl;;
  for (auto e : span )
    std::cout<<e<<std::endl;

  std::cout<<"---> tmp span(v):"<<std::endl;;
  for (auto e : NC::Span<double>(v) )
    std::cout<<e<<std::endl;

  std::cout<<"---> tmp span(v).subspan(3,2):"<<std::endl;;
  for (auto e : NC::Span<double>(v).subspan(3,2) )
    std::cout<<e<<std::endl;

  std::cout<<"---> span(v).subspan(2,2):"<<std::endl;;
  for (auto e : span.subspan(2,2) )
    std::cout<<e<<std::endl;

  std::cout<<"---> span(v).subspan(2,20):"<<std::endl;;
  for (auto e : span.subspan(2,20) )
    std::cout<<e<<std::endl;

  std::cout<<"---> span(v).subspan(2,5):"<<std::endl;;
  for (auto e : span.subspan(2,5) )
    std::cout<<e<<std::endl;

  std::cout<<"---> span(v).subspan(0,7):"<<std::endl;;
  for (auto e : span.subspan(0,7) )
    std::cout<<e<<std::endl;

  std::cout<<"---> span(v).subspan(0,20):"<<std::endl;;
  for (auto e : span.subspan(0,20) )
    std::cout<<e<<std::endl;

  std::cout<<"---> span(v).subspan(0,5):"<<std::endl;;
  for (auto e : span.subspan(0,5) )
    std::cout<<e<<std::endl;

  std::cout<<"---> tmp span(v).subspan(2,0):"<<std::endl;;
  for (auto e : NC::Span<double>(v).subspan(2,0) )
    std::cout<<e<<std::endl;

  std::cout<<"---> tmp span(v).first(0):"<<std::endl;;
  for (auto e : NC::Span<double>(v).first(0) )
    std::cout<<e<<std::endl;
  std::cout<<"---> tmp span(v).first(2):"<<std::endl;;
  for (auto e : NC::Span<double>(v).first(2) )
    std::cout<<e<<std::endl;
  std::cout<<"---> tmp span(v).first(7):"<<std::endl;;
  for (auto e : NC::Span<double>(v).first(7) )
    std::cout<<e<<std::endl;
  std::cout<<"---> tmp span(v).first(99):"<<std::endl;;
  for (auto e : NC::Span<double>(v).first(99) )
    std::cout<<e<<std::endl;

  std::cout<<"---> tmp span(v).last(0):"<<std::endl;;
  for (auto e : NC::Span<double>(v).last(0) )
    std::cout<<e<<std::endl;
  std::cout<<"---> tmp span(v).last(2):"<<std::endl;;
  for (auto e : NC::Span<double>(v).last(2) )
    std::cout<<e<<std::endl;
  std::cout<<"---> tmp span(v).last(7):"<<std::endl;;
  for (auto e : NC::Span<double>(v).last(7) )
    std::cout<<e<<std::endl;
  std::cout<<"---> tmp span(v).last(99):"<<std::endl;;
  for (auto e : NC::Span<double>(v).last(99) )
    std::cout<<e<<std::endl;

  //Reverse iterators:
  std::cout<<"---> reversed iterators of span(v):"<<std::endl;;
  auto sss = NC::Span<double>(v);
  for (auto it = sss.rbegin(); it != sss.rend(); ++it)
    std::cout<<*it<<std::endl;

  std::cout<<"---> const reversed iterators of span(v):"<<std::endl;;
  auto sss2 = NC::Span<const double>(v);
  for (auto it = sss2.rbegin(); it != sss2.rend(); ++it)
    std::cout<<*it<<std::endl;

  std::cout<<"---> tmp span(<empty>).last(2):"<<std::endl;;
  for (auto e : NC::Span<double>().last(2) )
    std::cout<<e<<std::endl;
  std::cout<<"---> tmp span(<empty>).last(0):"<<std::endl;;
  for (auto e : NC::Span<double>().last(0) )
    std::cout<<e<<std::endl;

  std::cout<<"---> tmp span(<empty>).first(2):"<<std::endl;;
  for (auto e : NC::Span<double>().first(2) )
    std::cout<<e<<std::endl;
  std::cout<<"---> tmp span(<empty>).first(0):"<<std::endl;;
  for (auto e : NC::Span<double>().first(0) )
    std::cout<<e<<std::endl;

  std::cout<<"---> tmp span(<empty>).subspan(0,1):"<<std::endl;;
  for (auto e : NC::Span<double>().subspan(0,1) )
    std::cout<<e<<std::endl;


  std::cout<<"---> span(span(v)):"<<std::endl;;
  auto span2 = span;
  for (auto e : span2 )
    std::cout<<e<<std::endl;

  std::cout<<"---> span(move(span(v))):"<<std::endl;;
  auto span3 = std::move(span);
  for (auto e : span3 )
    std::cout<<e<<std::endl;

  std::cout<<"---> array a:"<<std::endl;;
  for (auto e : a )
    std::cout<<e<<std::endl;

  NC::Span<double> spana(a);
  std::cout<<"---> span(a):"<<std::endl;;
  for (auto e : spana )
    std::cout<<e<<std::endl;

  std::cout<<"---> c-array ca:"<<std::endl;;
  for (auto e : ca )
    std::cout<<e<<std::endl;

  std::cout<<"---> span(ca):"<<std::endl;;

  NC::Span<double> spanca(ca);
  for (auto e : spanca )
    std::cout<<e<<std::endl;

  //assignments:
  std::cout<<"---> span = span (=v):"<<std::endl;;
  spanca = span;
  for (auto e : spanca )
    std::cout<<e<<std::endl;
  spanca = v;
  std::cout<<"---> span = v:"<<std::endl;;
  for (auto e : spanca )
    std::cout<<e<<std::endl;
  spanca = ca;
  std::cout<<"---> span = ca:"<<std::endl;;
  for (auto e : spanca )
    std::cout<<e<<std::endl;
  spanca = a;
  std::cout<<"---> span = a:"<<std::endl;;
  for (auto e : spanca )
    std::cout<<e<<std::endl;

  //Check various function calls compiles:

  std::cout<<"Through fct1 call (v):"<<std::endl; testfct1(v);
  std::cout<<"Through fct2 call (v):"<<std::endl; testfct2(v);
  std::cout<<"Through fct3 call (v):"<<std::endl; testfct3(v);
  std::cout<<"Through fct4 call (v):"<<std::endl; testfct4(v);
  std::cout<<"Through fct1 call (a):"<<std::endl; testfct1(a);
  std::cout<<"Through fct2 call (a):"<<std::endl; testfct2(a);
  std::cout<<"Through fct3 call (a):"<<std::endl; testfct3(a);
  std::cout<<"Through fct4 call (a):"<<std::endl; testfct4(a);
  std::cout<<"Through fct1 call (ca):"<<std::endl; testfct1(ca);
  std::cout<<"Through fct2 call (ca):"<<std::endl; testfct2(ca);
  std::cout<<"Through fct3 call (ca):"<<std::endl; testfct3(ca);
  std::cout<<"Through fct4 call (ca):"<<std::endl; testfct4(ca);
  NC::Span<double> nonc_span(v);
  std::cout<<"Through fct1 call (nonc_span):"<<std::endl; testfct1(nonc_span);
  std::cout<<"Through fct2 call (nonc_span):"<<std::endl; testfct2(nonc_span);
  std::cout<<"Through fct3 call (nonc_span):"<<std::endl; testfct3(nonc_span);
  std::cout<<"Through fct4 call (nonc_span):"<<std::endl; testfct4(nonc_span);

  //only fct3 and fct4 can be called with const input:
  const std::vector<double> c_v = { 1.0, 1.0, 2.0, 3.0, 5.0, 8.0, 11.0 };
  const std::array<double,4> c_a = { 0.1, 0.2, 0.3, 0.4 };
  const double c_ca[3] = { 0.01, 0.02, 0.03 };
  NC::Span<const double> c_span(v);
  std::cout<<"Through fct3 call (c_v):"<<std::endl; testfct3(c_v);
  std::cout<<"Through fct4 call (c_v):"<<std::endl; testfct4(c_v);
  std::cout<<"Through fct3 call (c_a):"<<std::endl; testfct3(c_a);
  std::cout<<"Through fct4 call (c_a):"<<std::endl; testfct4(c_a);
  std::cout<<"Through fct3 call (c_ca):"<<std::endl; testfct3(c_ca);
  std::cout<<"Through fct4 call (c_ca):"<<std::endl; testfct4(c_ca);
  std::cout<<"Through fct3 call (c_span):"<<std::endl; testfct3(c_span);
  std::cout<<"Through fct4 call (c_span):"<<std::endl; testfct4(c_span);

  {
    std::vector<double> vvv;
    NC::Span<double> ssss;
    ssss = vvv;
  }

  return 0;
}

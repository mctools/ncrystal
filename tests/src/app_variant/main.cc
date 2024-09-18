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

#include "NCrystal/NCVariant.hh"
#include <iostream>
#include <vector>
namespace NC=NCrystal;

namespace {
  static long long nalive = 0;

  class MyPt {
    std::vector<double> m_v;
  public:
    MyPt( double a, double b ) { ++nalive; m_v = {a,b}; }

    MyPt( const MyPt& o ) noexcept { ++nalive; m_v = o.m_v; }

    MyPt( MyPt&& o ) noexcept { ++nalive; m_v = std::move(o.m_v); }
    ~MyPt() { --nalive; }
    void print(std::ostream& os) const {
      nc_assert_always(m_v.empty()||m_v.size()==2);
      if (m_v.empty())
        os << " (moved-from)";
      else
        os<<"("<<m_v.at(0)<<", "<<m_v.at(1)<<")";
    }
  };

  std::ostream& operator<<(std::ostream& os, const MyPt& a) { return a.print(os), os; }
}

void do_tests() {
  NC::Variant<double,MyPt> a;
  nc_assert_always(a.empty());
  std::cout<<"default init a"<<std::endl;
  std::cout<<"a has MyPt: "<<a.has_value<MyPt>()<<std::endl;
  std::cout<<"a has double: "<<a.has_value<double>()<<std::endl;
  std::cout<<"a assigned 17.123"<<std::endl;
  a = 17.123;
  std::cout<<"a has MyPt: "<<a.has_value<MyPt>()<<std::endl;
  std::cout<<"a has double: "<<a.has_value<double>()<<std::endl;
  std::cout<<"a get double: "<< a.get<double>()<<std::endl;
  std::cout<<"a emplace MyPt{1.0,5.5}"<<std::endl;
  a.emplace<MyPt>(1.0,5.5);
  nc_assert_always(!a.empty());
  std::cout<<"a has MyPt: "<<a.has_value<MyPt>()<<std::endl;
  std::cout<<"a has double: "<<a.has_value<double>()<<std::endl;
  std::cout<<"a.get<MyPt>() = "<<a.get<MyPt>()<<std::endl;

  std::cout<<"b = move from a"<<std::endl;
  auto b = std::move(a);
  nc_assert_always(a.empty());
  nc_assert_always(!b.empty());

  std::cout<<"a has MyPt: "<<a.has_value<MyPt>()<<std::endl;
  std::cout<<"a has double: "<<a.has_value<double>()<<std::endl;
  std::cout<<"b has MyPt: "<<b.has_value<MyPt>()<<std::endl;
  std::cout<<"b has double: "<<b.has_value<double>()<<std::endl;

  std::cout<<"b.get<MyPt>() = "<<b.get<MyPt>()<<std::endl;

  std::cout<<"c = copy of b"<<std::endl;
  auto c = b;
  nc_assert_always(a.empty());
  nc_assert_always(!b.empty());
  nc_assert_always(!c.empty());

  std::cout<<"a has MyPt: "<<a.has_value<MyPt>()<<std::endl;
  std::cout<<"a has double: "<<a.has_value<double>()<<std::endl;
  std::cout<<"b has MyPt: "<<b.has_value<MyPt>()<<std::endl;
  std::cout<<"b has double: "<<b.has_value<double>()<<std::endl;
  std::cout<<"c has MyPt: "<<c.has_value<MyPt>()<<std::endl;
  std::cout<<"c has double: "<<c.has_value<double>()<<std::endl;

  std::cout<<"b.get<MyPt>() = "<<b.get<MyPt>()<<std::endl;
  std::cout<<"c.get<MyPt>() = "<<c.get<MyPt>()<<std::endl;

  std::cout<<"d = move from c.get<MyPt>()"<<std::endl;
  MyPt d(std::move(c.get<MyPt>()));
  std::cout<<"c has MyPt: "<<c.has_value<MyPt>()<<std::endl;
  std::cout<<"c has double: "<<c.has_value<double>()<<std::endl;
  std::cout<<"c.get<MyPt>() = "<<c.get<MyPt>()<<std::endl;
  std::cout<<"d: "<<d<<std::endl;

  std::cout<<"c emplace double 12.34"<<std::endl;
  c.emplace<double>(12.34);
  std::cout<<"c has MyPt: "<<c.has_value<MyPt>()<<std::endl;
  std::cout<<"c has double: "<<c.has_value<double>()<<std::endl;
  std::cout<<"c.get<double>() = "<<c.get<double>()<<std::endl;

  std::cout<<"c = double 5.7"<<std::endl;
  c = 5.7;
  std::cout<<"c has MyPt: "<<c.has_value<MyPt>()<<std::endl;
  std::cout<<"c has double: "<<c.has_value<double>()<<std::endl;
  std::cout<<"c.get<double>() = "<<c.get<double>()<<std::endl;

  std::cout<<"c = move-from-(aa=MyPt(1.1,2.2))"<<std::endl;
  MyPt aa(1.1,2.2);
  c = std::move(aa);
  std::cout<<"c has MyPt: "<<c.has_value<MyPt>()<<std::endl;
  std::cout<<"c has double: "<<c.has_value<double>()<<std::endl;
  std::cout<<"c.get<MyPy>() = "<<c.get<MyPt>()<<std::endl;
  std::cout<<"aa = "<<aa<<std::endl;
}

int main() {
  do_tests();
  std::cout<<"MyObjects #constructors - #destructors: "<<nalive<<std::endl;
  nc_assert_always(nalive==0);
  return 0;
}

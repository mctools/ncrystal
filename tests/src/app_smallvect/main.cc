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

#include "NCrystal/core/NCSmallVector.hh"
#include "TestLib_fpe/FPE.hh"

#include <iostream>
#include <vector>

namespace {
  template<class TContainer>
  void print(const TContainer& c,const char * prefix = nullptr) {
    if (prefix)
      std::cout<<prefix;
    std::cout<<"[ capacity="<<c.capacity()<<" : ";
    for (auto& e : c)
      std::cout<<e<<", ";
    std::cout<<" ]"<<std::endl;
  }
}


#include <iostream>
namespace NC = NCrystal;

namespace {
  static long long nalive = 0;
  struct MyDblMoveOnly {
    double m_a;
    MyDblMoveOnly(double a) : m_a(a) { ++nalive; if ( a == 0.12324234345345 ) NCRYSTAL_THROW(BadInput,"argh"); }
    MyDblMoveOnly(MyDblMoveOnly&& o) noexcept : m_a(o.m_a) { ++nalive; o.m_a = -998; }// if ( m_a == 0.12324234345345 ) NCRYSTAL_THROW(BadInput,"argh"); }
    ~MyDblMoveOnly() { --nalive; }
  };

  struct MyDbl {
    double m_a;
    MyDbl() noexcept : m_a(-1.111111) { ++nalive; }
    MyDbl(double a) : m_a(a) { ++nalive; if ( a == 0.12324234345345 ) NCRYSTAL_THROW(BadInput,"argh"); }
    MyDbl(MyDbl&& o) noexcept : m_a(o.m_a) { ++nalive; o.m_a = -999.0; }// if ( m_a == 0.12324234345345 ) NCRYSTAL_THROW(BadInput,"argh"); }
    MyDbl(const MyDbl& o) : m_a(o.m_a) { ++nalive; if ( m_a == 0.12324234345345 ) NCRYSTAL_THROW(BadInput,"argh"); }
    ~MyDbl() { --nalive; }
  };

  class MyPt {
    std::vector<double> m_v;
  public:
    MyPt( double a, double b ) { ++nalive; m_v = {a,b}; }
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


  class MyThrowingPt {
    //Throws on copy if first value of other object is 123.0 or 234.0;
    std::vector<double> m_v;
  public:
    MyThrowingPt( double a, double b )
    { ++nalive;
      m_v = {a,b};
      if ( a == 123.0 )
        throw std::runtime_error("MyThrowingPt(a,b) constructor");
    }
    MyThrowingPt( MyThrowingPt&& o ) noexcept
    {
      //must be noexcept or we won't compile
      ++nalive;
      m_v = std::move(o.m_v);
    }
    MyThrowingPt( const MyThrowingPt& o )
    {
      bool do_throw = (!o.m_v.empty()&&o.m_v.at(0)==234.0);
      m_v = o.m_v;
      if (do_throw)
        throw std::runtime_error("MyThrowingPt copy-constructor");
      ++nalive;//only increment if we don't throw (only full-constructed objects are deleted)
    }
    ~MyThrowingPt() { --nalive; }
    void print(std::ostream& os) const {
      nc_assert_always(m_v.empty()||m_v.size()==2);
      if (m_v.empty())
        os << " (moved-from)";
      else
        os<<"("<<m_v.at(0)<<", "<<m_v.at(1)<<")";
    }
  };

  std::ostream& operator<<(std::ostream& os, const MyDblMoveOnly& a) { return os << a.m_a; }
  std::ostream& operator<<(std::ostream& os, const MyDbl& a) { return os << a.m_a; }
  std::ostream& operator<<(std::ostream& os, const MyPt& a) { return a.print(os), os; }
  std::ostream& operator<<(std::ostream& os, const MyThrowingPt& a) { return a.print(os), os; }

  template <std::size_t alignment>
  inline bool is_aligned(void * ptr) noexcept {
    std::size_t max = 1u;
    return std::align(alignment, 1u, ptr, max);
  }

}

void do_tests() {

  {
    NC::SmallVector<MyDblMoveOnly,3> a;
    print(a);
    a.push_back(1.0);
    print(a);
    a.push_back(2.0);
    print(a);
    std::cout<<"shrink_to_fit"<<std::endl;
    a.shrink_to_fit();
    print(a);
    a.push_back(3.0);
    print(a);
    a.push_back(4.0);
    print(a);
    a.push_back(5.0);
    print(a);
    std::cout<<"shrink_to_fit"<<std::endl;
    a.shrink_to_fit();
    print(a);
    std::cout<<"reserve_hint(17)"<<std::endl;
    a.reserve_hint(17);
    print(a);
    a.push_back(6.0);
    print(a);
    std::cout<<"shrink_to_fit"<<std::endl;
    a.shrink_to_fit();
    print(a);
    std::cout<<"shrink_to_fit"<<std::endl;
    a.shrink_to_fit();
    print(a);
    std::cout<<"clear"<<std::endl;
    a.clear();
    print(a);
    a.push_back(1.0);
    print(a);
    a.push_back(2.0);
    print(a);
    a.push_back(3.0);
    print(a);
    a.push_back(4.0);
    print(a);
    a.push_back(5.0);
    print(a);
    a.push_back(6.0);
    print(a);
    a.emplace_back(7.0);
    print(a);
    a.push_back(8.0);
    print(a);
    a.push_back(9.0);
    print(a);
    a.push_back(10.0);
    print(a);
    //  a = { 1.0, 2.0, 3.0, 4.0 };
    print(a);
    std::cout<<"----- b"<<std::endl;
    NC::SmallVector<MyDbl,2> b;
    print(b);
    b = { 1.0, 2.0, 3.0, 4.0 };
    print(b);
    b = { 1.0, 2.0 };
    print(b);

    std::cout<<"----- c"<<std::endl;
    std::vector<double> v = { 10.0, 20.0, 3.0, 4.0, 5.0 };

    NC::SmallVector<MyDblMoveOnly,1> c( NC::SVAllowCopy,v.begin(),v.end());
    print(c,"c");
    print(a,"a");

    std::cout<<"c.setByMove(..a..)"<<std::endl;
    c.setByMove(a.begin(),a.end());
    print(c,"c");
    print(a,"a");
    NC::SmallVector<MyDbl,2> d = std::move(b);
    d = { 8.,7.,6.,5. };
    print(d,"d");
    std::cout<<"b.setByCopy(..d..)"<<std::endl;
    b.setByCopy(d.begin(),d.end());
    print(d,"d");
    print(b,"b");

    struct alignas(512) MyBigAlign {
      std::vector<double> v;
    };

    NC::SmallVector<MyBigAlign,3> bigalign;
    nc_assert_always(bigalign.isLocalStorage());
    bigalign.emplace_back(MyBigAlign{{1.0,3.0,4.0}});
    nc_assert_always(is_aligned<512>(bigalign.begin()));
    nc_assert_always(bigalign.isLocalStorage());
    bigalign.emplace_back(MyBigAlign{{1.0,3.0,4.0}});
    bigalign.emplace_back(MyBigAlign{{1.0,3.0,4.0}});
    nc_assert_always(bigalign.isLocalStorage());
    bigalign.emplace_back(MyBigAlign{{1.0,3.0,4.0}});
    nc_assert_always(!bigalign.isLocalStorage());
    bigalign.emplace_back(MyBigAlign{{1.0,3.0,4.0}});
    nc_assert_always(is_aligned<512>(bigalign.begin()));
    nc_assert_always(!bigalign.isLocalStorage());
    nc_assert_always(bigalign.isFastAccess());

    nc_assert_always( (NC::SmallVector<MyBigAlign,3>::isFastAccess()==true) );
    nc_assert_always( (NC::SmallVector<MyDbl,3>::isFastAccess()==true) );

    nc_assert_always( (NC::SmallVector<MyBigAlign,3,NC::SVMode::LOWFOOTPRINT>::isFastAccess()==true) );
    nc_assert_always( (NC::SmallVector<float,3,NC::SVMode::LOWFOOTPRINT>::isFastAccess()==false) );
    std::cout<<"sizeof(MyBigAlign) = "<<sizeof(MyBigAlign)<<std::endl;
    std::cout<<"alignof(MyBigAlign) = "<<alignof(MyBigAlign)<<std::endl;

    NC::SmallVector<double,7> aa = { 2.0, 6.0, 8.0, 9.0, 2.1, 5., 4., 2. };
    NC::SmallVector<double,3> aaa( NC::SVAllowCopy, aa.begin(), aa.end() );
    print(aa,"aa");
    print(aaa,"aaa");
    std::cout<<"sort aa"<<std::endl;
    std::sort(aa.begin(),aa.end());
    print(aa,"aa");
    print(aaa,"aaa");
    std::cout<<"sort aaa[0:-1]"<<std::endl;
    std::sort(aaa.begin(),std::prev(aaa.end()));
    print(aa,"aa");
    print(aaa,"aaa");

    NC::SmallVector<double,7,NC::SVMode::LOWFOOTPRINT> aaaa = { 2.0, 6.0, 8.0, 9.0, 2.1, 5., 4., 2. };
    print(aaaa,"aaaa");
    std::cout<<"sort aaaa"<<std::endl;
    std::sort(aaaa.begin(),aaaa.end());
    print(aaaa,"aaaa");

    NC::SmallVector<MyPt,2> vpt;
    NC::SmallVector<MyPt,2> vpt2;
    print(vpt,"vpt");
    nc_assert_always(vpt.capacity()==2);
    std::cout<<"emplace_back(4,5)"<<std::endl; vpt.emplace_back(4,5);
    print(vpt,"vpt");
    std::cout<<"reserve_hint(17)"<<std::endl; vpt.reserve_hint(17);
    print(vpt,"vpt");
    nc_assert_always(vpt.capacity()==2);

    std::cout<<"move-assign to vpt2"<<std::endl;
    vpt2 = std::move(vpt);
    print(vpt,"vpt");
    print(vpt2,"vpt2");
    nc_assert_always(vpt.capacity()==2);

    std::cout<<"push_back(40,50)"<<std::endl; vpt.push_back(MyPt(40,50));
    print(vpt,"vpt");
    nc_assert_always(vpt.capacity()==2);

    std::cout<<"push_back(41,51)"<<std::endl; vpt.push_back(MyPt(41,51));
    print(vpt,"vpt");
    nc_assert_always(vpt.capacity()==2);

    std::cout<<"push_back(42,52)"<<std::endl; vpt.push_back(MyPt(42,52));
    print(vpt,"vpt");
    nc_assert_always(vpt.capacity()==4);

    std::cout<<"push_back(1,2)"<<std::endl; vpt.push_back(MyPt(1,2));
    print(vpt,"vpt");
    std::cout<<"reserve_hint(17)"<<std::endl; vpt.reserve_hint(17);
    print(vpt,"vpt");
    nc_assert_always(vpt.capacity()==17);

    std::cout<<"move-assign to vpt2"<<std::endl;
    vpt2 = std::move(vpt);
    print(vpt,"vpt");
    print(vpt2,"vpt2");
    nc_assert_always(vpt.capacity()==2);
    nc_assert_always(vpt2.capacity()==17);

    std::cout<<"vpt2.emplace_back 13 times"<<std::endl;
    for ( unsigned i = 0; i<13; ++i )
      vpt2.emplace_back(100+i,200+i);
    print(vpt2,"vpt2");
    std::cout<<"vpt2.emplace_back"<<std::endl;
    vpt2.emplace_back(0.123,0.456);
    print(vpt2,"vpt2");
    std::cout<<"vpt2.shrink_to_fit"<<std::endl;
    vpt2.shrink_to_fit();
    print(vpt,"vpt");
    print(vpt2,"vpt2");
    std::cout<<"vpt.swap(vpt2)"<<std::endl;
    vpt.swap(vpt2);
    print(vpt,"vpt");
    print(vpt2,"vpt2");

    while (!vpt.empty()) {
      std::cout<<"vpt.pop_back()"<<std::endl;
      vpt.pop_back();
      print(vpt,"vpt");
    }
  }
  nc_assert_always(nalive==0);
  {
    std::cout<<"---------------------------------------"<<std::endl;
    NC::SmallVector<MyDbl,3> u = { 1, 2, 3, 4, 5, 6, 7, 8 };
    print(u,"u (SmallVector<MyDbl,NSMALL=3>) ");
    nc_assert_always((std::size_t)nalive==u.size());
    //  vpt.resize(10);
    std::cout<<"  Other(SVAllowCopy, u)"<<std::endl;
    NC::SmallVector<MyDbl,3>( NC::SVAllowCopy, u );
    print(u,"u");
    nc_assert_always((std::size_t)nalive==u.size());
    std::cout<<"  u.resize(8)"<<std::endl;
    u.resize(8);
    print(u,"u");
    nc_assert_always((std::size_t)nalive==u.size());
    std::cout<<"  u.resize(4)"<<std::endl;
    u.resize(4);
    print(u,"u");
    nc_assert_always((std::size_t)nalive==u.size());
    std::cout<<"  u.resize(2)"<<std::endl;
    u.resize(2);
    print(u,"u");
    nc_assert_always((std::size_t)nalive==u.size());
    std::cout<<"  u.resize(3)"<<std::endl;
    u.resize(3);
    print(u,"u");
    std::cout<<" nalive="<<nalive<<std::endl;
    nc_assert_always((std::size_t)nalive==u.size());

    std::cout<<"  u.resize(10)"<<std::endl;
    u.resize(10);
    print(u,"u");

    std::cout<<"  u.resize(0)"<<std::endl;
    u.resize(0);
    print(u,"u");

    std::cout<<"  u = {1,2}"<<std::endl;
    u = {1,2};
    print(u,"u");

    std::cout<<"  u.resize(5)"<<std::endl;
    u.resize(5);
    print(u,"u");
  }
  {
    std::cout<<"---------------------------------------"<<std::endl;
    {  NC::SmallVector<int,2> u{ 2 }; print(u,"SmallVector<int,2>{ 2 } = "); }
    {  NC::SmallVector<int,2> u{ 2, 4 }; print(u,"SmallVector<int,2>{ 2, 4 } = "); }
    {  NC::SmallVector<int,2> u{ 2, 4, 7 }; print(u,"SmallVector<int,2>{ 2, 4, 7 } = "); }
    //{  NC::SmallVector<int,2> u{ NC::SVCountConstruct, 2 }; print(u,"SmallVector<int,2>{ SVCountConstruct, 2 } = "); }
    {  NC::SmallVector<int,2> u{ NC::SVCountConstruct, 2, 4 }; print(u,"SmallVector<int,2>{ SVCountConstruct, 2, 4 } = "); }
    {  NC::SmallVector<int,2> u{ 5 }; print(u,"SmallVector<int,2>{ 5 } = "); }
    {  NC::SmallVector<int,2> u{ 5, 4 }; print(u,"SmallVector<int,2>{ 5, 4 } = "); }
    {  NC::SmallVector<int,2> u{ 5, 4, 7 }; print(u,"SmallVector<int,2>{ 5, 4, 7 } = "); }
    //{  NC::SmallVector<int,2> u{ NC::SVCountConstruct, 5 }; print(u,"SmallVector<int,2>{ SVCountConstruct, 5 } = "); }
    {  NC::SmallVector<int,2> u{ NC::SVCountConstruct, 5, 4 }; print(u,"SmallVector<int,2>{ SVCountConstruct, 5, 4 } = "); }
    {  NC::SmallVector<MyDbl,2> u{ 2 }; print(u,"SmallVector<MyDbl,2>{ 2 } = "); }
    {  NC::SmallVector<MyDbl,2> u{ 2, 4 }; print(u,"SmallVector<MyDbl,2>{ 2, 4 } = "); }
    {  NC::SmallVector<MyDbl,2> u{ 2, 4, 7 }; print(u,"SmallVector<MyDbl,2>{ 2, 4, 7 } = "); }
    //    {  NC::SmallVector<MyDbl,2> u{ NC::SVCountConstruct, 2 }; print(u,"SmallVector<MyDbl,2>{ SVCountConstruct, 2 } = "); }
    {  NC::SmallVector<MyDbl,2> u{ NC::SVCountConstruct, 2, 4 }; print(u,"SmallVector<MyDbl,2>{ SVCountConstruct, 2, 4 } = "); }
    {  NC::SmallVector<MyDbl,2> u{ 5 }; print(u,"SmallVector<MyDbl,2>{ 5 } = "); }
    {  NC::SmallVector<MyDbl,2> u{ 5, 4 }; print(u,"SmallVector<MyDbl,2>{ 5, 4 } = "); }
    {  NC::SmallVector<MyDbl,2> u{ 5, 4, 7 }; print(u,"SmallVector<MyDbl,2>{ 5, 4, 7 } = "); }
    //{  NC::SmallVector<MyDbl,2> u{ NC::SVCountConstruct, 5 }; print(u,"SmallVector<MyDbl,2>{ SVCountConstruct, 5 } = "); }
    {  NC::SmallVector<MyDbl,2> u{ NC::SVCountConstruct, 5, 4 }; print(u,"SmallVector<MyDbl,2>{ SVCountConstruct, 5, 4 } = "); }
  }
  {
    std::cout<<"---------------------------------------"<<std::endl;
    NC::SmallVector<MyDbl,4> u{1,2,3}, v{4,5,6,7,8,9};
    print(u,"u");
    print(v,"v");
    std::cout<<"v = std::move(u)"<<std::endl;
    v = std::move(u);
    print(u,"u");
    print(v,"v");
  }
  {
    std::cout<<"---------------------------------------"<<std::endl;
    NC::SmallVector<MyDbl,4> u{1,2,3}, v{4,5,6,7,8,9};
    print(u,"u");
    print(v,"v");
    std::cout<<"u = std::move(v)"<<std::endl;
    u = std::move(v);
    print(u,"u");
    print(v,"v");
    std::cout<<"v.push_back(17)"<<std::endl;
    v.push_back(17);
    print(u,"u");
    print(v,"v");
  }

  {
    std::cout<<"---------------------------------------"<<std::endl;

    //two vectors with heap storage:
    NC::SmallVector<MyDbl,1> w{1,2,3}, x{4,5,6,7,8,9};
    print(w,"w");
    print(x,"x");
    auto wdata(w.data()), xdata(x.data());
    std::swap(w,x);//should be a simple pointer swap.
    nc_assert_always(wdata = x.data());
    nc_assert_always(xdata = w.data());
    print(w,"w");
    print(x,"x");

    //two vectors with local storage:
    NC::SmallVector<MyDbl,10> y{1,2,3}, z{4,5,6,7,8,9};
    print(y,"y");
    print(z,"z");
    std::swap(y,z);
    print(y,"y");
    print(z,"z");

    //One with local, one with heap:
    NC::SmallVector<MyDbl,4> u{1,2,3}, v{4,5,6,7,8,9};
    print(u,"u");
    print(v,"v");
    std::cout<<"u.swap(v)"<<std::endl;
    u.swap(v);
    print(u,"u");
    print(v,"v");
    //Also test std::swap, ADL swap:
    std::cout<<"std::swap(u,v)"<<std::endl;
    std::swap(u,v);
    print(u,"u");
    print(v,"v");
    //TODO: Try all small/large combinations
    [&u,&v](){
      std::cout<<"ADL: swap(u,v)"<<std::endl;
      swap(u,v);
      print(u,"u");
      print(v,"v");
    }();
    [&u,&v](){
      std::cout<<"using std::swap; swap(u,v)"<<std::endl;
      using std::swap;
      swap(u,v);
      print(u,"u");
      print(v,"v");
    }();
  }

  NC::SmallVector<NC::SmallVector<MyDbl,7,NC::SVMode::LOWFOOTPRINT>,4> u;
  u.resize(7);
  u.push_back({2,4,5,6,7});
  nc_assert(u.size()==8);
  u.resize(3);
  nc_assert(u.size()==3);

  //Unit test that a push_back/emplace_back a throwing TValue constructor does
  //not mess up the SmallVector state:

  for ( unsigned i : { 3, 4, 5, 6, 7, 8, 9, 10, 11 } ) {
    std::cout<<"--------------------------------------- TEST THROWING COPY CONSTRUCTOR"<<std::endl;
    NC::SmallVector<MyThrowingPt,4> s;
    while ( s.size() < i )
      s.push_back({s.size() + 1.0,s.size() +1.1});
    print(s,"s");
    MyThrowingPt point_which_throws_on_copy{234.0,52.1};
    std::cout<<"now push_back with "<<point_which_throws_on_copy<<" (should throw)"<<std::endl;
    try {
      s.push_back(point_which_throws_on_copy);
    } catch ( std::runtime_error& e ) {
      std::cout<<"--> got exception: "<< e.what() << std::endl;
    }
    print(s,"s");
  }

  {
    NC::SmallVector_IC<double,2> v_copyable = { 1.1, 2.2 };
    decltype(v_copyable) v_copyable2(v_copyable);
    decltype(v_copyable) v_copyable3;
    v_copyable3 = v_copyable;
  }
}


void do_tests_IC() {
  {
    NC::SmallVector_IC<MyDblMoveOnly,3> a;
    print(a);
    a.push_back(1.0);
    print(a);
    a.push_back(2.0);
    print(a);
    std::cout<<"shrink_to_fit"<<std::endl;
    a.shrink_to_fit();
    print(a);
    a.push_back(3.0);
    print(a);
    a.push_back(4.0);
    print(a);
    a.push_back(5.0);
    print(a);
    std::cout<<"shrink_to_fit"<<std::endl;
    a.shrink_to_fit();
    print(a);
    std::cout<<"reserve_hint(17)"<<std::endl;
    a.reserve_hint(17);
    print(a);
    a.push_back(6.0);
    print(a);
    std::cout<<"shrink_to_fit"<<std::endl;
    a.shrink_to_fit();
    print(a);
    std::cout<<"shrink_to_fit"<<std::endl;
    a.shrink_to_fit();
    print(a);
    std::cout<<"clear"<<std::endl;
    a.clear();
    print(a);
    a.push_back(1.0);
    print(a);
    a.push_back(2.0);
    print(a);
    a.push_back(3.0);
    print(a);
    a.push_back(4.0);
    print(a);
    a.push_back(5.0);
    print(a);
    a.push_back(6.0);
    print(a);
    a.emplace_back(7.0);
    print(a);
    a.push_back(8.0);
    print(a);
    a.push_back(9.0);
    print(a);
    a.push_back(10.0);
    print(a);
    //  a = { 1.0, 2.0, 3.0, 4.0 };
    print(a);
    std::cout<<"----- b"<<std::endl;
    NC::SmallVector_IC<MyDbl,2> b;
    print(b);
    b = { 1.0, 2.0, 3.0, 4.0 };
    print(b);
    b = { 1.0, 2.0 };
    print(b);

    std::cout<<"----- c"<<std::endl;
    std::vector<double> v = { 10.0, 20.0, 3.0, 4.0, 5.0 };

    NC::SmallVector_IC<MyDblMoveOnly,1> c( NC::SVAllowCopy,v.begin(),v.end());
    print(c,"c");
    print(a,"a");

    std::cout<<"c.setByMove(..a..)"<<std::endl;
    c.setByMove(a.begin(),a.end());
    print(c,"c");
    print(a,"a");
    NC::SmallVector_IC<MyDbl,2> d = std::move(b);
    d = { 8.,7.,6.,5. };
    print(d,"d");
    std::cout<<"b.setByCopy(..d..)"<<std::endl;
    b.setByCopy(d.begin(),d.end());
    print(d,"d");
    print(b,"b");

    struct alignas(512) MyBigAlign {
      std::vector<double> v;
    };

    NC::SmallVector_IC<MyBigAlign,3> bigalign;
    nc_assert_always(bigalign.isLocalStorage());
    bigalign.emplace_back(MyBigAlign{{1.0,3.0,4.0}});
    nc_assert_always(is_aligned<512>(bigalign.begin()));
    nc_assert_always(bigalign.isLocalStorage());
    bigalign.emplace_back(MyBigAlign{{1.0,3.0,4.0}});
    bigalign.emplace_back(MyBigAlign{{1.0,3.0,4.0}});
    nc_assert_always(bigalign.isLocalStorage());
    bigalign.emplace_back(MyBigAlign{{1.0,3.0,4.0}});
    nc_assert_always(!bigalign.isLocalStorage());
    bigalign.emplace_back(MyBigAlign{{1.0,3.0,4.0}});
    nc_assert_always(is_aligned<512>(bigalign.begin()));
    nc_assert_always(!bigalign.isLocalStorage());
    nc_assert_always(bigalign.isFastAccess());

    nc_assert_always( (NC::SmallVector_IC<MyBigAlign,3>::isFastAccess()==true) );
    nc_assert_always( (NC::SmallVector_IC<MyDbl,3>::isFastAccess()==true) );

    nc_assert_always( (NC::SmallVector_IC<MyBigAlign,3,NC::SVMode::LOWFOOTPRINT>::isFastAccess()==true) );
    nc_assert_always( (NC::SmallVector_IC<float,3,NC::SVMode::LOWFOOTPRINT>::isFastAccess()==false) );
    std::cout<<"sizeof(MyBigAlign) = "<<sizeof(MyBigAlign)<<std::endl;
    std::cout<<"alignof(MyBigAlign) = "<<alignof(MyBigAlign)<<std::endl;

    NC::SmallVector_IC<double,7> aa = { 2.0, 6.0, 8.0, 9.0, 2.1, 5., 4., 2. };
    NC::SmallVector_IC<double,3> aaa( NC::SVAllowCopy, aa.begin(), aa.end() );
    print(aa,"aa");
    print(aaa,"aaa");
    std::cout<<"sort aa"<<std::endl;
    std::sort(aa.begin(),aa.end());
    print(aa,"aa");
    print(aaa,"aaa");
    std::cout<<"sort aaa[0:-1]"<<std::endl;
    std::sort(aaa.begin(),std::prev(aaa.end()));
    print(aa,"aa");
    print(aaa,"aaa");

    NC::SmallVector_IC<double,7,NC::SVMode::LOWFOOTPRINT> aaaa = { 2.0, 6.0, 8.0, 9.0, 2.1, 5., 4., 2. };
    print(aaaa,"aaaa");
    std::cout<<"sort aaaa"<<std::endl;
    std::sort(aaaa.begin(),aaaa.end());
    print(aaaa,"aaaa");

    NC::SmallVector_IC<MyPt,2> vpt;
    NC::SmallVector_IC<MyPt,2> vpt2;
    print(vpt,"vpt");
    nc_assert_always(vpt.capacity()==2);
    std::cout<<"emplace_back(4,5)"<<std::endl; vpt.emplace_back(4,5);
    print(vpt,"vpt");
    std::cout<<"reserve_hint(17)"<<std::endl; vpt.reserve_hint(17);
    print(vpt,"vpt");
    nc_assert_always(vpt.capacity()==2);

    std::cout<<"move-assign to vpt2"<<std::endl;
    vpt2 = std::move(vpt);
    print(vpt,"vpt");
    print(vpt2,"vpt2");
    nc_assert_always(vpt.capacity()==2);

    std::cout<<"push_back(40,50)"<<std::endl; vpt.push_back(MyPt(40,50));
    print(vpt,"vpt");
    nc_assert_always(vpt.capacity()==2);

    std::cout<<"push_back(41,51)"<<std::endl; vpt.push_back(MyPt(41,51));
    print(vpt,"vpt");
    nc_assert_always(vpt.capacity()==2);

    std::cout<<"push_back(42,52)"<<std::endl; vpt.push_back(MyPt(42,52));
    print(vpt,"vpt");
    nc_assert_always(vpt.capacity()==4);

    std::cout<<"push_back(1,2)"<<std::endl; vpt.push_back(MyPt(1,2));
    print(vpt,"vpt");
    std::cout<<"reserve_hint(17)"<<std::endl; vpt.reserve_hint(17);
    print(vpt,"vpt");
    nc_assert_always(vpt.capacity()==17);

    std::cout<<"move-assign to vpt2"<<std::endl;
    vpt2 = std::move(vpt);
    print(vpt,"vpt");
    print(vpt2,"vpt2");
    nc_assert_always(vpt.capacity()==2);
    nc_assert_always(vpt2.capacity()==17);

    std::cout<<"vpt2.emplace_back 13 times"<<std::endl;
    for ( unsigned i = 0; i<13; ++i )
      vpt2.emplace_back(100+i,200+i);
    print(vpt2,"vpt2");
    std::cout<<"vpt2.emplace_back"<<std::endl;
    vpt2.emplace_back(0.123,0.456);
    print(vpt2,"vpt2");
    std::cout<<"vpt2.shrink_to_fit"<<std::endl;
    vpt2.shrink_to_fit();
    print(vpt,"vpt");
    print(vpt2,"vpt2");
    std::cout<<"vpt.swap(vpt2)"<<std::endl;
    vpt.swap(vpt2);
    print(vpt,"vpt");
    print(vpt2,"vpt2");

    while (!vpt.empty()) {
      std::cout<<"vpt.pop_back()"<<std::endl;
      vpt.pop_back();
      print(vpt,"vpt");
    }
  }
  nc_assert_always(nalive==0);
  {
    std::cout<<"---------------------------------------"<<std::endl;
    NC::SmallVector_IC<MyDbl,3> u = { 1, 2, 3, 4, 5, 6, 7, 8 };
    print(u,"u (SmallVector_IC<MyDbl,NSMALL=3>) ");
    nc_assert_always((std::size_t)nalive==u.size());
    //  vpt.resize(10);
    std::cout<<"  Other(SVAllowCopy, u)"<<std::endl;
    NC::SmallVector_IC<MyDbl,3>( NC::SVAllowCopy, u );
    print(u,"u");
    nc_assert_always((std::size_t)nalive==u.size());
    std::cout<<"  u.resize(8)"<<std::endl;
    u.resize(8);
    print(u,"u");
    nc_assert_always((std::size_t)nalive==u.size());
    std::cout<<"  u.resize(4)"<<std::endl;
    u.resize(4);
    print(u,"u");
    nc_assert_always((std::size_t)nalive==u.size());
    std::cout<<"  u.resize(2)"<<std::endl;
    u.resize(2);
    print(u,"u");
    nc_assert_always((std::size_t)nalive==u.size());
    std::cout<<"  u.resize(3)"<<std::endl;
    u.resize(3);
    print(u,"u");
    std::cout<<" nalive="<<nalive<<std::endl;
    nc_assert_always((std::size_t)nalive==u.size());

    std::cout<<"  u.resize(10)"<<std::endl;
    u.resize(10);
    print(u,"u");

    std::cout<<"  u.resize(0)"<<std::endl;
    u.resize(0);
    print(u,"u");

    std::cout<<"  u = {1,2}"<<std::endl;
    u = {1,2};
    print(u,"u");

    std::cout<<"  u.resize(5)"<<std::endl;
    u.resize(5);
    print(u,"u");
  }
  {
    std::cout<<"---------------------------------------"<<std::endl;
    {  NC::SmallVector_IC<int,2> u{ 2 }; print(u,"SmallVector_IC<int,2>{ 2 } = "); }
    {  NC::SmallVector_IC<int,2> u{ 2, 4 }; print(u,"SmallVector_IC<int,2>{ 2, 4 } = "); }
    {  NC::SmallVector_IC<int,2> u{ 2, 4, 7 }; print(u,"SmallVector_IC<int,2>{ 2, 4, 7 } = "); }
    // {  NC::SmallVector_IC<int,2> u{ NC::SVCountConstruct, 2 }; print(u,"SmallVector_IC<int,2>{ SVCountConstruct, 2 } = "); }
    {  NC::SmallVector_IC<int,2> u{ NC::SVCountConstruct, 2, 4 }; print(u,"SmallVector_IC<int,2>{ SVCountConstruct, 2, 4 } = "); }
    {  NC::SmallVector_IC<int,2> u{ 5 }; print(u,"SmallVector_IC<int,2>{ 5 } = "); }
    {  NC::SmallVector_IC<int,2> u{ 5, 4 }; print(u,"SmallVector_IC<int,2>{ 5, 4 } = "); }
    {  NC::SmallVector_IC<int,2> u{ 5, 4, 7 }; print(u,"SmallVector_IC<int,2>{ 5, 4, 7 } = "); }
    // {  NC::SmallVector_IC<int,2> u{ NC::SVCountConstruct, 5 }; print(u,"SmallVector_IC<int,2>{ SVCountConstruct, 5 } = "); }
    {  NC::SmallVector_IC<int,2> u{ NC::SVCountConstruct, 5, 4 }; print(u,"SmallVector_IC<int,2>{ SVCountConstruct, 5, 4 } = "); }
    {  NC::SmallVector_IC<MyDbl,2> u{ 2 }; print(u,"SmallVector_IC<MyDbl,2>{ 2 } = "); }
    {  NC::SmallVector_IC<MyDbl,2> u{ 2, 4 }; print(u,"SmallVector_IC<MyDbl,2>{ 2, 4 } = "); }
    {  NC::SmallVector_IC<MyDbl,2> u{ 2, 4, 7 }; print(u,"SmallVector_IC<MyDbl,2>{ 2, 4, 7 } = "); }
    // {  NC::SmallVector_IC<MyDbl,2> u{ NC::SVCountConstruct, 2 }; print(u,"SmallVector_IC<MyDbl,2>{ SVCountConstruct, 2 } = "); }
    {  NC::SmallVector_IC<MyDbl,2> u{ NC::SVCountConstruct, 2, 4 }; print(u,"SmallVector_IC<MyDbl,2>{ SVCountConstruct, 2, 4 } = "); }
    {  NC::SmallVector_IC<MyDbl,2> u{ 5 }; print(u,"SmallVector_IC<MyDbl,2>{ 5 } = "); }
    {  NC::SmallVector_IC<MyDbl,2> u{ 5, 4 }; print(u,"SmallVector_IC<MyDbl,2>{ 5, 4 } = "); }
    {  NC::SmallVector_IC<MyDbl,2> u{ 5, 4, 7 }; print(u,"SmallVector_IC<MyDbl,2>{ 5, 4, 7 } = "); }
    // {  NC::SmallVector_IC<MyDbl,2> u{ NC::SVCountConstruct, 5 }; print(u,"SmallVector_IC<MyDbl,2>{ SVCountConstruct, 5 } = "); }
    {  NC::SmallVector_IC<MyDbl,2> u{ NC::SVCountConstruct, 5, 4 }; print(u,"SmallVector_IC<MyDbl,2>{ SVCountConstruct, 5, 4 } = "); }
  }
  {
    std::cout<<"---------------------------------------"<<std::endl;
    NC::SmallVector_IC<MyDbl,4> u{1,2,3}, v{4,5,6,7,8,9};
    print(u,"u");
    print(v,"v");
    std::cout<<"v = std::move(u)"<<std::endl;
    v = std::move(u);
    print(u,"u");
    print(v,"v");
  }
  {
    std::cout<<"---------------------------------------"<<std::endl;
    NC::SmallVector_IC<MyDbl,4> u{1,2,3}, v{4,5,6,7,8,9};
    print(u,"u");
    print(v,"v");
    std::cout<<"u = std::move(v)"<<std::endl;
    u = std::move(v);
    print(u,"u");
    print(v,"v");
    std::cout<<"v.push_back(17)"<<std::endl;
    v.push_back(17);
    print(u,"u");
    print(v,"v");
  }

  {
    std::cout<<"---------------------------------------"<<std::endl;

    //two vectors with heap storage:
    NC::SmallVector_IC<MyDbl,1> w{1,2,3}, x{4,5,6,7,8,9};
    print(w,"w");
    print(x,"x");
    auto wdata(w.data()), xdata(x.data());
    std::swap(w,x);//should be a simple pointer swap.
    nc_assert_always(wdata = x.data());
    nc_assert_always(xdata = w.data());
    print(w,"w");
    print(x,"x");

    //two vectors with local storage:
    NC::SmallVector_IC<MyDbl,10> y{1,2,3}, z{4,5,6,7,8,9};
    print(y,"y");
    print(z,"z");
    std::swap(y,z);
    print(y,"y");
    print(z,"z");

    //One with local, one with heap:
    NC::SmallVector_IC<MyDbl,4> u{1,2,3}, v{4,5,6,7,8,9};
    print(u,"u");
    print(v,"v");
    std::cout<<"u.swap(v)"<<std::endl;
    u.swap(v);
    print(u,"u");
    print(v,"v");
    //Also test std::swap, ADL swap:
    std::cout<<"std::swap(u,v)"<<std::endl;
    std::swap(u,v);
    print(u,"u");
    print(v,"v");
    //TODO: Try all small/large combinations
    [&u,&v](){
      std::cout<<"ADL: swap(u,v)"<<std::endl;
      swap(u,v);
      print(u,"u");
      print(v,"v");
    }();
    [&u,&v](){
      std::cout<<"using std::swap; swap(u,v)"<<std::endl;
      using std::swap;
      swap(u,v);
      print(u,"u");
      print(v,"v");
    }();
  }

  NC::SmallVector_IC<NC::SmallVector_IC<MyDbl,7,NC::SVMode::LOWFOOTPRINT>,4> u;
  u.resize(7);
  u.push_back({2,4,5,6,7});
  nc_assert(u.size()==8);
  u.resize(3);
  nc_assert(u.size()==3);

  //Unit test that a push_back/emplace_back a throwing TValue constructor does
  //not mess up the SmallVector_IC state:

  for ( unsigned i : { 3, 4, 5, 6, 7, 8, 9, 10, 11 } ) {
    std::cout<<"--------------------------------------- TEST THROWING COPY CONSTRUCTOR"<<std::endl;
    NC::SmallVector_IC<MyThrowingPt,4> s;
    while ( s.size() < i )
      s.push_back({s.size() + 1.0,s.size() +1.1});
    print(s,"s");
    MyThrowingPt point_which_throws_on_copy{234.0,52.1};
    std::cout<<"now push_back with "<<point_which_throws_on_copy<<" (should throw)"<<std::endl;
    try {
      s.push_back(point_which_throws_on_copy);
    } catch ( std::runtime_error& e ) {
      std::cout<<"--> got exception: "<< e.what() << std::endl;
    }
    print(s,"s");
  }

  {
    NC::SmallVector_IC<double,2> v_copyable = { 1.1, 2.2 };
    decltype(v_copyable) v_copyable2(v_copyable);
    decltype(v_copyable) v_copyable3;
    v_copyable3 = v_copyable;
  }
}


int main() {
  NCTests::catch_fpe();
  do_tests();
  do_tests_IC();
  std::cout<<"MyObjects #constructors - #destructors: "<<nalive<<std::endl;
  nc_assert_always(nalive==0);
  return 0;
}

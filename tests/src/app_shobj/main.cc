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

#include "NCrystal/core/NCDefs.hh"

#include <memory>
#include <type_traits>
#include <iostream>

namespace NC = NCrystal;

namespace {


  struct Dummy {
    int a = 17;
    Dummy(int aa) : a(aa) {}
  };

  NC::shared_obj<const Dummy> createDummy() { return NC::makeSO<Dummy>(123); }

  struct DummyBase {
    virtual void print() const = 0;
    virtual ~DummyBase() = default;
  };
  struct DummyDerived : DummyBase {
    void print() const override { std::cout<<"I am DummyDerived"<<std::endl; }
  };

  void fct1(Dummy& d){
    std::cout<<"fct1: "<<d.a<<std::endl;
  }
  void fct1b(const Dummy& d){
    std::cout<<"fct1b: "<<d.a<<std::endl;
  }
  void fct2(Dummy* d){
    std::cout<<"fct2: "<<d->a<<std::endl;
  }
  void fct2b(const Dummy* d){
    std::cout<<"fct2b: "<<d->a<<std::endl;
  }
  void fct3(std::shared_ptr<Dummy> d) {
    std::cout<<"fct3: "<<d->a<<std::endl;
  }
  void fct4(std::shared_ptr<const Dummy> d) {
    std::cout<<"fct4: "<<d->a<<std::endl;
  }
  void fct5(NC::shared_obj<const Dummy> d) {
    std::cout<<"fct5: "<<d->a<<std::endl;
  }
  void fct5b(NC::shared_obj<Dummy> d) {
    std::cout<<"fct5b: "<<d->a<<std::endl;
  }

  NC::shared_obj<Dummy> fct6() { return std::make_shared<Dummy>(117); }
  std::shared_ptr<Dummy> fct7() { return NC::makeSO<Dummy>(118); }
}

int main()
{
  NC::shared_obj<Dummy> d(new Dummy{12});
  NC::shared_obj<const Dummy> dc(d);
  fct1(*d);
  fct1(d);
  fct1b(*d);
  fct1b(d);
  fct2(d.get());
  fct2((Dummy*)d);//requires explicit cast for safety
  fct2b(d.get());
  fct2b((const Dummy*)d);//requires explicit cast for safety
  fct3(d);
  fct4(d);
  fct5(d);
  fct5b(d);
  auto d2 = NC::makeSO<Dummy>(14);
  fct1(*d2);
  fct1(d2);
  fct1b(*d2);
  fct1b(d2);
  fct2(d2.get());
  fct2((Dummy*)d2);//requires explicit cast for safety
  fct2b(d2.get());
  fct2b((const Dummy*)d2);//requires explicit cast for safety
  fct3(d2);
  fct4(d2);
  fct5(d2);
  fct5b(d2);

  auto d3 = NC::makeSO<const Dummy>(15);
  //  fct1(*d3);
  //fct1(d3);
  fct1b(*d3);
  fct1b(d3);
  //  fct2(d3.get());
  // fct2(d3);
  // fct3(d3);
  fct4(d3);
  fct5(d3);
  //fct5b(d3);

  auto sp1 = std::make_shared<Dummy>(20);
  fct5(sp1);
  fct5b(sp1);
  auto sp2 = std::make_shared<const Dummy>(21);
  fct5(sp2);

  fct5(fct6());
  fct5(fct7());
  //todo: verify exceptions:
  //  NC::shared_obj<const Dummy>(nullptr);

  NC::shared_obj<const DummyBase> dbase = NC::makeSO<DummyDerived>();
  NC::shared_obj<const DummyBase> dbase2{ NC::makeSO<const DummyDerived>()};
  NC::shared_obj<DummyBase> dbase3 = NC::makeSO<DummyDerived>();
  NC::shared_obj<DummyBase> dbase4{ NC::makeSO<DummyDerived>()};
  //should not compile: NC::shared_obj<const DummyBase> dbase3{ NC::makeSO<const Dummy>()};
  dbase->print();
  dbase2->print();
  dbase3->print();
  dbase4->print();

  NC::optional_shared_obj<const DummyBase> opt = NC::makeSO<DummyDerived>();

  NC::shared_obj<const DummyBase> dbaseb = std::make_shared<DummyDerived>();
  NC::shared_obj<const DummyBase> dbaseb2{ std::make_shared<const DummyDerived>()};
  NC::shared_obj<DummyBase> dbaseb3 = std::make_shared<DummyDerived>();
  NC::shared_obj<DummyBase> dbaseb4{ std::make_shared<DummyDerived>()};
  //will throw exception: NC::shared_obj<DummyBase> dbaseb5{ std::shared_ptr<DummyDerived>()};

  dbaseb->print();
  dbaseb2->print();
  dbaseb3->print();
  dbaseb4->print();

  struct A : public NC::EnableSharedFromThis<A> {
    int dummy = 17;
  };

  //First with OK object (kept in shared object):
  {
    auto spA = NC::makeSO<A>();
    (void)spA->dummy;//for clang
    A* a = spA.get();
    auto spA2 = a->shared_obj_from_this();
    const A* ac = a;
    auto spA3c = ac->shared_obj_from_this();
  }
  //Now with *not* OK object (not kept in shared object):
  {
    auto spA = std::make_unique<A>();
    A* a = spA.get();
    bool gotexpectederror(false);
    try {
#if nc_cplusplus >= 201703L
      auto spA2 = a->shared_obj_from_this();
#else
      //Pre-C++17 this was actually undefined behaviour!! So to avoid a
      //spuriously failing unit test, we emit the same error here:
      NCRYSTAL_THROW( LogicError, "shared_obj_from_this() does not work since this"
                      " instance is not managed by std::shared_ptr / NCrystal::shared_obj" );

#endif
    } catch ( NCrystal::Error::LogicError& e ) {
      std::cout<<"Got expected error: "<<e.what()<<std::endl;
      gotexpectederror = true;
    }
    nc_assert_always(gotexpectederror);
  }

  {
    auto dummy1 = createDummy();
    //This unsafe way no longer compiles:
    //  const Dummy* dummy_ptr_dangerous = createDummy();

    //UPDATE this now now longer compilers:
    //This unsafe way unfortunately still compiles (because we like to be able to pass shared_obj's directly to functions taking refs):
    //const Dummy& dummy_ref_dangerous = createDummy();
    //(void)dummy_ref_dangerous;
  }


  return 0;

}

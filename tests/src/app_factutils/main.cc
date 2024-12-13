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

#include "NCrystal/internal/fact_utils/NCFactoryUtils.hh"
namespace NC = NCrystal;
#include <vector>
#include <iostream>

namespace {
  //Dummy factory, integers (N,M) -> vector of values of length N and M apart:
  class MyFunkyFactory final : public NC::CachedFactoryBase<std::pair<int,int>,std::vector<int>,2> {
    const char* factoryName() const final { return "MyFunkyFactory"; }
    std::string keyToString( const key_type& key ) const final
    {
      int N = key.first;
      int M = key.second;
      std::ostringstream ss;
      ss<<"(N="<<N<<";M="<<M<<")";
      return ss.str();
    }
  protected:
    virtual ShPtr actualCreate( const key_type& key ) const override
    {
      int N = key.first;
      int M = key.second;
      std::cout<<"    ....(actualCreate invoked for N="<<N<<" M="<<M<<")"<<std::endl;
      std::vector<int> v;
      v.reserve(N);
      v.push_back(0);
      while ((int)v.size()<N)
        v.push_back(v.back()+M);
      return std::make_shared<value_type>(std::move(v));
    }
  };

  static MyFunkyFactory s_funkyfact;

  std::shared_ptr<const std::vector<int>> test(int N, int M) {
    // auto stats = s_funkyfact.currentStats();
    // std::cout<<" [Factory stats now: nstrongrefs="<<stats.nstrongrefs<<", nweakrefs="<<stats.nweakrefs<<"]"<<std::endl;
    auto v = s_funkyfact.create(std::make_pair(N,M));
    std::cout<<"("<<N<<","<<M<<") -> [ ";
    for (auto e: *v)
      std::cout<<e<<", ";
    std::cout<<"]"<<std::endl;
    return v;
  }
}


int main()
{
  test(2,3);
  test(2,3);
  test(2,4);
  test(2,3);
  test(4,2);
  test(2,4);
  auto v= test(2,3);
  std::cout<<"(keeping this result in external shptr)"<<std::endl;
  test(3,3);
  test(4,2);
  test(2,4);
  test(2,3);
  test(2,3);
  test(5,10);
  test(2,3);
  test(2,3);
  v.reset();
  std::cout<<"(reset external shptr)"<<std::endl;
  test(2,3);
  test(3,3);
  test(2,3);
  test(7,1);
  test(3,3);
  test(2,3);
  auto vv= test(2,3);
  std::cout<<"(keeping this result in external shptr)"<<std::endl;
  std::cout<<"(invoking global clearCaches)"<<std::endl;
  NC::clearCaches();//NB: this used not not clear the factory weak ptrs but now it does, hence the test output changed.
  auto vvv = test(2,3);
  test(3,3);
  test(7,1);
  test(7,3);
  test(2,3);
  std::cout<<"(invoking factory specific cleanup)"<<std::endl;//NB: this used to be completeCleanup!
  s_funkyfact.cleanup();
  test(2,3);
  test(3,3);
  test(7,1);
  test(7,3);
  test(2,3);
  test(2,3);



  //MT stuff:

  // std::vector<std::thread> workers;
  // for (int i = 0; i < 5; ++i) {
  //   workers.push_back(std::thread([i]()
  //                                 {
  //                                   std::cout << "thread function"<<i<<"\n";
  //                                 }));
  // }


}

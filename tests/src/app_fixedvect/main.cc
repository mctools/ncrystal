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

#include "NCrystal/core/NCTypes.hh"
#include <iostream>
namespace NC = NCrystal;

class V3 : public NC::StronglyTypedFixedVector<V3,double,3> {
public:
  using StronglyTypedFixedVector::StronglyTypedFixedVector;
};

class V3alt : public NC::StronglyTypedFixedVector<V3alt,double,3> {
public:
  using StronglyTypedFixedVector::StronglyTypedFixedVector;
};

void v3fct(const V3& v)
{
  std::cout << "v3fct: "<< v << std::endl;

}

int main() {

  std::cout<<V3{1.1,2.1,3.1}<<std::endl;
  V3 a{1.1,2.1,3.1};
  std::cout << a << std::endl;
  V3 a2 = a;
  V3 a3(a);
  V3 a4{a};
  std::cout << a2 << std::endl;
  std::cout << a3 << std::endl;
  std::cout << a4 << std::endl;
  V3 b = {1.1,2.1,3.1};
  std::cout << b << std::endl;
  v3fct(V3{1.1,2.1,3.1});
  v3fct({1.1,2.1,3.1});
  v3fct({1.1,2,3.1});
  V3alt& c = b.as<V3alt>();
  std::cout << c << std::endl;
  c[1] = 15.123;
  std::cout << b << std::endl;
  v3fct(c.as<V3>());

  NC::NeutronDirection& nd = b.as<NC::NeutronDirection>();
  nd.at(2)=12.0;
  std::cout << nd << std::endl;
  nd.rawArray()[1] = 6.0;
  std::cout << nd << std::endl;
  nd.array()[0] = 3.0;
  std::cout << nd << std::endl;
  std::cout << b << std::endl;

  return 0;
}

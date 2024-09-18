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

#include "NCrystal/internal/NCIter.hh"
#include "NCrystal/internal/NCSpan.hh"
#include <iostream>
#include <vector>
#include <set>
#include <string>
namespace NC = NCrystal;

int main () {
  int arr[] = { 0,1,2,3,4 };
  const std::vector<int> vec{ 17,6,7,8,9 };
  std::vector<int> ncvec{ 17,6,7,8,9 };
  std::set<std::string> strset{ "bla","blu","bleeee" };

  std::cout << "Const array enumeration:\n";
  for (auto&& e : NC::enumerate(arr))//changed "const auto& e" -> "auto&&e" for clang 12.0.0
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
  std::cout << "\n";
  std::cout << "array enumeration:\n";
  for (auto e : NC::enumerate(arr))
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
  std::cout << "\n";
  std::cout << "array enumeration:\n";
  for (auto&& e : NC::enumerate(arr))
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
  std::cout << "\n";

  std::cout << "Const vector enumeration:\n";
  for (auto&& e : NC::enumerate(vec))//changed "const auto& e" -> "auto&&e" for clang 12.0.0
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
  std::cout << "Const vector enumeration:\n";
  for (auto&& e : NC::enumerate(vec))//changed "const auto& e" -> "auto&&e" for clang 12.0.0
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
  std::cout << "Const vector enumeration:\n";
  for (auto&& e : NC::enumerate(vec))//changed "const auto& e" -> "auto&&e" for clang 12.0.0
    std::cout << "(" << e.idx << ", " << e.val << ")\n";

  std::cout << "Vector enumeration:\n";
  for (auto e : NC::enumerate(ncvec)) {
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
    e.val += 100;
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
  }
  std::cout << "Vector enumeration:\n";
  for (auto e : NC::enumerate(ncvec))
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
  std::cout << "Vector enumeration:\n";
  for (auto&& e : NC::enumerate(ncvec)) {
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
    e.val += 1000;
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
  }
  std::cout << "Vector enumeration:\n";
  for (auto e : NC::enumerate(ncvec))
    std::cout << "(" << e.idx << ", " << e.val << ")\n";

  std::cout << "Str set enumeration:\n";
  for (auto&& e : NC::enumerate(strset)) {
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
  }

  std::cout << "initializer list enumeration:\n";
  auto hackhack={"aaa","bbb","ccc"};//FIXME: We have to do it like this to extend the lifetime of the container in this case. Solution?
  for (auto&& e : NC::enumerate(hackhack)) {
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
  }
  std::cout << "std::string enumeration:\n";
  auto hackhack_str = std::string("abcdefghijkl");//FIXME: We have to do it like this to extend the lifetime of the std::string in this case. Solution?
  for (auto&& e : NC::enumerate(hackhack_str)) {
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
  }
  std::cout << "const char* enumeration:\n";
  for (auto&& e : NC::enumerate("abcdefghijkl")) {
    std::cout << "(" << e.idx << ", " << (e.val?e.val:'?') << ")\n";
  }
  std::cout << "span enumeration:\n";
  auto hackhack2=NC::Span<int>(ncvec);//FIXME: We have to do it like this to extend the lifetime of the span in this case. Solution?
  for (auto&& e : NC::enumerate(hackhack2)) {
    std::cout << "(" << e.idx << ", " << e.val << ")\n";
  }


  //ValRange<T> ncrange(T n) { return ValRange<T>(n); }
  std::cout<<"test ncrange"<<std::endl;
  for ( auto i : NC::ncrange(5) )
    std::cout<<i<<std::endl;
  std::cout<<"---"<<std::endl;
  for ( auto i : NC::ncrange(15,17) )
    std::cout<<i<<std::endl;
  std::cout<<"---"<<std::endl;
  for ( auto i : NC::ncrange(uint8_t(1),uint8_t(4)) ) {
    static_assert(std::is_same<decltype(i),uint8_t>::value,"");
    std::cout<<(int)i<<std::endl;
  }
  std::cout<<"---"<<std::endl;
  for ( auto i : NC::ncrange(-10,3) ) {
    std::cout<<i<<std::endl;
  }
  std::cout<<"---"<<std::endl;
  for ( auto i : NC::ncrange(-10,3) )
    std::cout<<i<<std::endl;

  std::cout<<"---"<<std::endl;
  for ( auto i : NC::ncrange((unsigned)5) )
    std::cout<<i<<std::endl;

  // std::cout<<"---"<<std::endl;
  // for ( auto i : NC::ncrange((uint16_t)2,(uint64_t)5) )
  //   std::cout<<i<<std::endl;

}

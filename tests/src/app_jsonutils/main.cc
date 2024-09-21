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

#include "NCrystal/internal/NCString.hh"
#include <iostream>
namespace NC = NCrystal;

namespace {
  void printObject( const char * t) { std::cout << t; }
  void printObject( const std::string& t) { std::cout << t; }
  void printObject( double t) { std::cout << t; }
  void printObject( float t) { std::cout << t; }
  template<class T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
  void printObject( T t ) {  std::cout << t; }

  template<typename TContainer, typename T = typename TContainer::value_type>
  void printObject( const TContainer& c)
  {
    std::cout <<"{";
    for ( auto & e : c ) {
      printObject(e);
      std::cout<<",";
    }
    std::cout<<"}";
  }

  template<class T>
  void test(const T& t )
  {
    std::cout << "Test object: >>>";
    printObject(t);
    std::cout << "<<<\n";
    std::cout << "As JSON: >>>";
    NC::streamJSON( std::cout, t );
    std::cout << "<<<\n";
  }
}

int main()
{
  test("Hello world");
  test("Hello world\\n");
  test(std::string("Hello world\\n"));
  test("Hello '\"{}\\world\n\tyo\xc3\x85\xc3\xa6\xc3\xb8");
  //                            ^^^^^^^^^^^^^^^^^^^^^^^^
  // The last three unicode chars above (two bytes each) are uppercase danish A
  // with ring lowercase danish ae lowercase danish o-slash

  test(17.5);
  test(NC::kPi);
  test(-17);
  test(1234567891234567ull);
  test((uint16_t)17);
  test((int16_t)-17);
  test(-17.5e3f);
  test(NC::VectS{"hello","world\n\tyo\xc3\x85\xc3\xa6\xc3\xb8","bla"});
  test(NC::VectS{});
  test(NC::VectD{1e-100,std::numeric_limits<double>::infinity()});
  test(NC::VectD{});
  test(std::vector<NC::VectS>{{"hello"},{"a","b","c"},{}});


  return 0;
}

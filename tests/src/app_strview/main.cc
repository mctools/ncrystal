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

#include "NCrystal/internal/utils/NCStrView.hh"
#include <iostream>

namespace NC = NCrystal;

namespace {
  void fctA( NC::StrView sv ) { std::cout<<sv<<std::endl; }

  template <class TVect>
  void printVect( const TVect& v )
  {
    std::cout<<'[';
    for ( auto i : NC::ncrange(v.size()) )
      std::cout<<(i?", ":" ") << '"'<< v.at(i)<<'"';
    std::cout<<" ]"<<std::endl;
  }
}



void test_splitting()
{
  auto test = [](NC::StrView sv )
  {
    std::cout<<"SV(\""<<sv<<"\")::split() = "; printVect(sv.split());
    std::cout<<"SV(\""<<sv<<"\")::split(';') = "; printVect(sv.split(';'));
    std::cout<<"SV(\""<<sv<<"\")::splitTrimmed(';') = "; printVect(sv.splitTrimmed(';'));
    std::cout<<"SV(\""<<sv<<"\")::splitTrimmedNoEmpty(';') = "; printVect(sv.splitTrimmedNoEmpty(';'));
  };
  test(" bla bla \r  \t \n\n yi ha   ");
  test(" bla bla ; ha   ");
  test(" bla bla ;;  ; ha   ");
  test("");
  test("  ");
  test(" ;  ");
  test("; ");
  test(" ;");
}

int main() {
  fctA( "test str lit" );
  std::string stdstr("test std::string");
  fctA( stdstr );
  stdstr += ".c_str()";
  fctA( stdstr.c_str() );
  //Does not compile (on purpose, deleted function!):  fctA( std::move(stdstr) );
  fctA( NC::StrView::make("test strview") );
  NC::StrView sv2 = NC::StrView::make("test strview2");
  fctA( sv2 );
  NC::StrView sv3 = NC::StrView::make("test strview3");
  fctA( std::move(sv3) );
  {
    char test_arr[13] = {'t','e','s','t',' ','c','h','a','r','[','1','3',']'};
    fctA( test_arr );
  }
  {
    char test_arr[14] = {'t','e','s','t',' ','c','h','a','r','[','1','4',']','\0'};
    fctA( test_arr );
  }
  {
    char test_arr[14] = {'t','e','s','t',' ','c','h','a','r',' ','p','t','r','\0'};
    const char * tp = &test_arr[0];
    fctA( NC::StrView(tp, 13) );
  }
  using SV = NC::StrView;
  constexpr auto lala = SV::make("bla");
  std::string lalastr("bla");
  nc_assert_always( lala == SV(lalastr) );
  nc_assert_always( SV::make("bla").to_string() == "bla" );
  nc_assert_always( SV::make("bla") == SV::make("bla") );
  nc_assert_always( SV::make("bla") != SV::make("bla2") );
  nc_assert_always( SV::make("bla") < SV::make("bla2") );
  nc_assert_always( SV::make("aa") < SV::make("bb") );
  nc_assert_always( SV::make("aa") < SV::make("aaa") );
  nc_assert_always( SV::make("A") < SV::make("a") );

  nc_assert_always( ! ( SV::make("bla") < SV::make("bla") ) );
  nc_assert_always( ! ( SV::make("atomdb") < SV::make("atomdb") ) );

  //Try with following char being different:
  auto atomdb1 = SV::make("atomdbA").substr(0,6);
  auto atomdb2 = SV::make("atomdbB").substr(0,6);
  nc_assert_always( ! ( atomdb1 != atomdb2 ) );
  nc_assert_always( atomdb1 == atomdb2 );
  nc_assert_always( ! ( atomdb1 < atomdb2 ) );
  nc_assert_always( ! ( atomdb2 < atomdb1 ) );


  nc_assert_always( SV::make(" \n \r \t aa     \r \t\t\t\t \n\n \r\n").trimmed() == SV::make("aa") );
  nc_assert_always( SV::make("abcd").substr(0) == SV::make("abcd") );
  nc_assert_always( SV::make("abcd").substr(1) == SV::make("bcd") );
  nc_assert_always( SV::make("abcd").substr(2) == SV::make("cd") );
  nc_assert_always( SV::make("abcd").substr(3) == SV::make("d") );
  nc_assert_always( SV::make("abcd").substr(4) == SV::make("") );
  nc_assert_always( SV::make("abcd").substr(5) == SV::make("") );
  nc_assert_always( SV::make("abcd").substr(0,2) == SV::make("ab") );
  nc_assert_always( SV::make("abcd").substr(0,4) == SV::make("abcd") );
  nc_assert_always( SV::make("abcd").substr(0,5) == SV::make("abcd") );
  nc_assert_always( SV::make("abcd").substr(5,1) == SV::make("") );

  nc_assert_always( SV::make("abcd").contains('a') );
  nc_assert_always( SV::make("abcd").contains('d') );
  nc_assert_always( !SV::make("abcd").contains('\0') );
  nc_assert_always( !SV::make("abcd").contains_any("efg") );
  nc_assert_always( SV::make("abcd").contains_any("efgc") );
  nc_assert_always( !SV::make("").contains_any("") );
  nc_assert_always( !SV::make("").contains_any("d") );

  nc_assert_always( SV::make("").find(' ') == SV::npos );
  nc_assert_always( SV::make("hello there").find(' ') == 5 );
  nc_assert_always( SV::make("hello there").find('e') == 1 );
  nc_assert_always( SV::make("abcd").find('d') == 3 );

  nc_assert_always( SV::make("").find(" ") == SV::npos );
  nc_assert_always( SV::make("hello there").find(" ") == 5 );
  nc_assert_always( SV::make("hello there").find("e") == 1 );
  nc_assert_always( SV::make("abcd").find("d") == 3 );

  nc_assert_always( SV::make("").find(SV::make(" ")) == SV::npos );
  nc_assert_always( SV::make("hello there").find(SV::make(" ")) == 5 );
  nc_assert_always( SV::make("hello there").find(SV::make("e")) == 1 );
  nc_assert_always( SV::make("abcd").find(SV::make("d")) == 3 );

  nc_assert_always( SV::make("abcd").find_first_of("dc") == 2 );
  nc_assert_always( SV::make("abcdcdcdcdcd").find_first_of("dc") == 2 );
  nc_assert_always( SV::make("akjbdslfjo;j vpwnsdncsndv").find_first_of("1232344545456") == SV::npos );

  {
    std::string tmp("WUHUYIR");
    for ( auto i : NC::ncrange(1500) ) {
      (void) i;
      tmp += "A";
    }
    tmp += "GGG";
    nc_assert_always( SV(tmp).contains_any("R") );
    nc_assert_always( SV(tmp).contains_any("A") );
    nc_assert_always( !SV(tmp).contains_any("B") );
    nc_assert_always( SV(tmp).contains_any("!Rt") );
    nc_assert_always( SV(tmp).contains_any("!Gt") );
    nc_assert_always( !SV(tmp).contains_any("!Bt") );

    nc_assert_always( SV(tmp).find_first_of("R") == 6);
    nc_assert_always( SV(tmp).find_first_of("A") == 7 );
    nc_assert_always( SV(tmp).find_first_of("B") == SV::npos );
    nc_assert_always( SV(tmp).find_first_of("!Rt") == 6 );
    nc_assert_always( SV(tmp).find_first_of("!Gt") == 1507 );
    nc_assert_always( SV(tmp).find_first_of("!Bt") == SV::npos );


  }

  test_splitting();

}

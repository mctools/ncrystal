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
#include "NCrystal/internal/NCMath.hh"
#include <iostream>
#include <limits>

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(__LINE__ && x)

void dump(const char* settings,
          const std::string& s,
          const std::vector<std::string>& result)
{
  std::cout << '"' << s << "\" splits to "<<result.size()<<" parts (settings: "<<settings<<"):"<<std::endl;
  for (std::size_t i = 0; i<result.size(); ++i)
    std::cout << "       -->  \"" << result.at(i) << '"' << std::endl;
}

void test_split()
{
  std::vector<std::string> parts;
  std::string s;

  s="";
  NC::split(parts,s);
  dump("no options",s,parts);

  s="";
  NC::split(parts,s,0,';');
  dump("sep=;",s,parts);

  s="   ";
  NC::split(parts,s,0,';');
  dump("sep=;",s,parts);

  s="   ";
  NC::split(parts,s);
  dump("no options",s,parts);

  s="  hej med dig \t lala\t";
  NC::split(parts,s);
  dump("no options",s,parts);

  s=";hej;med;;dig;";
  NC::split(parts,s,0,';');
  dump("sep=;",s,parts);

  s=";hej;med; dig;;;hejsa";
  NC::split(parts,s,0,';');
  dump("sep=;",s,parts);

  s="hej;med;dig;";
  NC::split(parts,s,1,';');
  dump("maxsplit 1;sep=;",s,parts);

  s="hej med dig";
  NC::split(parts,s,1);
  dump("maxsplit 1",s,parts);

  s=" hej med dig";
  NC::split(parts,s,1);
  dump("maxsplit 1",s,parts);

  s=" hej med dig ";
  NC::split(parts,s,1);
  dump("maxsplit 1",s,parts);

  s="temp=77.0";
  NC::split(parts,s,0,'=');
  dump("sep='='",s,parts);
}

void require_bad_dbl(const char * c) {
  bool bad = false;
  try {
    NC::str2dbl(c);
  } catch ( NC::Error::BadInput& e ) {
    bad = true;
  }
  REQUIRE(bad);
}
void require_bad_int(const char * c) {
  bool bad = false;
  try {
    NC::str2int(c);
  } catch ( NC::Error::BadInput& e ) {
    bad = true;
  }
  REQUIRE(bad);
}

void test_misc() {
  REQUIRE(!NC::isSimpleASCII("something\tsdf",NC::AllowTabs::No,NC::AllowNewLine::No));
  REQUIRE(!NC::isSimpleASCII("something\nsdf",NC::AllowTabs::No,NC::AllowNewLine::No));
  REQUIRE(!NC::isSimpleASCII("something\rsdf",NC::AllowTabs::No,NC::AllowNewLine::No));
  REQUIRE( NC::isSimpleASCII("something\t sdf",NC::AllowTabs::Yes,NC::AllowNewLine::No));
  REQUIRE( NC::isSimpleASCII("something\n sdf",NC::AllowTabs::No,NC::AllowNewLine::Yes));
  REQUIRE( NC::isSimpleASCII("something\r sdf",NC::AllowTabs::No,NC::AllowNewLine::Yes));
  REQUIRE( NC::isSimpleASCII("some\tthing\n sdf",NC::AllowTabs::Yes,NC::AllowNewLine::Yes));

  std::string t;
  t = "bla";NC::trim(t);REQUIRE(t=="bla");
  t = "\rbla";NC::trim(t);REQUIRE(t=="bla");
  t = "bla\t\n";NC::trim(t);REQUIRE(t=="bla");
  t = "\t\nb la\r";NC::trim(t);REQUIRE(t=="b la");
  t = "\t\nb l\t\na\r";NC::trim(t);REQUIRE(t=="b l\t\na");
  t = "b  la";NC::trim(t);REQUIRE(t=="b  la");
  t = " b l a ";NC::trim(t);REQUIRE(t=="b l a");

  REQUIRE(NC::startswith("hejsa","hej"));
  REQUIRE(!NC::startswith(" hejsa","hej"));
  REQUIRE(NC::startswith("hejsa",""));
  REQUIRE(NC::startswith("",""));
  REQUIRE(!NC::startswith("\n","\r"));
  REQUIRE(!NC::startswith("hej","hejsa"));

  REQUIRE(NC::contains("hejsa",'h'));
  REQUIRE(NC::contains("hejsa","h"));
  REQUIRE(NC::contains("hejsa","hej"));
  REQUIRE(NC::contains("hejsa","ejsa"));
  REQUIRE(!NC::contains("hejsa\n sdfsdf ds f$%^&^  \n lala","t"));
  REQUIRE(!NC::contains("hejsa\n sdfsdf ds f$%^&^  \n lala",'t'));
  REQUIRE(NC::contains("hejsa\n sdfsdf ds f$%^&^  \n lala","$"));
  REQUIRE(NC::contains("hejsa\n sdfsdf ds f$%^&^  \n lala",'$'));
  REQUIRE(!NC::contains("hejsa",'\0'));
  REQUIRE(!NC::contains("",'\0'));
  REQUIRE(NC::contains("",""));
  REQUIRE(NC::contains("sdfsdf",""));

  //Check if any of the chars in "needles" is present in the string (haystack):
  REQUIRE(NC::contains_any("sdfsdf","456457d"));
  REQUIRE(!NC::contains_any("","456457d"));
  REQUIRE(!NC::contains_any("sdfsd",""));
  REQUIRE(!NC::contains_any("hejsa","HEJSA"));
  REQUIRE(NC::contains_any("hej med dig!","'{}![]"));

  REQUIRE(NC::str2dbl("2.0")==2.0);
  REQUIRE(NC::str2dbl("2.0e-3")==2.0e-3);
  REQUIRE(NC::str2dbl("-2.0E-3")==-2.0e-3);
  REQUIRE(NC::str2dbl("-1E-3")==-1.0e-3);
  require_bad_dbl("2.0a");
  require_bad_dbl(" 2.0");
  //REQUIRE(NC::str2dbl(" 2.0")==2.0);//seems to work. Should we disallow?... disallowed now
  require_bad_dbl("2.0 ");
  require_bad_dbl("2e.0");
  require_bad_dbl("e-3");

  //Not sure if these should actually work or not, but they seem not to with gcc
  //6.3.1 so we put a test here for now. If they fail on some platforms, we
  //might have to add code in NC::str2dbl looking for strings like infinity
  //strings (to get consistent behaviour everywhere):

  REQUIRE(NC::ncisinf(NC::str2dbl("inf")));
  require_bad_dbl("\tinf");
  require_bad_dbl("\tINF");
  require_bad_dbl("  \t  inf");
  require_bad_dbl(" \t INF");
  require_bad_dbl("inf ");
  require_bad_dbl("INF ");
  require_bad_dbl("infinity");
  require_bad_dbl("INFINITY");

  REQUIRE(NC::str2int("-1")==-1);
  //DISALLOWED: REQUIRE(NC::str2int(" 2")==2);//seems to work. Should we disallow?
  REQUIRE(NC::str2int("34545")==34545);
  REQUIRE(NC::str2int("000")==0);
  require_bad_int("2.0");
  require_bad_int("1e2");
  require_bad_int("0.1");
  require_bad_int("2 ");
  require_bad_int(" 2");
}

int main(int,char**) {
  test_split();
  test_misc();

  return 0;
}

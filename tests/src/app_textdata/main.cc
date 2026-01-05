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

#include "NCrystal/text/NCTextData.hh"
#include "NCTestUtils/NCTestFindData.hh"
#include <iostream>
namespace NC = NCrystal;

int main()
{
  //Test nctest::find_test_data here, just to test it somewhere (asserts
  //internally that files are found):
  nctest::find_test_data("dummy_D.cif");
  nctest::find_test_data("refnc1","Al_sg225.ncmat");

  //The actual textdata test:
  static_assert(std::is_constructible<std::string,const char*>::value,"");
  static_assert(std::is_constructible<NC::DataSourceName,const char*>::value,"");
  static_assert(std::is_convertible<const char*,std::string>::value,"");
  static_assert(std::is_convertible<std::string,NC::DataSourceName>::value,"");

  auto TextDataFromRaw = [](const char *ptr,const char * dsn)// -> NC::TextData
  {
    NC::Optional<NC::DataSourceName> optDSN;
    if (dsn)
      optDSN = dsn;
    return NC::TextData( NC::RawStrData(NC::RawStrData::static_data_ptr_t(),ptr),
                         NC::TextData::DataType{"dummydatatype"},
                         optDSN );
  };

  auto doprint = []( const NC::TextData& td, const char * dsn = nullptr ) {
    std::cout<<">>>";
    if (dsn!=nullptr)
      std::cout<<" ("<<dsn<<")";
    std::cout<<":"<<std::endl;
    for ( auto l : td )
      std::cout<<"    >>>>>"<<l<<"<<<<<"<<std::endl;
  };
  auto test = [&doprint,&TextDataFromRaw](const char * data, const char * dsn = nullptr) {
    doprint(TextDataFromRaw(data,"/bla/bla.txt"), dsn);
  };
  auto testThrow = [&TextDataFromRaw](const char * data, const char * dsn = nullptr) {
    try {
      for ( auto& l : TextDataFromRaw(data,"/bla/bla.txt") ) { (void)l; }
    } catch ( NC::Error::BadInput& ) {
      std::cout<<">>> BadInput as expected: "<<dsn<<""<<std::endl;
      return;
    }
    nc_assert_always(0&&"did not get expected DataLoadError");
  };

  doprint(TextDataFromRaw("data with no descr",nullptr));


  test("bla bla bla\nblabla\r\nhello you\r\nblabla\n","both unix and dos line endings");

  test("bla\n","single line with single word (dos)");
  testThrow("bla\r","single line with single word (mac)");
  test("bla\r\n","single line with single word (win)");
  test("bla\0","single line with single word (no line ending)");
  test("","empty file");
  test("\0","empty file with extra null char");
  test("\n\0","single empty line (unix)");
  test("\r\n\0","single empty line (dos)");
  testThrow("\r\0","single empty line (mac)");
  test("\n \n\0","empty line and line with 1 space (dos)");
  test("\r\n \r\n\0","empty line and line with 1 space (win)");
  testThrow("\r \r\0","empty line and line with 1 space (mac)");
  test("\n \0","empty line and line with 1 space (dos + no final line ending");

  test("high bit set: \xE2\x98\xA0\n","char with high bit set (UTF8 multi-byte)");
  test("a tab\t\r\ntwo tabs\t\ttada.\n","string with TAB character (a valid non-newline non-null char without 4 high bits set)");

  test("\n\r\n ","three lines, mixed endings (dos+win+null)");
  testThrow("\n\r \n ","four lines, mixed endings (dos+mac+dos+null)");

  std::string::size_type lastcapacity = 0;

  for ( const auto& line :  TextDataFromRaw("123\n1234567\n1234567890123456789012345678901234567890\n123\n444",nullptr)) {
    std::cout<<"Got line >>>"<<line;
    std::cout<<"<<< [length="<<line.size()<<", buf capacity "<<"<not shown>"<<"]"<<std::endl;
    nc_assert_always(line.capacity() >= lastcapacity);
    lastcapacity = line.capacity();
    nc_assert_always(*(&line[0] + line.size()) == '\0');
  }
  std::string s100( 100, 'a' ), s200(s100+s100), s300(s200+s100);
  std::string largedata = "";
  largedata += s100+'\n';
  largedata += s200+'\n';
  largedata += s300+'\n';
  largedata += s200+'\n';
  largedata += s100+'\n';
  largedata += '\n';
  largedata += s300+'\n';
  largedata += s300+s300+s300+'\n';
  lastcapacity = 0;
  for ( const auto& line :  NC::TextData(NC::RawStrData(std::move(largedata)),
                                         NC::TextData::DataType{"dummy"})) {
    std::cout<<"Got line >>>";
    if ( line.size() < 30 ) {
      std::cout<<line;
    } else {
      std::cout<<line.substr(0,5)<<"..."<<line.substr(line.size()-6);
    }
    std::cout<<"<<< [length="<<line.size()<<", buf capacity "<<(line.capacity()<50?"<50":">=50")<<"]"<<std::endl;
    nc_assert_always(line.capacity() >= lastcapacity);
    lastcapacity = line.capacity();
    nc_assert_always(*(&line[0] + line.size()) == '\0');
  }

  std::cout<<largedata<<std::endl;

  {
    NC::TextData td = TextDataFromRaw("123\n1234567\n1234567890123456789012345678901234567890\n123\n444",nullptr);
    auto it = td.begin();
    auto itE = td.end();
    std::cout<<"Line 0: "<<*it<<std::endl;
    auto it2 = it;
    std::cout<<"Line 0 copy: "<<*it2<<std::endl;
    auto it3(it);
    std::cout<<"Line 0 copy: "<<*it3<<std::endl;
    auto it4(std::move(it));
    std::cout<<"Line 0 copy: "<<*it4<<std::endl;
    for ( auto it5 = std::move(it4); it5 != itE; ++it5 )
      std::cout<<"All lines: "<<*it5<<std::endl;
    ++it2;
    ++it3;
    nc_assert_always(it2==it3);
  }

  return 0;
}

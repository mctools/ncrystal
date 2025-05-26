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

#include "NCrystal/internal/cfgutils/NCCfgTypes.hh"
#include <iostream>

namespace NC = NCrystal;

namespace {
  using P = std::pair<std::string,std::string>;
  using PV = std::vector<P>;
  using KVMap = NC::Cfg::CfgKeyValMap;

  void test( PV data_holder )
  {
    KVMap::DecodedData data;
    for ( auto& e : data_holder )
      data.emplace_back(e.first,e.second);

    NC::Cfg::VarId dummy_varid = static_cast<NC::Cfg::VarId>( 10 );
    KVMap m( data, dummy_varid );

    std::cout<<"m: ";
    m.stream(std::cout);
    std::cout<<std::endl;
    m.streamJSON(std::cout);
    std::cout<<std::endl;
    auto varbuf = m.encoded();
    std::cout<<"m2: ";
    KVMap m2( varbuf );
    m2.stream(std::cout);
    std::cout<<std::endl;

    m2.streamJSON(std::cout);
    std::cout<<std::endl;


    KVMap::DecodedData data_sorted = data;
    std::stable_sort( data_sorted.begin(), data_sorted.end() );
    nc_assert_always( m2.decoded() == data_sorted );
    nc_assert_always( m2.varId() == dummy_varid );
    nc_assert_always( m.decoded() == data_sorted );
    nc_assert_always( m.varId() == dummy_varid );
    for ( auto& e : data ) {
      nc_assert_always( m.hasValue(e.first) );
      nc_assert_always( m.getValue(e.first) == e.second );
    }
    nc_assert_always( !m.hasValue("sddslkf234324cd3") );
    nc_assert_always( m.getValue("sddslkf234324cd3", "foo" )=="foo" );

    bool fail = false;
    try {
      m.getValue("sddslkf234324cd3");
    } catch ( NC::Error::BadInput& e ) {
      std::cout<<"Got expected NCBadInput error: "<<e.what()<<std::endl;
      fail = true;
    }
    if (!fail)
      throw std::runtime_error("Did not end in error as it should");
  }

  void test_bad( PV data_holder )
  {
    KVMap::DecodedData data;
    for ( auto& e : data_holder )
      data.emplace_back(e.first,e.second);
    try {
      KVMap m( data, NC::Cfg::VarId( 10 ) );
    } catch ( NC::Error::BadInput& e ) {
      std::cout<<"Got expected NCBadInput error: "<<e.what()<<std::endl;
      return;
    }
    throw std::runtime_error("Did not end in error as it should");
  }

  void test_lifetime()
  {
    NC::Cfg::VarId dummy_varid = static_cast<NC::Cfg::VarId>( 17 );
    std::string a("aaa"), b("bbb");
    KVMap::DecodedData dd;
    dd.emplace_back(a,b);
    KVMap mdd( dd, dummy_varid );
    //Overwrite original string content:
    a = "ccc";
    b = "ddd";
    std::cout<<"mdd: ";
    mdd.stream(std::cout);
    std::cout<<std::endl;
    nc_assert_always( mdd.decoded().size() == 1 );
    nc_assert_always( mdd.decoded().at(0).first == "aaa" );
    nc_assert_always( mdd.decoded().at(0).second == "bbb" );
  }
}

int main(int,char**) {
  test( {} );
  test( { P("foo", "bar"), P("a","b") } );
  test( { P("foo", "bar") } );
  test_bad( { P("foo", "bar"), P("foo","bar2") } );
  test_bad( { P("fo o", "bar") } );
  test_bad( { P("fo-o", "bar") } );
  test_bad( { P("foo", "bar\t") } );
  test( { P("foo", "bar 1.e24324 +- dsfsdf") } );
  test( { P("foo", "bar"),
        P("a","b"),
        P("a1","b"),
        P("a2","b"),
        P("a3","b"),
        P("a4","b"),
        P("a5","b"),
        P("a6","b"),
        P("a7","b"),
        P("a8","b"),
        P("a9","b"),
        P("a10","b"),
        P("a11","b"),
        P("a12","b"),
        P("a14","b"),
      } );
  test( { P("aaa", "aaaval"), P("bbb","bbbval") } );
  test( { P("bbb","bbbval"), P("aaa", "aaaval") } );

  test( { P("foo111111111111111111111111111111111111111111111111111"
            "22222222222222222222222222222222222222222222222222222222"
            "3333333333333333333333333333333333333333333333333333333333"
            "4444444444444444444444444444444444444444444444444444444444"
            "5555555555555555555555555555555555555555555555555555555555",
            "baraaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
            "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
            "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc")
    } );


  test_lifetime();

  return 0;
}

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

// Tests of the sgsym component (space groups and their symmetries).

#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include <iostream>
#include <sstream>
#include <set>

namespace NC = NCrystal;
#define REQUIRE(x) nc_assert_always(x)

////////////////////////////////////////////////////////////////////////////////
// Space group table
////////////////////////////////////////////////////////////////////////////////

//Tests the SpaceGroup and SpaceGroupHallNumber classes. The complete table of
//space group settings is printed, so the reference log also guards against
//accidental changes to it.

namespace test_table {

namespace {
  void expectBad( const char * s )
  {
    try {
      NC::SpaceGroup sg( s );
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  \"" << s << "\" -> BadInput: " << e.what() << std::endl;
      return;
    }
    nc_assert_always( false );
  }
  void expectBadHall( std::uint16_t hn )
  {
    try {
      NC::SpaceGroup sg{ NC::SpaceGroupHallNumber{ hn } };
    } catch ( NC::Error::BadInput& e ) {
      std::cout << "  Hall number " << hn << " -> BadInput: " << e.what()
                << std::endl;
      return;
    }
    nc_assert_always( false );
  }
}

void run()
{
  static_assert( sizeof( NC::SpaceGroup ) == 2, "" );
  static_assert( std::is_trivially_copyable<NC::SpaceGroup>::value, "" );
  //Constructors are explicit (no implicit conversions):
  static_assert( !std::is_convertible<const char*,NC::SpaceGroup>::value, "" );
  static_assert( !std::is_convertible<std::string,NC::SpaceGroup>::value, "" );
  static_assert( !std::is_convertible<NC::StrView,NC::SpaceGroup>::value, "" );
  static_assert( !std::is_convertible<NC::SpaceGroupHallNumber,
                                      NC::SpaceGroup>::value, "" );
  static_assert( std::is_constructible<NC::SpaceGroup,std::string&&>::value,
                 "" );

  std::cout << "All settings (Hall number, setting, Hall symbol, flags):"
            << std::endl;
  std::set<unsigned> numbers, numbers_multi;
  unsigned last_number = 0;
  for ( std::uint16_t hn = 1; hn <= 530; ++hn ) {
    const NC::SpaceGroup sg{ NC::SpaceGroupHallNumber{ hn } };
    REQUIRE( sg.hallNumber().get() == hn );
    REQUIRE( sg.number() >= 1 && sg.number() <= 230 );
    REQUIRE( sg.number() >= last_number );//sorted
    REQUIRE( sg.isDefaultSetting() == ( sg.number() != last_number ) );
    last_number = sg.number();
    numbers.insert( sg.number() );
    if ( sg.hasMultipleSettings() )
      numbers_multi.insert( sg.number() );
    REQUIRE( *sg.hallSymbol() );
    //Empty choice codes only for default settings:
    REQUIRE( sg.choice() != nullptr );
    REQUIRE( *sg.choice() || sg.isDefaultSetting() );//(*: non-empty?)
    //String round trip and stream output:
    const std::string str = sg.toString();
    REQUIRE( NC::SpaceGroup( str ) == sg );
    REQUIRE( NC::SpaceGroup( std::string( str ) ) == sg );//temporary
    REQUIRE( NC::SpaceGroup( NC::StrView( str ) ) == sg );
    REQUIRE( NC::SpaceGroup( str.c_str() ) == sg );
    //Hall number -> SpaceGroup -> string -> SpaceGroup -> Hall number:
    REQUIRE( NC::SpaceGroup( str ).hallNumber().get() == hn );
    std::ostringstream ss;
    ss << sg;
    REQUIRE( ss.str() == str );
    //Bare number gives default setting:
    const std::string nstr = std::to_string( sg.number() );
    REQUIRE( ( NC::SpaceGroup( nstr ) == sg ) == sg.isDefaultSetting() );
    REQUIRE( ( NC::SpaceGroup( nstr ) < sg ) == !sg.isDefaultSetting() );
    std::cout << "  " << hn << " " << sg << " \"" << sg.hallSymbol() << "\""
              << ( sg.isDefaultSetting() ? " [default]" : "" )
              << ( sg.hasMultipleSettings() ? "" : " [single]" ) << std::endl;
  }
  REQUIRE( numbers.size() == 230 );
  std::cout << "Numbers with multiple settings: " << numbers_multi.size()
            << std::endl;
  REQUIRE( numbers_multi.size() == 90 );

  std::cout << "Examples:" << std::endl;
  for ( auto s : { "227", "227:1", "227:2", "166", "166:H", "166:R", "14",
                   "14:b1", "14:c1", "62", "62:cab", "68:1cab", "70:2", "15:-b2",
                   "1", "225" } ) {
    const NC::SpaceGroup sg( s );
    std::cout << "  \"" << s << "\" -> " << sg << " (Hall number "
              << sg.hallNumber() << ", Hall symbol \"" << sg.hallSymbol()
              << "\")" << std::endl;
  }
  REQUIRE( NC::SpaceGroup( "227" ).hallNumber().get() == 525 );
  REQUIRE( NC::SpaceGroup( "227:2" ).hallNumber().get() == 526 );

  std::cout << "Invalid input:" << std::endl;
  for ( auto s : { "", "0", "231", "999", "1000", "-1", "+1", " 227", "227 ",
                   "abc", "227:", "227:3", "227:1:2", "166:h", "62:abc",
                   "225:1", "14:b", "3:b1", ":2", "22a" } )
    expectBad( s );
  expectBadHall( 0 );
  expectBadHall( 531 );
  std::cout << "All OK" << std::endl;
}

}

int main()
{
  std::cout << "==== Space group table ====" << std::endl;
  test_table::run();
  return 0;
}

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

#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCSGSymHallTable.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    const SGSymHallTableEntry& entry( std::uint16_t hn )
    {
      return sgsymHallTableEntry( SpaceGroupHallNumber{ hn } );
    }

    //Hall number of the first (i.e. default) setting of a given number (the
    //table is sorted by number):
    std::uint16_t firstHallNumber( unsigned number )
    {
      nc_assert( number >= 1 && number <= 230 );
      for ( std::uint16_t hn = 1; hn <= SpaceGroupHallNumber::max_value; ++hn )
        if ( entry( hn ).number == number )
          return hn;
      nc_assert_always( false );
      return 0;
    }

    [[noreturn]] void throwBadSpec( StrView s, const char * extra )
    {
      NCRYSTAL_THROW2(BadInput,"Invalid space group specification: \""<<s
                      <<"\" ("<<extra<<")");
    }
  }
}

NC::SpaceGroup::SpaceGroup( SpaceGroupHallNumber hn )
  : m_hall( DoValidate, hn )
{
}

NC::SpaceGroup::SpaceGroup( StrView s )
{
  //Parse "<number>" or "<number>:<choice>":
  const auto parts = s.split<2>( ':' );
  const char * expected = ( "expected \"<number>\" or \"<number>:<choice>\""
                            " with a number in 1..230" );
  if ( parts.size() > 2 || parts.front().empty()
       || parts.front().size() > 3 )
    throwBadSpec( s, expected );
  unsigned number = 0;
  for ( char c : parts.front() ) {
    if ( !( c >= '0' && c <= '9' ) )
      throwBadSpec( s, expected );
    number = 10 * number + static_cast<unsigned>( c - '0' );
  }
  if ( !( number >= 1 && number <= 230 ) )
    throwBadSpec( s, expected );
  std::uint16_t hn = firstHallNumber( number );
  if ( parts.size() == 2 ) {
    const StrView choice = parts.back();
    while ( hn <= SpaceGroupHallNumber::max_value
            && entry( hn ).number == number
            && choice != StrView( entry( hn ).choice,
                                  std::strlen( entry( hn ).choice ) ) )
      ++hn;
    if ( choice.empty() || hn > SpaceGroupHallNumber::max_value
         || entry( hn ).number != number ) {
      std::ostringstream ss;
      ss << "invalid setting choice for space group " << number;
      const auto hn0 = firstHallNumber( number );
      if ( !SpaceGroup( SpaceGroupHallNumber{ hn0 } ).hasMultipleSettings() ) {
        ss << ", which has only one setting";
      } else {
        ss << ", available choices are:";
        for ( auto h = hn0; ( h <= SpaceGroupHallNumber::max_value
                              && entry( h ).number == number ); ++h )
          ss << ' ' << ( *entry( h ).choice//i.e. not the empty string
                         ? entry( h ).choice : "(none)" );
      }
      throwBadSpec( s, ss.str().c_str() );
    }
  }
  m_hall = SpaceGroupHallNumber{ hn };
}

unsigned NC::SpaceGroup::number() const noexcept
{
  return entry( m_hall.get() ).number;
}

const char * NC::SpaceGroup::choice() const noexcept
{
  return entry( m_hall.get() ).choice;
}

const char * NC::SpaceGroup::hallSymbol() const noexcept
{
  return entry( m_hall.get() ).hall;
}

bool NC::SpaceGroup::isDefaultSetting() const noexcept
{
  const auto hn = m_hall.get();
  return hn == 1 || entry( hn - 1 ).number != entry( hn ).number;
}

bool NC::SpaceGroup::hasMultipleSettings() const noexcept
{
  const auto hn = m_hall.get();
  const auto n = entry( hn ).number;
  return ( ( hn > 1 && entry( hn - 1 ).number == n )
           || ( hn < SpaceGroupHallNumber::max_value
                && entry( hn + 1 ).number == n ) );
}

std::string NC::SpaceGroup::toString() const
{
  //NB: Only default settings have empty choice codes (e.g. orthorhombic
  //settings with standard axes):
  std::string res = std::to_string( number() );
  if ( *choice() ) {//i.e. choice code is not the empty string
    res += ':';
    res += choice();
  }
  return res;
}

std::ostream& NC::operator<<( std::ostream& os, const SpaceGroup& sg )
{
  return os << sg.toString();
}

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

NC::SGCrystalSystem NC::SpaceGroup::crystalSystem() const noexcept
{
  const unsigned n = number();
  if ( n <= 2 )
    return SGCrystalSystem::Triclinic;
  if ( n <= 15 )
    return SGCrystalSystem::Monoclinic;
  if ( n <= 74 )
    return SGCrystalSystem::Orthorhombic;
  if ( n <= 142 )
    return SGCrystalSystem::Tetragonal;
  if ( n <= 167 )
    return SGCrystalSystem::Trigonal;
  if ( n <= 194 )
    return SGCrystalSystem::Hexagonal;
  return SGCrystalSystem::Cubic;
}

NC::SGSettingInfo NC::SpaceGroup::settingInfo() const noexcept
{
  SGSettingInfo res{ SGHexFamilyAxes::NotApplicable,
                     SGUniqueAxis::NotApplicable,
                     SGCellChoice::NotApplicable,
                     SGOriginChoice::NotApplicable,
                     SGAxisPermutation::NotApplicable };
  //NB: The choice codes come from the hardwired table (fully verified by
  //tests), so anything unexpected would be a bug (hence nc_assert, since this
  //function is noexcept):
  const std::string c = choice();
  auto originFromDigit = []( char ch )
  {
    nc_assert( ch == '1' || ch == '2' );
    return ch == '1' ? SGOriginChoice::Choice1 : SGOriginChoice::Choice2;
  };
  switch ( crystalSystem() ) {
  case SGCrystalSystem::Triclinic:
    nc_assert( c.empty() );
    break;
  case SGCrystalSystem::Monoclinic:
    {
      //Codes like "b", "b1", "-c2":
      std::size_t i = 0;
      const bool minus = ( !c.empty() && c[0] == '-' );
      if ( minus )
        ++i;
      nc_assert( i < c.size() );
      const char ax = c[i++];
      nc_assert( ax == 'a' || ax == 'b' || ax == 'c' );
      res.uniqueAxis = ( ax == 'a'
                         ? ( minus ? SGUniqueAxis::minus_a : SGUniqueAxis::a )
                         : ax == 'b'
                         ? ( minus ? SGUniqueAxis::minus_b : SGUniqueAxis::b )
                         : ( minus ? SGUniqueAxis::minus_c : SGUniqueAxis::c ) );
      if ( i < c.size() ) {
        nc_assert( i + 1 == c.size() );
        const char cc = c[i];
        nc_assert( cc >= '1' && cc <= '3' );
        res.cellChoice = ( cc == '1' ? SGCellChoice::Choice1
                           : cc == '2' ? SGCellChoice::Choice2
                           : SGCellChoice::Choice3 );
      }
    }
    break;
  case SGCrystalSystem::Orthorhombic:
    {
      //Codes like "", "cab", "1", "2ba-c":
      std::string perm = c;
      if ( !perm.empty() && ( perm[0] == '1' || perm[0] == '2' ) ) {
        res.originChoice = originFromDigit( perm[0] );
        perm = perm.substr( 1 );
      }
      if ( perm.empty() || perm == "abc" )
        res.axisPermutation = SGAxisPermutation::abc;
      else if ( perm == "ba-c" )
        res.axisPermutation = SGAxisPermutation::ba_mc;
      else if ( perm == "cab" )
        res.axisPermutation = SGAxisPermutation::cab;
      else if ( perm == "-cba" )
        res.axisPermutation = SGAxisPermutation::mcba;
      else if ( perm == "bca" )
        res.axisPermutation = SGAxisPermutation::bca;
      else if ( perm == "a-cb" )
        res.axisPermutation = SGAxisPermutation::a_mcb;
      else
        nc_assert( false );
    }
    break;
  case SGCrystalSystem::Tetragonal:
  case SGCrystalSystem::Cubic:
    if ( !c.empty() ) {
      nc_assert( c.size() == 1 );
      res.originChoice = originFromDigit( c[0] );
    }
    break;
  case SGCrystalSystem::Trigonal:
  case SGCrystalSystem::Hexagonal:
    nc_assert( c.empty() || c == "H" || c == "R" );
    res.hexFamilyAxes = ( c == "R" ? SGHexFamilyAxes::Rhombohedral
                          : SGHexFamilyAxes::Hexagonal );
    break;
  }
  return res;
}

std::ostream& NC::operator<<( std::ostream& os, SGCrystalSystem cs )
{
  switch ( cs ) {
  case SGCrystalSystem::Triclinic: return os << "triclinic";
  case SGCrystalSystem::Monoclinic: return os << "monoclinic";
  case SGCrystalSystem::Orthorhombic: return os << "orthorhombic";
  case SGCrystalSystem::Tetragonal: return os << "tetragonal";
  case SGCrystalSystem::Trigonal: return os << "trigonal";
  case SGCrystalSystem::Hexagonal: return os << "hexagonal";
  case SGCrystalSystem::Cubic: return os << "cubic";
  }
  return os;
}

std::ostream& NC::operator<<( std::ostream& os, SGHexFamilyAxes v )
{
  switch ( v ) {
  case SGHexFamilyAxes::NotApplicable: return os << "n/a";
  case SGHexFamilyAxes::Hexagonal: return os << "hexagonal";
  case SGHexFamilyAxes::Rhombohedral: return os << "rhombohedral";
  }
  return os;
}

std::ostream& NC::operator<<( std::ostream& os, SGUniqueAxis v )
{
  switch ( v ) {
  case SGUniqueAxis::NotApplicable: return os << "n/a";
  case SGUniqueAxis::a: return os << "a";
  case SGUniqueAxis::b: return os << "b";
  case SGUniqueAxis::c: return os << "c";
  case SGUniqueAxis::minus_a: return os << "-a";
  case SGUniqueAxis::minus_b: return os << "-b";
  case SGUniqueAxis::minus_c: return os << "-c";
  }
  return os;
}

std::ostream& NC::operator<<( std::ostream& os, SGCellChoice v )
{
  switch ( v ) {
  case SGCellChoice::NotApplicable: return os << "n/a";
  case SGCellChoice::Choice1: return os << "1";
  case SGCellChoice::Choice2: return os << "2";
  case SGCellChoice::Choice3: return os << "3";
  }
  return os;
}

std::ostream& NC::operator<<( std::ostream& os, SGOriginChoice v )
{
  switch ( v ) {
  case SGOriginChoice::NotApplicable: return os << "n/a";
  case SGOriginChoice::Choice1: return os << "1";
  case SGOriginChoice::Choice2: return os << "2";
  }
  return os;
}

std::ostream& NC::operator<<( std::ostream& os, SGAxisPermutation v )
{
  switch ( v ) {
  case SGAxisPermutation::NotApplicable: return os << "n/a";
  case SGAxisPermutation::abc: return os << "abc";
  case SGAxisPermutation::ba_mc: return os << "ba-c";
  case SGAxisPermutation::cab: return os << "cab";
  case SGAxisPermutation::mcba: return os << "-cba";
  case SGAxisPermutation::bca: return os << "bca";
  case SGAxisPermutation::a_mcb: return os << "a-cb";
  }
  return os;
}

std::ostream& NC::operator<<( std::ostream& os, const SGSettingInfo& si )
{
  return os << "{hexFamilyAxes=" << si.hexFamilyAxes
            << ", uniqueAxis=" << si.uniqueAxis
            << ", cellChoice=" << si.cellChoice
            << ", originChoice=" << si.originChoice
            << ", axisPermutation=" << si.axisPermutation << "}";
}

#ifndef NCrystal_SpaceGroup_hh
#define NCrystal_SpaceGroup_hh

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

#include "NCrystal/internal/sgsym/NCSpaceGroupHallNumber.hh"
#include "NCrystal/internal/utils/NCStrView.hh"

namespace NCRYSTAL_NAMESPACE {

  //Immutable class representing one of the 530 space group settings of ITVB
  //2001 Table A1.4.2.7, i.e. a space group number (1..230) along with a
  //setting (e.g. an origin choice, a monoclinic unique axis and cell choice,
  //an orthorhombic axis permutation, or hexagonal/rhombohedral axes).
  //
  //String representations are "<number>:<choice>", using the choice codes of
  //ITVB (e.g. "227:2", "166:R", "14:c1", "62:cab"), or just "<number>" when
  //the choice code is empty (numbers with a single setting, and the default
  //setting of orthorhombic numbers). When constructing from a string, a bare
  //number selects the default (i.e. first listed) setting.

  class SpaceGroup final {
  public:
    //Constructors (throw BadInput in case of invalid input). The string
    //constructors accept StrView, std::string (incl. temporaries, which a
    //StrView can not refer to), and C-strings (incl. literals):
    explicit SpaceGroup( SpaceGroupHallNumber );
    explicit SpaceGroup( StrView );
    explicit SpaceGroup( const std::string& s ) : SpaceGroup( StrView( s ) ) {}
    explicit SpaceGroup( const char * s ) : SpaceGroup( StrView( s ) ) {}

    SpaceGroupHallNumber hallNumber() const noexcept { return m_hall; }
    unsigned number() const noexcept;//1..230
    const char * choice() const noexcept;//e.g. "2", "R", "c1", or "" (never
                                         //nullptr)
    const char * hallSymbol() const noexcept;//e.g. "-F 4vw 2vw 3"
    bool isDefaultSetting() const noexcept;//first listed setting of number?
    bool hasMultipleSettings() const noexcept;//number has other settings?
    std::string toString() const;//e.g. "227:2", "62:cab", "62", or "225"

    bool operator==( const SpaceGroup& o ) const noexcept;
    bool operator!=( const SpaceGroup& o ) const noexcept;
    bool operator<( const SpaceGroup& o ) const noexcept;//Hall number order
  private:
    SpaceGroupHallNumber m_hall;
  };

  std::ostream& operator<<( std::ostream&, const SpaceGroup& );

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline bool SpaceGroup::operator==( const SpaceGroup& o ) const noexcept
  {
    return m_hall == o.m_hall;
  }
  inline bool SpaceGroup::operator!=( const SpaceGroup& o ) const noexcept
  {
    return m_hall != o.m_hall;
  }
  inline bool SpaceGroup::operator<( const SpaceGroup& o ) const noexcept
  {
    return m_hall < o.m_hall;
  }
}

#endif

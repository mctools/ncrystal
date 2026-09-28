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

//TODO: Migrate crystalSystem(), checkAndCompleteLattice[Angles](),
//isRhombohedralSpaceGroup() and usesRhombohedralAxes() from NCLatticeUtils.hh
//to the sgsym component (where they naturally belong, cf. SGCrystalSystem and
//SGCellConstraints). The CrystalSystem enum name can then be reused.

namespace NCRYSTAL_NAMESPACE {

  //Crystal system of a space group (NB: named SGCrystalSystem, since the name
  //CrystalSystem is currently taken by an enum in NCLatticeUtils.hh):
  enum class SGCrystalSystem { Triclinic, Monoclinic, Orthorhombic, Tetragonal,
                               Trigonal, Hexagonal, Cubic };

  //Description of the setting of a space group, derived from the ITVB choice
  //code. Fields which do not apply to a given space group are NotApplicable.

  //Axes of the hexagonal crystal family (i.e. trigonal and hexagonal crystal
  //systems). The R space groups can be described in either (Hexagonal is the
  //standard setting), all other space groups of the family always use
  //Hexagonal:
  enum class SGHexFamilyAxes { NotApplicable, Hexagonal, Rhombohedral };

  //Unique axis of monoclinic space groups. The standard settings use b, while
  //the other values are only for non-standard settings. The minus_ values
  //correspond to ITVB choice codes like "-b" (same unique axis, but with the
  //two other axes interchanged, i.e. a reversed cell orientation):
  enum class SGUniqueAxis { NotApplicable, a, b, c, minus_a, minus_b, minus_c };

  //Cell choice of monoclinic space groups with centring or glide planes. The
  //standard settings use Choice1, while Choice2 and Choice3 are only for
  //non-standard settings:
  enum class SGCellChoice { NotApplicable, Choice1, Choice2, Choice3 };

  //Origin choice of the 24 space groups which have two (both are standard
  //settings in ITVB):
  enum class SGOriginChoice { NotApplicable, Choice1, Choice2 };

  //Axis permutations of orthorhombic space groups (ITVB notation abc, ba-c,
  //cab, -cba, bca, a-cb, where "m" here stands for the minus sign). The
  //standard settings use abc, while the other values are only for
  //non-standard settings:
  enum class SGAxisPermutation { NotApplicable, abc, ba_mc, cab, mcba, bca,
                                 a_mcb };

  struct SGSettingInfo {
    SGHexFamilyAxes hexFamilyAxes;
    SGUniqueAxis uniqueAxis;
    SGCellChoice cellChoice;
    SGOriginChoice originChoice;
    SGAxisPermutation axisPermutation;
  };

  //Output (using ITVB spellings, e.g. "ba-c" or "-b", and "n/a" for
  //NotApplicable):
  std::ostream& operator<<( std::ostream&, SGCrystalSystem );
  std::ostream& operator<<( std::ostream&, SGHexFamilyAxes );
  std::ostream& operator<<( std::ostream&, SGUniqueAxis );
  std::ostream& operator<<( std::ostream&, SGCellChoice );
  std::ostream& operator<<( std::ostream&, SGOriginChoice );
  std::ostream& operator<<( std::ostream&, SGAxisPermutation );
  std::ostream& operator<<( std::ostream&, const SGSettingInfo& );

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

    SGCrystalSystem crystalSystem() const noexcept;
    SGSettingInfo settingInfo() const noexcept;

    //NB: Symmetry operations etc. are available via SGSymmetry::get(..).

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

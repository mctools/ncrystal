#ifndef NCrystal_EqRefl_hh
#define NCrystal_EqRefl_hh

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

#include "NCrystal/core/NCTypes.hh"

//Class EqRefl provides symmetry-equivalent reflections for a given space group
//number, by providing a list of all (h,k,l) indices symmetry-equivalent to a
//given (h,k,l) index. Or more precisely, half of the list, since the symmetry
//(h,k,l)->(-h,-k,-l) is trivial and would be a waste of memory to expand
//explicitly.
//
//Note that while two equivalent (h,k,l) indices will always have identical
//values of d-spacing and structure factors (F^2), the opposite is not in
//general true. So (h,k,l) groupings based on (d,F^2) values might contain more
//than one symmetry-equivalency group in each (d,F^2) group.
//
//TODO: SCBragg currently assumes Info objects to either have
//demi-normals or contain equivalency groupings. This is currently true, but
//seems a bit fragile (e.g. the NCMAT factory produces (d,F^2) groups, not
//equivalency groups - but does fortunately provide demi-normals). A more robust
//situation would be desirable (either all factories should produce hkl lists
//grouped by equivalency, or SCBragg code should become more intelligent).

namespace NCRYSTAL_NAMESPACE {

  class EqRefl
  {
  public:

    //Initialise with space group number in 1..230:
    EqRefl( int spacegroup );

    class EquivReflList {
      //List with half the HKL points in a given family of symmetric equivalent
      //HKL points (only half, since (h,k,l) and (-h,-k,-l) are trivially in the
      //same group). This uses non-dynamic storage to keep malloc overhead low,
      //and thus takes up roughly 100B of storage (it is thus mostly suitable
      //for short-lived usage on the stack).
    public:
      const HKL* begin() const noexcept;
      const HKL& front() const noexcept { return *begin(); }
      const HKL* end() const noexcept;
      std::size_t size() const noexcept;
      bool contains( HKL ) const noexcept;//is (h,k,l) or (-h,-k,-l) in the list?

      //Empty (for later assignment):
      EquivReflList();

      //Copy only:
      EquivReflList( const EquivReflList& );
      EquivReflList& operator=( const EquivReflList& );

    private:
      HKL m_data[24];
      HKL * m_end;
      friend class EqRefl;
      void add(int, int, int);
    };

    EquivReflList getEquivalentReflections(int h, int k, int l) const;
    EquivReflList getEquivalentReflections( const HKL& hkl ) const
    {
      return getEquivalentReflections( hkl.h, hkl.k, hkl.l );
    }

    //Same as getEquivalentReflections(..).front() but slightly faster (for when
    //only that value is needed):
    HKL getEquivalentReflectionsRepresentativeValue( const HKL& hkl ) const;
    HKL getEquivalentReflectionsRepresentativeValue( int h, int k, int l ) const;

  private:
    struct Helper;
    friend struct Helper;
    EquivReflList (*m_calc) (int,int,int) = nullptr;
  };
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  inline const HKL* EqRefl::EquivReflList::begin() const noexcept { return &m_data[0]; }
  inline const HKL* EqRefl::EquivReflList::end() const noexcept { return m_end; }
  inline std::size_t EqRefl::EquivReflList::size() const noexcept { return std::distance(begin(),end()); }
  inline EqRefl::EquivReflList::EquivReflList() { m_end = &m_data[0]; }

  inline EqRefl::EquivReflList::EquivReflList( const EquivReflList& o )
  {
    m_end = &m_data[0];
    for ( auto& e : o )
      *m_end++ = e;
  }

  inline EqRefl::EquivReflList& EqRefl::EquivReflList::operator=( const EquivReflList& o )
  {
    m_end = &m_data[0];
    for ( auto& e : o )
      *m_end++ = e;
    return *this;
  }

  inline void EqRefl::EquivReflList::add(int h, int k, int l)
  {
    //Only ever consider the form of (h,k,l) and (-h,-k,-l) which has first
    //non-zero coordinate positive:
    HKL a(h,k,l);
    auto am = a.flipped();
    *m_end++ =  am < a ? am : a;
  }

  inline EqRefl::EquivReflList EqRefl::getEquivalentReflections(int h, int k, int l) const
  {
    auto raw = m_calc(h,k,l);
    //Sort and discard duplicates:
    std::sort((HKL*)&raw.m_data[0],raw.m_end);
    raw.m_end = std::unique((HKL*)&raw.m_data[0],raw.m_end);
    return raw;
  }

  inline HKL EqRefl::getEquivalentReflectionsRepresentativeValue( int h, int k, int l ) const
  {
    auto raw = m_calc(h,k,l);
    auto itLow = raw.begin();
    auto itE = raw.end();
    for ( auto it = std::next(itLow); it!=itE; ++it ) {
      if ( *it < *itLow )
        itLow = it;
    }
    return *itLow;
  }

  inline HKL EqRefl::getEquivalentReflectionsRepresentativeValue( const HKL& e ) const
  {
    return getEquivalentReflectionsRepresentativeValue( e.h, e.k, e.l );
  }

  inline bool EqRefl::EquivReflList::contains( HKL a ) const noexcept
  {
    //At worst 24 entries, often just a few handfuls => use simple linear search.
    HKL am(-a.h,-a.k,-a.l);
    if ( am < a )
      a = am;
    for ( auto& e : *this ) {
      if ( e == a )
        return true;
    }
    return false;
  }

}

#endif

#ifndef NCrystal_CompositionUtils_hh
#define NCrystal_CompositionUtils_hh

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

#include "NCrystal/interfaces/NCInfo.hh"

namespace NCRYSTAL_NAMESPACE {

  //Utilities for breaking down the basic composition of a material into
  //elements and isotopes. If an element only occurs as a natural element (A=0)
  //and has no specific isotopes, that element will be returned as a natural
  //isotope unless the ForceIsotopes flag is chosen. In any case, it is best if
  //a NaturalAbundanceProvider is provided, since if a given Z value has a mix
  //of isotopes and the natural elements, the natural element must always be
  //broken up. This behaviour is chosen in order to facilitate the needs of the
  //Geant4 interface (where we want to avoid two G4Elements's with the same Z
  //value), and might be further customised for other interfaces if needed at
  //some point. Note however, that it is allowed to not provide a
  //NaturalAbundanceProvider (i.e. set it to nullptr), which will then only give
  //an exception if it is actually needed.
  //
  //Both full and light-weight breakdowns are available. The former rely solely
  //on standard containers, while the latter use custom classes with a lower
  //memory footprint.

  namespace CompositionUtils {

    typedef std::pair<unsigned,std::vector<std::pair<unsigned,double>>> FullElementBreakdown;//( Z, [(A,fraction),...] )
    typedef std::function< std::vector<std::pair<unsigned,double>> ( unsigned ) > NaturalAbundanceProvider;//Z -> [(A,fraction),...]
    typedef std::vector<FullElementBreakdown> FullBreakdown;
    enum ForceIsotopesChoice { ForceIsotopes, PreferNaturalElements };
    NCRYSTAL_API FullBreakdown createFullBreakdown( const Info::Composition&,
                                                    const NaturalAbundanceProvider&,
                                                    ForceIsotopesChoice = PreferNaturalElements );


    NCRYSTAL_API std::string fullBreakdownToJSON( const FullBreakdown& );

    class NCRYSTAL_API ElementBreakdownLW : private MoveOnly {
      //Struct for keeping isotope breakdown of a particular element. Care is taken
      //to keep memory usage low and to have it easily sortable.
    public:
      ElementBreakdownLW(const FullElementBreakdown&);
      unsigned Z() const { return m_ZAN >> 24; }
      unsigned nIsotopes() const { return m_ZAN & 0x3FFF; }//0 for natural elements
      bool isNaturalElement() const { return firstA()==0; }
      //A(idx) requires idx < nIsotopes (thus it can not be called for natural elements)
      unsigned A(unsigned idx) const { nc_assert(idx<nIsotopes()&&(idx==0||m_other!=nullptr)); return idx ? m_other[idx-1].second : firstA(); }
      bool operator<(const ElementBreakdownLW& o) const;
      double fraction(unsigned idx) const;
      std::string description(unsigned precision=6) const;

      ElementBreakdownLW( ElementBreakdownLW&& o );
      ElementBreakdownLW& operator=( ElementBreakdownLW&& o );
      bool valid() const;//objects are always valid, unless has been moved-from.
    private:
      uint32_t m_ZAN = 0;//lowest 14 bits: nisotopes. Middle 10 bits: A. Highest 8 bits: Z.
      std::unique_ptr<std::pair<double,uint16_t>[]> m_other;//Other isotopes and their fractions.
      bool cmpOthers(const ElementBreakdownLW&) const;
      double calcFirstFraction() const;
      unsigned firstA() const { return ( m_ZAN >> 14 ) & 0x3FF; }//0 for natural elements
    };

    typedef std::vector<std::pair<double,ElementBreakdownLW>> LWBreakdown;//(fraction, element)
    NCRYSTAL_API std::string breakdownToStr( const LWBreakdown&, unsigned precision=6 );
    NCRYSTAL_API LWBreakdown createLWBreakdown( const Info::Composition& a,
                                                const NaturalAbundanceProvider& b,
                                                ForceIsotopesChoice c = PreferNaturalElements );
  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace CompositionUtils {
    inline bool ElementBreakdownLW::operator<(const ElementBreakdownLW& o) const {
      //TK: I am not 100% sure if it a bug or not, but it seems that when using
      //LWBreakdown objects as keys in a std::map with gcc 7.5.0, this
      //comparison operator gets called on objects that have already been
      //moved-from. Thus it might be that we are dealing with invalid objects
      //here, which is why we had to reimplement the move operators below to
      //ensure all moved-from objects have m_ZAN==0. The comparison below should
      //be safe and deterministic for invalid objects as well.
      //
      //nc_assert(valid()); nc_assert(o.valid());
      //
      //First sort by ZAN field:
      if ( m_ZAN != o.m_ZAN )
        return m_ZAN < o.m_ZAN;
      //Ok, ZAN fields are identical, so must have same number of isotopes.
      //Compare compositions:
      return cmpOthers(o);
    }
    inline double ElementBreakdownLW::fraction(unsigned idx) const {
      nc_assert( valid() );
      const unsigned nc = nIsotopes();
      nc_assert( idx < nc );
      if ( idx > 0 ) {
        nc_assert(m_other!=nullptr);
        return m_other[idx-1].first;
      }
      return nc==1 ? 1.0 : calcFirstFraction();
    }
    inline ElementBreakdownLW::ElementBreakdownLW( ElementBreakdownLW&& o )
      : m_ZAN(o.m_ZAN), m_other(std::move(o.m_other))
    {
      //Make sure even a moved-from object fullfills bool(m_other)==bool(N()>1)
      o.m_ZAN = 0;
      o.m_other.reset();
      nc_assert(!o.valid());
      nc_assert(valid());
    }
    inline ElementBreakdownLW& ElementBreakdownLW::operator=( ElementBreakdownLW&& o )
    {
      nc_assert(o.valid());
      m_ZAN = o.m_ZAN;
      m_other = std::move(o.m_other);
      //Make sure even a moved-from object fullfills bool(m_other)==bool(N()>1)
      o.m_ZAN = 0;
      o.m_other.reset();
      nc_assert(!o.valid());
      nc_assert(valid());
      return *this;
    }
    inline bool ElementBreakdownLW::valid() const
    {
      return m_ZAN!=0;
    }
  }
}

#endif

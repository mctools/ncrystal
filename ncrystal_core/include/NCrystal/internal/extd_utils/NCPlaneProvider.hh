#ifndef NCrystal_PlaneProvider_hh
#define NCrystal_PlaneProvider_hh

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

#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/phys_utils/NCEqRefl.hh"
#include "NCrystal/internal/utils/NCSpan.hh"
#include "NCrystal/interfaces/NCInfoTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  class Info;

  //Utilities for extracting expanded lists of plane normals or Miller (hkl)
  //indices from HKLLists, usign respectively the PlaneProvider or the
  //ExpandHKLHelper classes. Note that in both cases, only half of the entries
  //are returned, since (h,k,l)->(-h,-k,-l) will always have similar fsquared,
  //dspacing and a plane normal which is simply a single flip.

  class PlaneProvider {
  public:

    //Helper class providing loops over crystal (demi) planes, each loop giving
    //dspacing, F-squared and a demi normal.

    PlaneProvider();
    virtual ~PlaneProvider();

    //Method used to loop over planes (returns NullOpt when iteration is done):
    struct Plane {
      double dspacing, fsq;
      Vector demi_normal;
    };
    virtual Optional<Plane> getNextPlane() = 0;

    //Rewind the looping to prepare for a new loop with getNextPlane (does not
    //need to be called for the first loop):
    virtual void prepareLoop() = 0;

    //Whether or not it is safe to call getNextPlane and prepareLoop. If calling
    //anyway, exceptions will be thrown by conforming implementations. A false
    //return value here usually indicates incomplete information for normals to
    //be provided:
    virtual bool canProvide() const = 0;
  };

  //Creates standard plane provider from Info object, which will attempt various
  //means of producing the HKL normals (preferring actual deminormals if
  //available, then expanded hkl info and finally falling back to attempting
  //their construction based on space group and multiplicity info):

  //Version partaking in lifetime management of Info:
  std::unique_ptr<PlaneProvider> createStdPlaneProvider( shared_obj<const Info> );

  //Version in which caller guarantees Info object will remain alive as long
  //as any plane provider methods are called:
  std::unique_ptr<PlaneProvider> createStdPlaneProvider( const Info* );

  class ExpandHKLHelper {
  public:

    //Helper class for extracting HKL list from a given HKLInfo object:

    ExpandHKLHelper( const Info& );//convenience, simply digs out spacegroup
    ExpandHKLHelper( int spacegroup );

    //Use the following function to iterate over all (or rather, half of) HKL
    //entries in a given HKLInfo object. The returned span is invalidated once
    //the fullHKLList method is called again or the HKLInfo objects or
    //ExpandHKLHelper itself are modified:
    Span<const HKL> expand( const HKLInfo& );

    bool canExpand( HKLInfoType ) const;

  private:
    Optional<EqRefl> m_sym;
    Optional<EqRefl::EquivReflList> m_list;
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline ExpandHKLHelper::ExpandHKLHelper( int spacegroup )
  {
    nc_assert( spacegroup>=0 && spacegroup <= 230 );
    if ( spacegroup )
      m_sym.emplace(spacegroup);
  }

  inline Span<const HKL> ExpandHKLHelper::expand( const HKLInfo& hi )
  {
    switch ( hi.type() ) {
    case HKLInfoType::SymEqvGroup:
      nc_assert(m_sym.has_value()&&"ExpandHKLHelper not constructed with"
                " spacegroup and requested to expand an SymEqvGroup entry");
      m_list = m_sym.value().getEquivalentReflections( hi.hkl );
      nc_assert(hi.multiplicity == m_list.value().size()*2);
      return { m_list.value().begin(), m_list.value().end() };
    case HKLInfoType::ExplicitHKLs:
      nc_assert( hi.explicitValues->list.has_value<std::vector<HKL>>() );
      return hi.explicitValues->list.get<std::vector<HKL>>();
    case HKLInfoType::ExplicitNormals:
    case HKLInfoType::Minimal:
      break;
    };
    nc_assert_always(false);
    return {nullptr,nullptr};//should not happen
  }

  inline bool ExpandHKLHelper::canExpand( HKLInfoType hitype ) const
  {
    switch ( hitype ) {
    case HKLInfoType::SymEqvGroup:
      return m_sym.has_value();
    case HKLInfoType::ExplicitHKLs:
      return true;
    case HKLInfoType::ExplicitNormals:
    case HKLInfoType::Minimal:
      break;
    };
    return false;
  }

}

#endif

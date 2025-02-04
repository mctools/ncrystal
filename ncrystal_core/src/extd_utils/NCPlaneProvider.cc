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

#include "NCrystal/internal/extd_utils/NCPlaneProvider.hh"
#include "NCrystal/internal/extd_utils/NCOrientUtils.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"
#include "NCrystal/internal/phys_utils/NCEqRefl.hh"

namespace NC = NCrystal;

NC::PlaneProvider::PlaneProvider() = default;
NC::PlaneProvider::~PlaneProvider() = default;

namespace NCRYSTAL_NAMESPACE {

  namespace {

    class PlaneProviderStd_Unable final  : public PlaneProvider {
    public:
      bool canProvide() const override { return false; }
      void prepareLoop() override
      {
        NCRYSTAL_THROW(MissingInfo,"Insufficient information to provide reflection plane normals.");
      }
      Optional<Plane> getNextPlane() override
      {
        NCRYSTAL_THROW(LogicError,"Do not call getNextPlane() without first checking canProvide() and calling prepareLoop().");
        return NullOpt;
      }
    };

    class PlaneProviderStd_AbleButEmpty final  : public PlaneProvider {
    public:
      //Special case, we can always provide all the normals of an empty hkl list
      //- there simply won't be any.
      bool canProvide() const override { return true; }
      void prepareLoop() override {}
      Optional<Plane> getNextPlane() override { return NullOpt; }
    };

    class PlaneProviderStd_Normals final : public PlaneProvider {
      OptionalInfoPtr m_strongRef;
      double m_dsp, m_fsq;
      HKLList::const_iterator m_it, m_itB, m_itE;
      std::vector<HKLInfo::Normal>::const_iterator m_it_inner, m_it_innerE;
    public:
      PlaneProviderStd_Normals( const Info * info, OptionalInfoPtr iptr )
        : PlaneProvider(), m_strongRef(iptr)
      {
        nc_assert_always( info );
        nc_assert_always( info->hasHKLInfo() );
        nc_assert_always( info->hklInfoType() == HKLInfoType::ExplicitNormals );
        auto& l = info->hklList();
        m_it = m_itB = l.begin();
        m_itE = l.end();
        prepareLoop();
      }

      void prepareLoopInner()
      {
        if ( m_it == m_itE )
          return;
        nc_assert( m_it->explicitValues != nullptr );
        nc_assert(m_it->explicitValues->list.has_value<std::vector<HKLInfo::Normal>>());
        auto& v = m_it->explicitValues->list.get<std::vector<HKLInfo::Normal>>();
        m_it_inner = v.begin();
        m_it_innerE = v.end();
        m_dsp = m_it->dspacing;
        m_fsq = m_it->fsquared;
      }

      bool canProvide() const override { return true; }

      void prepareLoop() override
      {
        m_it = m_itB;
        prepareLoopInner();
      }

      Optional<Plane> getNextPlane() override
      {
        if ( m_it_inner == m_it_innerE ) {
          if ( ++m_it == m_itE )
            return NullOpt;
          prepareLoopInner();
          return getNextPlane();
        }
        return Plane{ m_dsp, m_fsq, (m_it_inner++)->as<Vector>() };
      }
    };


    class PlaneProviderStd_HKL final : public PlaneProvider {
      OptionalInfoPtr m_strongRef;
      double m_dsp, m_fsq;
      ExpandHKLHelper m_hklExpander;
      RotMatrix m_reci_lattice;
      HKLList::const_iterator m_it, m_itB, m_itE;
      const HKL * m_it_inner;
      const HKL * m_it_innerE;
    public:
      PlaneProviderStd_HKL( const Info * info, OptionalInfoPtr iptr )
        : PlaneProvider(),
          m_strongRef(iptr),
          m_hklExpander( [&info](){
            nc_assert(info->hasStructureInfo());
            nc_assert_always( info );
            nc_assert_always( info->hasHKLInfo() );
            nc_assert_always( isOneOf(info->hklInfoType(),HKLInfoType::SymEqvGroup,HKLInfoType::ExplicitHKLs) );
            return info->getStructureInfo().spacegroup;
          }() ),
          m_reci_lattice( getReciprocalLatticeRot( info->getStructureInfo() ) )
      {
        nc_assert( m_hklExpander.canExpand( info->hklInfoType() ) );
        auto& l = info->hklList();
        m_it = m_itB = l.begin();
        m_itE = l.end();
        prepareLoop();
      }

      void prepareLoopInner()
      {
        if ( m_it == m_itE )
          return;
        nc_assert( isOneOf( m_it->type(), HKLInfoType::SymEqvGroup, HKLInfoType::ExplicitHKLs) );
        auto v = m_hklExpander.expand( *m_it );
        m_it_inner = v.begin();
        m_it_innerE = v.end();
        m_dsp = m_it->dspacing;
        m_fsq = m_it->fsquared;
      }

      bool canProvide() const override { return true; }

      void prepareLoop() override
      {
        m_it = m_itB;
        prepareLoopInner();
      }

      Optional<Plane> getNextPlane() override
      {
        if ( m_it_inner == m_it_innerE ) {
          if ( ++m_it == m_itE )
            return NullOpt;
          prepareLoopInner();
          return getNextPlane();
        }
        Plane p{ m_dsp,
                 m_fsq,
                 m_reci_lattice * Vector( m_it_inner->h,
                                          m_it_inner->k,
                                          m_it_inner->l ) };
        p.demi_normal.normalise();
        ++m_it_inner;
        return p;
      }
    };

    std::unique_ptr<PlaneProvider> actual_createStdPlaneProvider( const Info* info, OptionalInfoPtr iptr )
    {
      auto unable = []() { return std::make_unique<PlaneProviderStd_Unable>(); };

      if ( !info->hasHKLInfo() )
        return unable();

      if ( info->hklList().empty() ) {
        //special case, no matter the hkl info type it is always possible to
        //turn an empty hkl list into a valid but empty list of normals! This is
        //particularly important to avoid a spurious failure, since we (at least
        //for now?) always use a type of HKLInfoType::Minimal when the hkl list
        //is valid but empty.
        return std::make_unique<PlaneProviderStd_AbleButEmpty>();
      }

      auto hitype = info->hklInfoType();
      switch( hitype ) {
      case HKLInfoType::SymEqvGroup:
        if ( !info->hasStructureInfo() || info->getStructureInfo().spacegroup == 0 )
          return unable();
        return std::make_unique<PlaneProviderStd_HKL>( info, std::move(iptr) );
      case HKLInfoType::ExplicitHKLs:
        if ( !info->hasStructureInfo() )
          return unable();
        return std::make_unique<PlaneProviderStd_HKL>( info, std::move(iptr) );
      case HKLInfoType::ExplicitNormals:
        return std::make_unique<PlaneProviderStd_Normals>( info, std::move(iptr) );
      case HKLInfoType::Minimal:
        return unable();
      };
      return unable();
    }
  }
}

std::unique_ptr<NC::PlaneProvider> NC::createStdPlaneProvider( InfoPtr info)
{
  auto rawinfo = info.get();
  return actual_createStdPlaneProvider( rawinfo, std::move(info) );
}

std::unique_ptr<NC::PlaneProvider> NC::createStdPlaneProvider(const Info* info)
{
  nc_assert(info!=nullptr);
  return actual_createStdPlaneProvider( info, nullptr );
}

NC::ExpandHKLHelper::ExpandHKLHelper( const Info& info )
  : ExpandHKLHelper( info.hasStructureInfo() ? info.getStructureInfo().spacegroup : 0 )
{
}

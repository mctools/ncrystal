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

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  //As good a place as any for all of these SmallVector asserts (since we do not
  //have NCSmallVector.cc):
  static_assert(SmallVector_IC<int,1>::mode==SVMode::FASTACCESS,"");
  static_assert(SmallVector_IC<int,1,SVMode::FASTACCESS>::mode==SVMode::FASTACCESS,"");
  static_assert(SmallVector_IC<int,1,SVMode::LOWFOOTPRINT>::mode==SVMode::LOWFOOTPRINT,"");
  static_assert(SmallVector<int,1>::mode==SVMode::FASTACCESS,"");
  static_assert(SmallVector<int,1,SVMode::FASTACCESS>::mode==SVMode::FASTACCESS,"");
  static_assert(SmallVector<int,1,SVMode::LOWFOOTPRINT>::mode==SVMode::LOWFOOTPRINT,"");
  static_assert(std::is_copy_constructible<SmallVector_IC<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(std::is_copy_constructible<SmallVector_IC<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(std::is_copy_assignable<SmallVector_IC<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(std::is_copy_assignable<SmallVector_IC<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(!std::is_copy_constructible<SmallVector<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(!std::is_copy_constructible<SmallVector<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(!std::is_copy_assignable<SmallVector<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(!std::is_copy_assignable<SmallVector<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(std::is_move_constructible<SmallVector_IC<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(std::is_move_constructible<SmallVector_IC<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(std::is_move_assignable<SmallVector_IC<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(std::is_move_assignable<SmallVector_IC<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(std::is_move_constructible<SmallVector<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(std::is_move_constructible<SmallVector<int,1,SVMode::LOWFOOTPRINT>>::value,"");
  static_assert(std::is_move_assignable<SmallVector<int,1,SVMode::FASTACCESS>>::value,"");
  static_assert(std::is_move_assignable<SmallVector<int,1,SVMode::LOWFOOTPRINT>>::value,"");
}

NC::DataSourceName::DataSourceName()
  : m_str( [](){ static auto s_def = makeSO<std::string>(); return s_def; }() )
{
  static_assert( EnergyDomain::null().isNull(), "" );
}

std::ostream& NC::operator<<( std::ostream& os, const OrientDir& od )
{
  if ( od.crystal.has_value<CrystalAxis>() ) {
    auto& v = od.crystal.get<CrystalAxis>();
    os << "@crys:" << fmt(v[0]) << ',' << fmt(v[1]) << ',' << fmt(v[2]);
  } else if ( od.crystal.has_value<HKLPoint>() ) {
    auto& v = od.crystal.get<HKLPoint>();
    os << "@crys_hkl:" << fmt(v[0]) << ',' << fmt(v[1]) << ',' << fmt(v[2]);
  } else {
    os << "@crys:<MISSING>";
  }
  os << "@lab:" << fmt(od.lab[0])<< ',' << fmt(od.lab[1])<< ',' << fmt(od.lab[2]);
  return os;
}

std::ostream& NC::operator<< (std::ostream& os, const DensityState& ds)
{
  if ( ds.type == DensityState::Type::SCALEFACTOR ) {
    os << fmt(ds.value)<<"x";
  } else if ( ds.type == DensityState::Type::DENSITY ) {
    os << fmt(ds.value)<<"gcm3";
  } else {
    nc_assert( ds.type == DensityState::Type::NUMBERDENSITY );
    os << fmt(ds.value)<<"perAa3";
  }
  return os;
}

std::ostream& NC::operator<<(std::ostream& os, const UCNMode& ucn )
{
  switch( ucn.mode ) {
  case UCNMode::Mode::Refine:
    os << "refine";
    break;
  case UCNMode::Mode::Remove:
    os << "remove";
    break;
  case UCNMode::Mode::Only:
    os << "only";
    break;
  };
  if ( ucn.threshold != UCNMode::default_threshold() ) {
    //encode in best unit (same logic as in NCCfgVars.hh):
    const double val = ucn.threshold.dbl();
    os << ':';
    if ( val >= 1e-9 && val < 1000e-9 ) {
      os << fmt(val * 1e9) << "neV";
    } else if ( val >= 1e-3 && val < 1.0 ) {
      os << fmt(val * 1e3) << "meV";
    } else {
      os << fmt(val);
    }
  }
  return os;
}

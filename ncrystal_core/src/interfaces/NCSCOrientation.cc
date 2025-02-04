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

#include "NCrystal/interfaces/NCSCOrientation.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"

namespace NC = NCrystal;

void NC::SCOrientation::setPrimaryDirection( const OrientDir& odir)
{
  try {
    precheckLatticeOrientDir( odir );
  } catch (NC::Error::BadInput&e) {
    NCRYSTAL_THROW2(BadInput,"Problem with primary direction: "<<e.what());
  }
  if ( m_dir2data.has_value() )
    precheckLatticeOrientDef( odir, m_dir2data.value().first, m_dir2data.value().second );
  m_dir1data = odir;
}

void NC::SCOrientation::setSecondaryDirection( const OrientDir& odir, double dirtol )
{
  try {
    precheckLatticeOrientDir( odir );
  } catch (NC::Error::BadInput&e) {
    NCRYSTAL_THROW2(BadInput,"Problem with secondary direction: "<<e.what());
  }
  if ( m_dir1data.has_value() )
    precheckLatticeOrientDef( m_dir1data.value(), odir, dirtol );
  m_dir2data.emplace( odir, dirtol );
}

void NC::SCOrientation::stream( std::ostream &os ) const
{
  auto fmtodir = [&os]( const OrientDir& odir )
  {
    auto& c  = odir.crystal;
    auto& l  = odir.lab;
    std::ostringstream ss;
    if ( c.has_value<CrystalAxis>() ) {
      auto& v = c.get<CrystalAxis>();
      os << "@crys:" << fmt(v[0]) << ',' << fmt(v[1]) << ',' << fmt(v[2]);
    } else if ( c.has_value<HKLPoint>() ) {
      auto& v = c.get<HKLPoint>();
      os << "@crys_hkl:" << fmt(v[0]) << ',' << fmt(v[1]) << ',' << fmt(v[2]);
    } else {
      os << "@crys:UNSET";
    }
    os << "@lab:" << fmt(l[0]) << ',' << fmt(l[1]) << ',' << fmt(l[2]);
  };
  const char * unsetfmt = "@crys:UNSET@lab:UNSET";
  os << "SCOrientation(dir1=";
  if ( m_dir1data.has_value() )
    fmtodir( m_dir1data.value() );
  else
    os << unsetfmt;
  os << ";dir2=";
  if ( m_dir2data.has_value() )
    fmtodir( m_dir2data.value().first );
  else
    os << unsetfmt;
  if ( m_dir2data.has_value() )
    os << ";dirtol=" << fmt(m_dir2data.value().second) << ")";
}

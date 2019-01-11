////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2019 NCrystal developers                                   //
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

#include "NCrystal/NCAbsOOV.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCDefs.hh"
#include "NCMath.hh"

NCrystal::AbsOOV::AbsOOV(const Info*ci)
  : Absorption("NCAbsOOV")
{
  nc_assert_always(ci);

  if ( ! ci->hasXSectAbsorption() )
    NCRYSTAL_THROW(MissingInfo,"Info object does not contain absorption cross-section.");

  m_c = sqrt(const_ekin_2200m_s) * ci->getXSectAbsorption();
  validate();
}

NCrystal::AbsOOV::~AbsOOV()
{
}

double NCrystal::AbsOOV::crossSection(double ekin, const double (&)[3] ) const
{
  return ekin ? m_c / std::sqrt(ekin) : std::numeric_limits<double>::infinity();
}

double NCrystal::AbsOOV::crossSectionNonOriented( double ekin ) const
{
  return ekin ? m_c / std::sqrt(ekin) : std::numeric_limits<double>::infinity();
}

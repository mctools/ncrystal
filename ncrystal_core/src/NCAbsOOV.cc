////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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

#include "NCrystal/internal/NCAbsOOV.hh"
#include "NCrystal/internal/NCMath.hh"

namespace NC = NCrystal;

NC::AbsOOV::AbsOOV( const NC::Info& info )
  : AbsOOV( [&info]()
  {
    if ( ! info.hasXSectAbsorption() )
      NCRYSTAL_THROW(MissingInfo,"Info object does not contain absorption cross-section.");
    return info.getXSectAbsorption();
  }() )
{
}

namespace NCrystal {
  namespace {
    constexpr double sqrt_const_ekin_2200m_s = constexpr_sqrt(const_ekin_2200m_s.dbl());
  }
}

NC::AbsOOV::AbsOOV( SigmaAbsorption sigabs)
  : m_c( sqrt_const_ekin_2200m_s * sigabs.dbl() )
{
}


NC::CrossSect NC::AbsOOV::crossSectionIsotropic(CachePtr&, NeutronEnergy ekin ) const
{
  return CrossSect{ ekin.dbl() ? m_c / std::sqrt(ekin.dbl()) : kInfinity };
}

std::shared_ptr<NC::ProcImpl::Process> NC::AbsOOV::createMerged( const Process& oraw ) const
{
  auto o = dynamic_cast<const AbsOOV*>(&oraw);
  if (!o)
    return nullptr;
  auto result = std::make_shared<AbsOOV>( SigmaAbsorption{1.0} );
  result->m_c = this->m_c + o->m_c;
  return result;
}

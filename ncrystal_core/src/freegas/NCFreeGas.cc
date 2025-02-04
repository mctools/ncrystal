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

#include "NCrystal/internal/freegas/NCFreeGas.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/phys_utils/NCFreeGasUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;

struct NC::FreeGas::Impl {

  Impl( Temperature t,
        AtomMass target_mass_amu,
        SigmaFree sigma )
    : m_xsprovider(t, target_mass_amu, sigma),
      m_temperature(DoValidate,t),
      m_target_mass_amu(DoValidate,target_mass_amu)
  {
  }

  FreeGasXSProvider m_xsprovider;
  Temperature m_temperature;
  AtomMass m_target_mass_amu;
};

NC::FreeGas::FreeGas( Temperature t,
                      AtomMass target_mass_amu,
                      SigmaFree sigma )
  : m_impl(t,target_mass_amu,sigma)
{
}

NC::FreeGas::FreeGas( Temperature t,
                      AtomMass target_mass_amu,
                      SigmaBound sb )
  : FreeGas( t, target_mass_amu, sb.free(target_mass_amu) )
{
}

NC::FreeGas::FreeGas( Temperature t, const AtomData& ad )
  : FreeGas( t, ad.averageMassAMU(), ad.freeScatteringXS() )
{
}

NC::FreeGas::~FreeGas() = default;

NC::CrossSect NC::FreeGas::crossSectionIsotropic(CachePtr&, NeutronEnergy ekin ) const
{
  return CrossSect{ m_impl->m_xsprovider.crossSection(ekin) };
}

NC::ScatterOutcomeIsotropic NC::FreeGas::sampleScatterIsotropic(CachePtr&, RNG& rng, NeutronEnergy ekin ) const
{
  double delta_ekin, mu;
  std::tie(delta_ekin,mu) = FreeGasSampler(ekin,m_impl->m_temperature,m_impl->m_target_mass_amu).sampleDeltaEMu(rng);
  return { NeutronEnergy{ncmax(0.0,ekin.get()+delta_ekin)}, CosineScatAngle{mu} };
}

NC::Optional<std::string> NC::FreeGas::specificJSONDescription() const
{
  auto sigmafree = m_impl->m_xsprovider.sigmaFree();
  std::ostringstream ss;
  {
    std::ostringstream tmp;
    tmp << "sigma_free="<<sigmafree<<";T="<<m_impl->m_temperature<<";M="<<m_impl->m_target_mass_amu;
    streamJSONDictEntry( ss, "summarystr", tmp.str(), JSONDictPos::FIRST );
  }
  streamJSONDictEntry( ss, "sigma_free", sigmafree.dbl() );
  streamJSONDictEntry( ss, "temperature", m_impl->m_temperature.dbl() );
  streamJSONDictEntry( ss, "atom_mass", m_impl->m_target_mass_amu.dbl(), JSONDictPos::LAST );
  return ss.str();
}

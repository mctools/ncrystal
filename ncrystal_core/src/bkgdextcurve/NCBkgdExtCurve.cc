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

#include "NCrystal/internal/bkgdextcurve/NCBkgdExtCurve.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"

namespace NC = NCrystal;

NC::BkgdExtCurve::BkgdExtCurve( shared_obj<const Info> ci )
  : m_ci(std::move(ci))
{
  if (!m_ci->providesNonBraggXSects())
    NCRYSTAL_THROW(MissingInfo,"BkgdExtCurve: Passed Info object lacks NonBraggXSects needed for cross sections.");
}

NC::BkgdExtCurve::~BkgdExtCurve() = default;

NC::CrossSect NC::BkgdExtCurve::crossSectionIsotropic(NC::CachePtr&, NC::NeutronEnergy ekin ) const
{
  return m_ci->xsectScatNonBragg(ekin);
}

NC::ScatterOutcomeIsotropic NC::BkgdExtCurve::sampleScatterIsotropic( CachePtr&,
                                                                      RNG& rng,
                                                                      NeutronEnergy ekin ) const
{
  //Elastic, isotropic.
  return { ekin, randIsotropicScatterMu( rng ) };
}

NC::ScatterOutcome NC::BkgdExtCurve::sampleScatter( CachePtr&, RNG& rng, NeutronEnergy ekin, const NeutronDirection& ) const
{
  //Elastic, isotropic.
  return { ekin, randIsotropicNeutronDirection(rng) };
}

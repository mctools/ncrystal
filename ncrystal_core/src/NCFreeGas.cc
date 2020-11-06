////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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

#include "NCrystal/internal/NCFreeGas.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include "NCrystal/internal/NCFreeGasUtils.hh"

namespace NC = NCrystal;

struct NC::FreeGas::Impl {

  Impl( double temp_kelvin,
        double target_mass_amu,
        double sigma,
        SigmaType sigma_type )
    : m_xsprovider(temp_kelvin, target_mass_amu,
                  SigmaFree{sigma_type==SigmaType::FREE
                            ? sigma
                            : std::pow(target_mass_amu/(const_neutron_atomic_mass+target_mass_amu),2)*sigma}),
      m_temperature(temp_kelvin),
      m_target_mass_amu(target_mass_amu)
  {
  }

  FreeGasXSProvider m_xsprovider;
  double m_temperature, m_target_mass_amu;

};

NC::FreeGas::FreeGas( double temp_kelvin,
                      double target_mass_amu,
                      double sigma,
                      SigmaType sigma_type )
  : ScatterIsotropic("FreeGas"), m_impl(temp_kelvin,target_mass_amu,sigma,sigma_type)
{
  validate();
}

NC::FreeGas::~FreeGas() = default;

double NC::FreeGas::crossSection(double ekin, const double (&)[3] ) const
{
  return m_impl->m_xsprovider.crossSection(ekin);
}

double NC::FreeGas::crossSectionNonOriented(double ekin) const
{
  return m_impl->m_xsprovider.crossSection(ekin);
}

void NC::FreeGas::generateScatteringNonOriented( double ekin, double& angle, double& delta_ekin ) const
{
  double mu;
  std::tie(delta_ekin,mu) = FreeGasSampler(ekin,m_impl->m_temperature,m_impl->m_target_mass_amu).sampleDeltaEMu(*getRNG());
  angle = std::acos(mu);
}

void NC::FreeGas::generateScattering( double ekin, const double (&indir)[3],
                                      double (&outdir)[3], double& delta_ekin ) const
{
  RandomBase * rng = getRNG();
  double mu;
  std::tie(delta_ekin,mu) = FreeGasSampler(ekin,m_impl->m_temperature,m_impl->m_target_mass_amu).sampleDeltaEMu(*rng);
  randDirectionGivenScatterMu( rng, mu, indir, outdir );
}

NC::FreeGas::FreeGas( double temp_kelvin, const AtomData& ad )
  : FreeGas( temp_kelvin,
             ad.averageMassAMU(),
             ad.scatteringXS().val,
             SigmaType::BOUND )
{
}


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

#include "NCrystal/internal/NCElIncScatter.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/internal/NCElIncXS.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include "NCrystal/internal/NCDebyeMSD.hh"
#include "NCrystal/internal/NCSpan.hh"
namespace NC = NCrystal;

NC::ElIncScatter::~ElIncScatter() = default;

NC::ElIncScatter::ElIncScatter( const Info* ci )
  : ScatterIsotropic("ElIncScatter")
{
  nc_assert_always(ci);
  if ( !ci->hasAtomInfo() )
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks AtomInfo information"
                   " (elastic-incoherent model only works with crystalline materials).");

  auto atominfos = span<const AtomInfo>(&*ci->atomInfoBegin(),&*ci->atomInfoEnd());

  if ( !ci->hasAtomMSD() ) {
    if ( !ci->hasTemperature() )
      NCRYSTAL_THROW(MissingInfo,"Passed Info object contains neither atomic mean-square-displacements"
                     " (MSD), nor material temperature which is needed for determination of MSDs.");
    if ( !ci->hasAnyDebyeTemperature() )
      NCRYSTAL_THROW(MissingInfo,"Passed Info object contains neither atomic mean-square-displacements"
                     " (MSD), nor Debye temperature info which is needed for determination of MSDs.");
  }

  VectD msd, bixs, scale;
  msd.reserve(atominfos.size());
  bixs.reserve(atominfos.size());
  scale.reserve(atominfos.size());

  unsigned ntot(0);
  for ( const auto& ai : atominfos )
    ntot += ai.number_per_unit_cell;

  for ( const auto& ai : atominfos ) {
    scale.push_back(double(ai.number_per_unit_cell)/ntot);
    bixs.push_back(ai.atom.data().incoherentXS().val);
    if (ai.mean_square_displacement) {
      msd.push_back(ai.mean_square_displacement);
    } else {
      //Fall-back to calculating MSDs from the isotropic Debye model. Eventually
      //we would like to avoid this here, and make sure this is done on the Info
      //object itself.
      double debyeTemp = ai.debye_temp ? ai.debye_temp : ci->getGlobalDebyeTemperature();
      double temperature = ci->getTemperature();
      double atomMass = ai.atom.data().averageMassAMU();
      nc_assert(debyeTemp>0.0&&temperature>0.0&&atomMass>0.0);
      msd.push_back( debyeIsotropicMSD( debyeTemp, temperature, atomMass ) );
    }
  }

  m_elincxs = std::make_unique<ElIncXS>( msd, bixs, scale );
}

NC::ElIncScatter::ElIncScatter( const VectD& elements_meanSqDisp,
                                const VectD& elements_boundincohxs,
                                const VectD& elements_scale )
  : ScatterIsotropic("ElIncScatter")
{
  m_elincxs = std::make_unique<ElIncXS>( elements_meanSqDisp,
                                         elements_boundincohxs,
                                         elements_scale );
}

double NC::ElIncScatter::crossSectionNonOriented(double ekin) const
{
  return m_elincxs->evaluate(ekin);
}

void NC::ElIncScatter::generateScatteringNonOriented( double ekin, double& angle, double& delta_ekin ) const
{
  delta_ekin = 0.0;
  double mu = m_elincxs->sampleMu( getRNG(), ekin );
  nc_assert( mu >= -1.0 && mu <= 1.0 );
  angle = std::acos(mu);
}

void NC::ElIncScatter::generateScattering( double ekin, const double (&indir)[3],
                                           double (&outdir)[3], double& delta_ekin ) const
{
  delta_ekin = 0.0;
  RandomBase* rng = getRNG();
  double mu = m_elincxs->sampleMu( rng, ekin );
  nc_assert( mu >= -1.0 && mu <= 1.0 );
  randDirectionGivenScatterMu( rng, mu, indir, outdir );
}

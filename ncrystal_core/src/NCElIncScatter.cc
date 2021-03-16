
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

#include "NCrystal/internal/NCElIncScatter.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/internal/NCElIncXS.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include "NCrystal/internal/NCDebyeMSD.hh"
#include "NCrystal/internal/NCSpan.hh"
namespace NC = NCrystal;

NC::ElIncScatter::~ElIncScatter() = default;

NC::ElIncScatter::ElIncScatter( const Info& ci )
{
  if ( !ci.hasAtomInfo() )
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks AtomInfo information"
                   " (elastic-incoherent model only works with crystalline materials).");

  auto atominfos = Span<const AtomInfo>(&*ci.atomInfoBegin(),&*ci.atomInfoEnd());

  if ( !ci.hasAtomMSD() ) {
    if ( !ci.hasTemperature() )
      NCRYSTAL_THROW(MissingInfo,"Passed Info object contains neither atomic mean-square-displacements"
                     " (MSD), nor material temperature which is needed for determination of MSDs.");
    if ( !ci.hasDebyeTemperature() )
      NCRYSTAL_THROW(MissingInfo,"Passed Info object contains neither atomic mean-square-displacements"
                     " (MSD), nor Debye temperature info which is needed for determination of MSDs.");
  }

  VectD msd, bixs, scale;
  msd.reserve(atominfos.size());
  bixs.reserve(atominfos.size());
  scale.reserve(atominfos.size());

  unsigned ntot(0);
  for ( const auto& ai : atominfos )
    ntot += ai.numberPerUnitCell();

  for ( const auto& ai : atominfos ) {
    scale.push_back(double(ai.numberPerUnitCell())/ntot);
    bixs.push_back(ai.atomData().incoherentXS().get());
    if ( ai.msd().has_value() ) {
      msd.push_back( ai.msd().value() );
    } else {
      //Fall-back to calculating MSDs from the isotropic Debye model. Eventually
      //we would like to avoid this here, and make sure this is done on the Info
      //object itself.
      nc_assert_always( ai.debyeTemp().has_value() );
      auto debyeTemp = ai.debyeTemp().value();
      auto temperature = ci.getTemperature();
      auto atomMass = ai.atomData().averageMassAMU();
      nc_assert(debyeTemp.get()>0.0&&temperature.get()>0.0&&atomMass.get()>0.0);
      msd.push_back( debyeIsotropicMSD( debyeTemp, temperature, atomMass ) );
    }
  }

  m_elincxs = std::make_unique<ElIncXS>( msd, bixs, scale );
}

NC::ElIncScatter::ElIncScatter( const VectD& elements_meanSqDisp,
                                const VectD& elements_boundincohxs,
                                const VectD& elements_scale )
{
  m_elincxs = std::make_unique<ElIncXS>( elements_meanSqDisp,
                                         elements_boundincohxs,
                                         elements_scale );
}

NC::CrossSect NC::ElIncScatter::crossSectionIsotropic( CachePtr&, NeutronEnergy ekin ) const
{
  return CrossSect{ m_elincxs->evaluate( ekin ) };
}

NC::ScatterOutcomeIsotropic NC::ElIncScatter::sampleScatterIsotropic( CachePtr&, RNG& rng, NeutronEnergy ekin ) const
{
  double mu = m_elincxs->sampleMu( rng, ekin );
  nc_assert( mu >= -1.0 && mu <= 1.0 );
  return { ekin, CosineScatAngle{mu} };
}

NC::ElIncScatter::ElIncScatter( std::unique_ptr<ElIncXS> p )
  : m_elincxs(std::move(p))
{
}

std::shared_ptr<NC::ProcImpl::Process> NC::ElIncScatter::createMerged( const Process& oraw ) const
{
  auto optr = dynamic_cast<const ElIncScatter*>(&oraw);
  if (!optr)
    return nullptr;
  auto& o = *optr;
  nc_assert( m_elincxs != nullptr );
  nc_assert( o.m_elincxs != nullptr );
  return std::make_shared<ElIncScatter>(std::make_unique<ElIncXS>( *m_elincxs, 1.0,
                                                                   *o.m_elincxs, 1.0 ));
}

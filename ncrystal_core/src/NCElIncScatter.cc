
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
#include "NCrystal/internal/NCVDOSEval.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/internal/NCElIncXS.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include "NCrystal/internal/NCDebyeMSD.hh"
#include "NCrystal/internal/NCSpan.hh"
#include "NCrystal/internal/NCDebyeMSD.hh"
#include <iostream>

namespace NC = NCrystal;

NC::ElIncScatter::~ElIncScatter() = default;

NC::ElIncScatter::ElIncScatter( const Info& info, ElIncScatterCfg cfg )
{
  nc_assert_always( !(cfg.scale_factor<=0.0) );
  nc_assert_always( cfg.use_sigma_incoherent || cfg.use_sigma_coherent );
  if ( !info.hasTemperature() )
    NCRYSTAL_THROW(MissingInfo,"Info object passed to ElIncScatter lacks temperature.");

  VectD msd, bixs, scale;
  msd.reserve(info.getComposition().size());
  bixs.reserve(msd.size());
  scale.reserve(msd.size());

  auto getSigma = [&cfg](const AtomData& ad)
  {
    return (  cfg.use_sigma_incoherent ? ad.incoherentXS().get() : 0.0 )
      + ( cfg.use_sigma_coherent ? ad.coherentXS().get() : 0.0 );
  };


  //Prefer initialising via atom infos (in this case we require all to contribute):
  if ( info.hasAtomInfo() ) {

    unsigned ntot(0);
    for ( const auto& ai : info.getAtomInfos() )
      ntot += ai.numberPerUnitCell();

    for ( auto& ai : info.getAtomInfos() ) {
      //scale factor + cross section:
      scale.push_back(double(ai.numberPerUnitCell())*cfg.scale_factor/ntot);
      bixs.push_back( getSigma(ai.atomData()) );
      //msd:
      if ( ai.msd().has_value() ) {
        msd.push_back( ai.msd().value() );
      } else {
        //Fall-back to calculating MSDs from the isotropic Debye model. Eventually
        //we would like to avoid this here, and make sure this is done on the Info
        //object itself.
        if ( !ai.debyeTemp().has_value() || ! info.hasTemperature() )
          NCRYSTAL_THROW(MissingInfo,"Info object passed to ElIncScatter has AtomInfo object without "
                         "mean-square-displacements (MSD), and there is not enough information to"
                         " estimate one (a Debye temperature + material temperature is required).");
        auto debyeTemp = ai.debyeTemp().value();
        auto temperature = info.getTemperature();
        auto atomMass = ai.atomData().averageMassAMU();
        nc_assert(debyeTemp.get()>0.0&&temperature.get()>0.0&&atomMass.get()>0.0);
        msd.push_back( debyeIsotropicMSD( debyeTemp, temperature, atomMass ) );
      }
    }
  } else {
    //Try to initialise via dyninfo sections (ok if some, but not all, have missing info):
    for ( auto& di : info.getDynamicInfoList() ) {
      auto di_vdos = dynamic_cast<const DI_VDOS*>(di.get());
      auto di_vdosdebye = dynamic_cast<const DI_VDOSDebye*>(di.get());
      Optional<double> msd_value;
      if ( di_vdos ) {
        msd_value = VDOSEval( di_vdos->vdosData() ).getMSD();
      } else if ( di_vdosdebye ) {
        if ( !info.hasTemperature() )
          NCRYSTAL_THROW(MissingInfo,"Requested to evaluate atomic mean-squared-displacements for a material without a temperature value.");
        msd_value = debyeIsotropicMSD( di_vdosdebye->debyeTemperature(),
                                       info.getTemperature(),
                                       di_vdosdebye->atomData().averageMassAMU() );
      }
      if ( msd_value.has_value() ) {
        msd.push_back( msd_value.value() );
        scale.push_back( di->fraction() * cfg.scale_factor );
        bixs.push_back( getSigma( di->atomData() ) );
      }
    }
    if ( msd.empty() )
      NCRYSTAL_THROW(MissingInfo,"Info object passed to ElIncScatter lacks information to create Debye-Waller factors.");
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

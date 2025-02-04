
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

#include "NCrystal/internal/elincscatter/NCElIncScatter.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/phys_utils/NCElIncXS.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/phys_utils/NCDebyeMSD.hh"
#include "NCrystal/internal/utils/NCSpan.hh"
#include "NCrystal/internal/phys_utils/NCDebyeMSD.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <sstream>

namespace NC = NCrystal;

NC::ElIncScatter::~ElIncScatter() = default;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    struct ElIncInfo {
      VectD msd, bixs, scale;
    };

    enum class ForCapabilityCheckOnly{ Yes, No };
    Optional<ElIncInfo> extractInfo( const Info& info,
                                     const ElIncScatterCfg& cfg,
                                     ForCapabilityCheckOnly fcco = ForCapabilityCheckOnly::No )
    {
      const bool capOnly = ( fcco == ForCapabilityCheckOnly::Yes );
      if ( capOnly ) {
        if ( cfg.scale_factor == 0.0 )
          return NullOpt;
        if ( ! ( cfg.use_sigma_incoherent || cfg.use_sigma_coherent ) )
          return NullOpt;
        if ( !info.hasTemperature() )
          return NullOpt;
      } else {
        nc_assert_always( !(cfg.scale_factor<=0.0) );
        nc_assert_always( cfg.use_sigma_incoherent || cfg.use_sigma_coherent );
        if ( !info.hasTemperature() )
          NCRYSTAL_THROW(MissingInfo,"Info object passed to ElIncScatter lacks temperature.");
      }

      ElIncInfo res;
      SmallVector<double,6> msd;
      SmallVector<double,6> bixs;
      SmallVector<double,6> scale;

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
            if ( !ai.debyeTemp().has_value() ) {
              if ( capOnly )
                return NullOpt;
              else
                NCRYSTAL_THROW(MissingInfo,"Info object passed to ElIncScatter has AtomInfo object without "
                               "mean-square-displacements (MSD), and there is not enough information to"
                               " estimate one (a Debye temperature + material temperature is required).");
            }
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
        if ( msd.empty() ) {
          if ( capOnly )
            return NullOpt;
          else
            NCRYSTAL_THROW(MissingInfo,"Info object passed to ElIncScatter lacks information to create Debye-Waller factors.");
        }
      }

      if ( capOnly ) {
        //Return non-NullOpt if at least one component has actual contribution:
        for ( auto i : ncrange(bixs.size()) )
          if ( bixs.at(i) * scale.at(i) > 0 )
            return res;
        return NullOpt;
      } else {
        //Copy and return actual results.
        res.msd.reserve( msd.size() );
        res.bixs.reserve( bixs.size() );
        res.scale.reserve( scale.size() );
        res.msd.insert(res.msd.end(), msd.begin(), msd.end() );
        res.bixs.insert(res.bixs.end(), bixs.begin(), bixs.end() );
        res.scale.insert(res.scale.end(), scale.begin(), scale.end() );
        return res;
      }
    }

    class CacheElInc : public CacheBase {
    public:
      void invalidateCache() override { data.ekin.dbl() = -1.0; }
      CacheElInc() { data.ekin.dbl() = -1.0; }
      ElIncXS::EPointAnalysis data;
    };

  }
}

bool NC::ElIncScatter::hasSufficientInfo( const Info& info, const ElIncScatterCfg& cfg )
{
  auto res = extractInfo( info, cfg, ForCapabilityCheckOnly::Yes );
  return res.has_value();
}


NC::ElIncScatter::ElIncScatter( const Info& info, const ElIncScatterCfg& cfg )
{
  auto res = extractInfo( info, cfg );
  if (!res.has_value())
    NCRYSTAL_THROW(MissingInfo,"Info object passed to ElIncScatter lacks information to create Debye-Waller factors.");

  m_elincxs = std::make_unique<ElIncXS>( std::move(res.value().msd),
                                         std::move(res.value().bixs),
                                         std::move(res.value().scale) );
}

NC::ElIncScatter::ElIncScatter( const VectD& elements_meanSqDisp,
                                const VectD& elements_boundincohxs,
                                const VectD& elements_scale )
{
  m_elincxs = std::make_unique<ElIncXS>( elements_meanSqDisp,
                                         elements_boundincohxs,
                                         elements_scale );
}

NC::CrossSect NC::ElIncScatter::crossSectionIsotropic( CachePtr& cp, NeutronEnergy ekin ) const
{
  auto& cache = accessCache<CacheElInc>(cp);
  if ( cache.data.ekin != ekin )
    cache.data = m_elincxs->analyseEnergyPoint( ekin );
  return cache.data.getXS();
}

NC::ScatterOutcomeIsotropic NC::ElIncScatter::sampleScatterIsotropic( CachePtr& cp, RNG& rng, NeutronEnergy ekin ) const
{
  auto& cache = accessCache<CacheElInc>(cp);
  if ( cache.data.ekin != ekin )
    cache.data = m_elincxs->analyseEnergyPoint( ekin );
  return { ekin, cache.data.sampleMu( *m_elincxs, rng ) };
}

NC::ElIncScatter::ElIncScatter( std::unique_ptr<ElIncXS> p )
  : m_elincxs(std::move(p))
{
}

std::shared_ptr<NC::ProcImpl::Process> NC::ElIncScatter::createMerged( const Process& oraw,
                                                                       double scale_self,
                                                                       double scale_other ) const
{
  auto optr = dynamic_cast<const ElIncScatter*>(&oraw);
  if (!optr)
    return nullptr;
  auto& o = *optr;
  nc_assert( m_elincxs != nullptr );
  nc_assert( o.m_elincxs != nullptr );
  return std::make_shared<ElIncScatter>(std::make_unique<ElIncXS>( *m_elincxs, scale_self,
                                                                   *o.m_elincxs, scale_other ));
}

NC::Optional<std::string> NC::ElIncScatter::specificJSONDescription() const
{
  std::ostringstream ss;
  auto xs_lowE = m_elincxs->evaluate( NeutronEnergy{ 0.0 } );
  auto nelem = m_elincxs == nullptr ? 0 : m_elincxs->nElements();
  {
    std::ostringstream tmp;
    tmp << "nelements="<<nelem;
    tmp << ";max_contrib="<<xs_lowE;
    streamJSONDictEntry( ss, "summarystr", tmp.str(), JSONDictPos::FIRST );
  }
  streamJSONDictEntry( ss, "sigma_lowE_limit", xs_lowE.dbl() );
  streamJSONDictEntry( ss, "nelements", nelem, JSONDictPos::LAST );
  return ss.str();
}

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
void NC::ElIncScatter::evalManyXSIsotropic( CachePtr&, const double* ekin, std::size_t N, double* out_xs ) const
{
  m_elincxs->evaluateMany( Span<const double>( ekin, ekin + N ),
                           Span<double>( out_xs, out_xs + N ) );
}

std::pair<NC::CrossSect,NC::ScatterOutcomeIsotropic>
NC::ElIncScatter::evalXSAndSampleScatterIsotropic( CachePtr& cp, RNG& rng, NeutronEnergy ekin ) const
{
  auto& cache = accessCache<CacheElInc>(cp);
  if ( cache.data.ekin != ekin )
    cache.data = m_elincxs->analyseEnergyPoint( ekin );
  return {
    cache.data.getXS(),
    { ekin, cache.data.sampleMu( *m_elincxs, rng ) }
  };
}

#endif

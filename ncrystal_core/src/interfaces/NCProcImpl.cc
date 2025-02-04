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

#include "NCrystal/interfaces/NCProcImpl.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;
namespace NCPI = NCrystal::ProcImpl;

NC::ScatterOutcome NCPI::ScatterIsotropicMat::sampleScatter( CachePtr& cp,
                                                             RNG& rng,
                                                             NeutronEnergy ekin,
                                                             const NeutronDirection& indir ) const
{
  auto outcome_isotropic = sampleScatterIsotropic(cp,rng,ekin);
  auto outdir = randNeutronDirectionGivenScatterMu( rng, outcome_isotropic.mu.get(), indir.as<Vector>() );
  return { outcome_isotropic.ekin, outdir };
}

NC::CrossSect NCPI::ScatterAnisotropicMat::crossSectionIsotropic( CachePtr&, NeutronEnergy ) const
{
  NCRYSTAL_THROW(LogicError,"Process::crossSectionIsotropic can only be called for isotropic materials.");
  return CrossSect{0.0};
}

NC::ScatterOutcomeIsotropic NCPI::ScatterAnisotropicMat::sampleScatterIsotropic( CachePtr&,
                                                                                 RNG&,
                                                                                 NeutronEnergy ) const
{
  NCRYSTAL_THROW(LogicError,"Process::sampleScatterIsotropic can only be called for isotropic materials.");
  return { NeutronEnergy{0.0}, CosineScatAngle{0.0} };
}

NC::ScatterOutcome NCPI::AbsorptionIsotropicMat::sampleScatter( CachePtr&,
                                                                RNG&,
                                                                NeutronEnergy,
                                                                const NeutronDirection& ) const
{
  NCRYSTAL_THROW(LogicError,"Process::sampleScatter can not be called for an absorption process.");
  return { NeutronEnergy{0.0}, NeutronDirection{0,0,1} };
}

NC::ScatterOutcomeIsotropic NCPI::AbsorptionIsotropicMat::sampleScatterIsotropic( CachePtr&,
                                                                                  RNG&,
                                                                                  NeutronEnergy ) const
{
  NCRYSTAL_THROW(LogicError,"Process::sampleScatterIsotropic can not be called for an absorption process.");
  return { NeutronEnergy{0.0}, CosineScatAngle{0.0} };
}

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
bool NCPI::AbsorptionIsotropicMat::isPureElasticScatter() const
{
  NCRYSTAL_THROW(LogicError,"Process::isPureElasticScatter can not be called for an absorption process.");
  return false;
}

std::pair<NC::CrossSect,NC::ScatterOutcome>
NCPI::AbsorptionIsotropicMat::evalXSAndSampleScatter( CachePtr&, RNG&, NeutronEnergy ekin, const NeutronDirection& ndir ) const
{
  NCRYSTAL_THROW(LogicError,"Process::evalXSAndSampleScatter can not be called for an absorption process.");
  return { CrossSect{0.0}, ScatterOutcome{ ekin, ndir } };
}

std::pair<NC::CrossSect,NC::ScatterOutcomeIsotropic>
NCPI::AbsorptionIsotropicMat::evalXSAndSampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ekin ) const
{
  NCRYSTAL_THROW(LogicError,"Process::evalXSAndSampleScatterIsotropic can not be called for an absorption process.");
  return { CrossSect{0.0}, ScatterOutcomeIsotropic::noScat(ekin) };
}
#endif

NC::ScatterOutcomeIsotropic NCPI::NullProcess::sampleScatterIsotropic( CachePtr&,
                                                                       RNG&,
                                                                       NeutronEnergy ekin ) const
{
  return { ekin, CosineScatAngle{1.0} };
}

NC::ScatterOutcome NCPI::NullProcess::sampleScatter( CachePtr&,
                                                     RNG&,
                                                     NeutronEnergy ekin,
                                                     const NeutronDirection& dir ) const
{
  return { ekin, dir };
}

namespace NCRYSTAL_NAMESPACE {
  namespace ProcImpl {

    class CacheProcComp final : public CacheBase {
    public:
      void invalidateCache() override { key_ekin = NeutronEnergy{-1.0}; }

      unsigned nHistory = 0;
      NeutronEnergy key_ekin = NeutronEnergy{-1.0};
      NeutronDirection key_dir = NeutronDirection{0.,0.,0.};//only used for anisotropic materials
      double tot_xs = -1.0;
      struct ComponentCache {
        CachePtr cachePtr;
        EnergyDomain domain;
      };
      SmallVector<ComponentCache,6> componentCache;
      SmallVector<double,6> componentXSectCommul;

      void reset(unsigned nhist,const ProcComposition::ComponentList& comps) {
        nHistory = nhist;
        key_ekin = NeutronEnergy{-1.0};
        key_dir = NeutronDirection{0.,0.,0.};
        tot_xs = -1.0;
        componentCache.clear();
        componentCache.reserve_hint(comps.size());
        for ( auto e : comps )
          componentCache.push_back({{nullptr},e.process->domain()});
        componentXSectCommul.clear();
        componentXSectCommul.resize(comps.size(),0.0);
      }
      CacheProcComp() { reset(nHistory,{}); }
    };

    class ProcComposition::Impl {
    public:
      static CacheProcComp& initAndAccessCache( const ProcComposition* THIS,
                                                CachePtr& cacheptr )
      {
        auto& cache = THIS->accessCache<CacheProcComp>(cacheptr);
        if ( cache.nHistory != THIS->m_nHistory ) {
          //m_components was modified since cache object was created.
          if ( THIS->m_components.empty() )
            NCRYSTAL_THROW(CalcError,"Attempting to use ProcComposition which has no components (if"
                           " intended to be vanishing use a NullProcess component instead).");
          cache.reset(THIS->m_nHistory,THIS->m_components);
        }
        nc_assert(cache.componentCache.size()==THIS->m_components.size());
        nc_assert(cache.componentXSectCommul.size()==THIS->m_components.size());
        return cache;
      }

      static CacheProcComp& updateCacheIsotropic( const ProcComposition* THIS,
                                                  CachePtr& cacheptr,
                                                  NeutronEnergy ekin )
      {
        nc_assert( THIS->m_materialType == MaterialType::Isotropic );
        nc_assert(ekin.dbl()>=0.0);
        nc_assert(THIS->m_domain.contains(ekin) );
        auto& cache = initAndAccessCache(THIS,cacheptr);
        nc_assert(cache.key_dir.as<Vector>().isStrictNullVector());//no mixing between anisotropic and isotropic cache.

        //Compare cached ekin to provided value.
        if ( cache.key_ekin == ekin )
          return cache;

        //Try a bit more FP-sensible cache checking, in case 80-bit registers
        //are somehow messing up stuff (although it is rather unlikely that they
        //will given the non-inlined source of ekin):
        if ( floateq(cache.key_ekin.dbl(),ekin.dbl(),1e-15,0.0) )
          return cache;

        //Ok, cache was not valid!
        cache.key_ekin = NeutronEnergy{-1.0};//put to invalid value while
                                             //updating for exception safety.

        unsigned ncomp = THIS->m_components.size();
        cache.tot_xs = 0.0;
        for ( unsigned i = 0; i < ncomp; ++ i ) {
          auto comp = THIS->m_components[i];
          auto& compCache = cache.componentCache[i];
          CrossSect xs = ( compCache.domain.contains(ekin)
                           ? comp.process->crossSectionIsotropic(compCache.cachePtr,ekin)
                           : CrossSect{0.0} );
          cache.componentXSectCommul[i] = ( cache.tot_xs += ( comp.scale * xs.dbl() ) );
        }

        //All ok:
        cache.key_ekin = ekin;
        return cache;
      }

      static CacheProcComp& updateCacheAnisotropic( const ProcComposition* THIS,
                                                    CachePtr& cacheptr,
                                                    NeutronEnergy ekin,
                                                    const NeutronDirection& dir )
      {
        nc_assert( THIS->m_materialType == MaterialType::Anisotropic );
        nc_assert(ekin.dbl()>=0.0);
        nc_assert(THIS->m_domain.contains(ekin) );
        auto& cache = initAndAccessCache(THIS,cacheptr);

        //Compare cached ekin to provided value.
        if ( cache.key_ekin == ekin && cache.key_dir == dir )
          return cache;

        //Try a bit more FP-sensible cache checking, in case 80-bit registers
        //are somehow messing up stuff (although it is rather unlikely that they
        //will given the non-inlined source of ekin):
        auto cmpfloat = [](double a,double b){ return floateq(a,b,1e-15,0.0); };
        if ( cmpfloat( cache.key_ekin.dbl(),ekin.dbl() )
             && cmpfloat( cache.key_dir[0],dir[0] )
             && cmpfloat( cache.key_dir[1],dir[1] )
             && cmpfloat( cache.key_dir[2],dir[2] ) )
          return cache;

        //Ok, cache was not valid!
        cache.key_ekin = NeutronEnergy{-1.0};//put to invalid value while
                                             //updating for exception safety.

        unsigned ncomp = THIS->m_components.size();
        cache.tot_xs = 0.0;
        for ( unsigned i = 0; i < ncomp; ++ i ) {
          auto comp = THIS->m_components[i];
          auto& compCache = cache.componentCache[i];
          CrossSect xs = ( compCache.domain.contains(ekin)
                           ? comp.process->crossSection(compCache.cachePtr,ekin,dir)
                           : CrossSect{0.0} );
          cache.componentXSectCommul[i] = ( cache.tot_xs += ( comp.scale * xs.get() ) );
        }

        //All ok:
        cache.key_ekin = ekin;
        cache.key_dir = dir;
        return cache;
      }

    };
  }
}

NCPI::ProcComposition::ProcComposition( ComponentList components,
                                        ProcessType processType )
  : m_processType(processType),
    m_materialType(MaterialType::Isotropic)
{
  addComponents(std::move(components));
}

void NCPI::ProcComposition::addComponent( NCPI::ProcPtr process, double scale )
{
  if ( !process ) {
    NCRYSTAL_THROW(BadInput,"Trying to add nullptr component!");
  }
  if ( process->processType() != m_processType ) {
    NCRYSTAL_THROW2(BadInput,"Trying to add "<<process->processType()
                    <<" process to ProcComposition of "<<m_processType<<" processes");
  }
  if ( scale<0.0 || !std::isfinite(scale) )
    NCRYSTAL_THROW2(BadInput,"Trying to add component with invalid scale: "<<scale);

  if ( scale == 0.0 || process->isNull() )
    return;
  auto asproccomp = dynamic_cast<const ProcComposition*>(process.get());
  if (asproccomp) {
    if ( asproccomp == this )
      NCRYSTAL_THROW(BadInput,"It is not allowed to add a ProcComposition object as a component of itself");
    addComponents( {SVAllowCopy,asproccomp->components()}, scale );

    return;
  }
  ++m_nHistory;//record changes to m_components.
  auto expandDomain = [this]( EnergyDomain d )
  {
    if ( d.elow>=d.ehigh )
      return;//d is null range
    if ( m_domain.elow>=m_domain.ehigh ) {
      m_domain = d;
    } else {
      //m_domain and d are both non-null ranges:
      m_domain.elow = NeutronEnergy{ std::min<double>(m_domain.elow.get(),d.elow.get()) };
      m_domain.ehigh = NeutronEnergy{ std::max<double>(m_domain.ehigh.get(),d.ehigh.get()) };
    }
  };

  //Check if the process can be combined with existing components:
  for ( auto& e : m_components ) {
    if ( e.process.get() == process.get() ) {
      //We already have this process, simply increase the corresponding
      //scale (does not affect the domain!):
      e.scale += scale;
      return;
    }
    if ( true ) {
      //Check for processes which can be merged for efficiency (absorbing
      //non-unit scales in the process). This can for example be used to combine
      //multiple PowderBragg instances into one more efficient one which has
      //just one internal hkl list to search - this is potentially quite useful
      //for multiphase materials:
      auto merged_process = e.process->createMerged( *process, e.scale, scale );
      if ( merged_process != nullptr ) {
        e.process = std::move(merged_process);
        e.scale = 1.0;
        if (e.process->materialType()==MaterialType::Anisotropic)
          m_materialType = MaterialType::Anisotropic;
        //Update domain (assume the merging will only have expanded the domain!):
        expandDomain(e.process->domain());
        return;
      }
    }
  }
  //Add new component:
  if (process->materialType()==MaterialType::Anisotropic)
    m_materialType = MaterialType::Anisotropic;
  //Update domain (assume the merging will only have expanded the domain!):
  expandDomain(process->domain());
  //Add to list:
  m_components.push_back(Component{scale,std::move(process)});
}

void NCPI::ProcComposition::addComponents( NCPI::ProcComposition::ComponentList components, double scale )
{
  m_components.reserve_hint( m_components.size() + components.size() );
  for ( auto&& e : components )
    addComponent(std::move(e.process),e.scale*scale);
}

NC::CrossSect NCPI::ProcComposition::crossSection( CachePtr& cacheptr,
                                                   NeutronEnergy ekin,
                                                   const NeutronDirection& dir ) const
{
  if ( ! m_domain.contains(ekin) )
    return CrossSect{ 0.0 };
  auto& cache = ( m_materialType == MaterialType::Anisotropic
                  ? Impl::updateCacheAnisotropic( this, cacheptr, ekin, dir )
                  : Impl::updateCacheIsotropic( this, cacheptr, ekin ) );
  nc_assert( cache.tot_xs >= 0.0 );
  return CrossSect{ cache.tot_xs };
}

NC::CrossSect NCPI::ProcComposition::crossSectionIsotropic( CachePtr& cacheptr,
                                                            NeutronEnergy ekin ) const
{
  if (!m_domain.contains(ekin))
    return CrossSect{ 0.0 };
  nc_assert( m_materialType == MaterialType::Isotropic );
  auto& cache = Impl::updateCacheIsotropic( this, cacheptr, ekin );
  nc_assert( cache.tot_xs >= 0.0 );
  return CrossSect{cache.tot_xs};
}

NC::ScatterOutcome NCPI::ProcComposition::sampleScatter( CachePtr& cacheptr,
                                                         RNG& rng,
                                                         NeutronEnergy ekin,
                                                         const NeutronDirection& dir ) const
{
  if (!m_domain.contains(ekin))
    return { ekin, dir };//no effect when xs=0

  auto& cache = ( m_materialType == MaterialType::Anisotropic
                  ? Impl::updateCacheAnisotropic( this, cacheptr, ekin, dir )
                  : Impl::updateCacheIsotropic( this, cacheptr, ekin ) );
  auto ichoice = pickRandIdxByWeight( rng, cache.componentXSectCommul );
  return m_components[ichoice].process->sampleScatter(cache.componentCache[ichoice].cachePtr,rng,ekin,dir);
}

NC::ScatterOutcomeIsotropic NCPI::ProcComposition::sampleScatterIsotropic( CachePtr& cacheptr,
                                                                           RNG& rng,
                                                                           NeutronEnergy ekin ) const
{
  if (!m_domain.contains(ekin))
    return { ekin, CosineScatAngle{1.0} };//no effect when xs=0
  nc_assert( m_materialType == MaterialType::Isotropic );
  auto& cache = Impl::updateCacheIsotropic( this, cacheptr, ekin );
  auto ichoice = pickRandIdxByWeight( rng, cache.componentXSectCommul );
  return m_components[ichoice].process->sampleScatterIsotropic(cache.componentCache[ichoice].cachePtr,rng,ekin);
}

NCPI::ProcPtr NCPI::ProcComposition::combine( const ComponentList& components,
                                              ProcessType processType )
{
  return consumeAndCombine({SVAllowCopy,components},processType);
}

NCPI::ProcPtr NCPI::getGlobalNullScatter()
{
  static shared_obj<const Process> s_obj = makeSO<NullScatter>();
  return s_obj;
}

NCPI::ProcPtr NCPI::getGlobalNullAbsorption()
{
  static shared_obj<const Process> s_obj = makeSO<NullAbsorption>();
  return s_obj;
}

NCPI::ProcPtr NCPI::ProcComposition::consumeAndCombine( ComponentList&& components_raw,
                                                        ProcessType processType )
{
  auto isNullComponent = [](Component& c) { return c.process==nullptr || c.process->isNull() || c.scale<=0.0; };
  ComponentList components;
  for ( auto& e : components_raw ) {
    if ( isNullComponent(e) )
      continue;
    //Check if more than one identical component:
    auto uid = e.process->getUniqueID();
    bool merged(false);
    for ( auto& e0 : components ) {
      if ( e0.process->getUniqueID() == uid ) {
        merged = true;
        e0.scale += e.scale;
        break;
      }
    }
    if ( merged )
      continue;
    components.emplace_back( std::move( e ) );
  }

  if ( components.empty() ) {
    if ( processType == ProcessType::Scatter )
      return getGlobalNullScatter();
    else
      return getGlobalNullAbsorption();
  }

  //Next, simple cheap check for a single unscaled component:
  if ( components.size() == 1 && components.back().scale == 1.0 )
    return std::move(components.back().process);

  //Non-trivial case, add on ProcComposition:
  auto pc = makeSO<ProcImpl::ProcComposition>( std::move(components), processType );

  //The ProcComposition logic might have merged multiple components into one, so
  //we might again have ended up with just a single component:
  if ( pc->components().size() == 1 && pc->components().back().scale == 1.0 )
    return pc->components().back().process;

  //Alright, return ProcComposition object:
  return pc;
}

void NCPI::Process::initCachePtr(CachePtr& cp) const
{
  //We trigger the cache ptr setup in a brute-force way. The alternative is to
  //add a new required method for all physics models, e.g. prepareCachePtr(..),
  //which would be (slightly) inconvenient.
  cp.reset();
  for ( auto e : { 0.025, 0.0001, 1.0 } ) {
    for ( auto& dir : { NC::NeutronDirection{0.,0.,1.},
                        NC::NeutronDirection{0.,1.,0.},
                        NC::NeutronDirection{1.,0.,1.} } ) {
      crossSection(cp, NeutronEnergy{e}, dir );
      if (cp!=nullptr)
        return;
    }
  }
}

std::string NCPI::Process::jsonDescription() const
{
  std::ostringstream ss;
  streamJSONDictEntry( ss,"name", this->name(), JSONDictPos::FIRST );
  {
    std::ostringstream tmp; tmp << this->materialType();
    streamJSONDictEntry( ss,"materialType", tmp.str() );
  }
  {
    std::ostringstream tmp; tmp << this->processType();
    streamJSONDictEntry( ss,"processType", tmp.str() );
  }
  streamJSONDictEntry( ss,"isOriented", this->isOriented() );
  auto domain = this->domain();
  streamJSONDictEntry( ss,"domain", PairDD{ domain.elow.dbl(), domain.ehigh.dbl() } );
  streamJSONDictEntry( ss,"isNull", this->isNull() );
  auto specific_descr = this->specificJSONDescription();
  if ( specific_descr.has_value() ) {
    ss << ",\"specific\":" << specific_descr.value();
  } else {
    ss << ",\"specific\":{}";
  }
  streamJSONDictEntry( ss,"uid", this->getUniqueID().value, JSONDictPos::LAST );
  return ss.str();
}

NC::Optional<std::string> NCPI::ProcComposition::specificJSONDescription() const
{
  std::ostringstream ss;
  ss << "{\"summarystr\":\""<<m_components.size()<<" components, "<<(isOriented()?"oriented":"isotropic")<<"\"";
  ss << ",\"components\":[";
  bool first(true);
  for ( auto& c : m_components ) {
    if ( first)
      first = false;
    else
      ss << ',';
    ss << '[';
    streamJSON(ss,c.scale);
    ss << ',';
    ss << c.process->jsonDescription();
    ss << ']';
  }
  ss << "]}";
  return ss.str();
}

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE

bool NCPI::Process::isPureElasticScatter() const
{
  return false;
}

std::pair<NC::CrossSect,NC::ScatterOutcome>
NCPI::ScatterIsotropicMat::evalXSAndSampleScatter( CachePtr& cp, RNG& rng,
                                                   NeutronEnergy ekin,
                                                   const NeutronDirection& ndir ) const
{
#if 1
  //This avoids calling randNeutronDirectionGivenScatterMu, in case the process
  //has a more efficient implementation:
  auto xs = this->crossSectionIsotropic( cp, ekin );
  return { xs, this->sampleScatter( cp, rng, ekin, ndir ) };
#else
  //This could on the other hand in principle be more efficient if the process
  //had reimplemented evalXSAndSampleScatterIsotropic but not
  //evalXSAndSampleScatter. Given that high-quality processes might anyway
  //override evalXSAndSampleScatter if it makes sense, we opt to avoid the call
  //to randNeutronDirectionGivenScatterMu, and thus do NOT use the following
  //code:
  auto res = this->evalXSAndSampleScatterIsotropic( cp, rng, ekin );
  return {
    res.first,
    { res.second.ekin,
      randNeutronDirectionGivenScatterMu( rng, res.second.mu, indir ) }
  };
#endif
}

std::pair<NC::CrossSect,NC::ScatterOutcome>
NCPI::ScatterAnisotropicMat::evalXSAndSampleScatter( CachePtr& cp, RNG& rng,
                                                     NeutronEnergy ekin,
                                                     const NeutronDirection& ndir ) const
{
  auto xs = this->crossSection( cp, ekin, ndir );
  return { xs, this->sampleScatter( cp, rng, ekin, ndir ) };
}

std::pair<NC::CrossSect,NC::ScatterOutcomeIsotropic>
NCPI::ScatterIsotropicMat::evalXSAndSampleScatterIsotropic( CachePtr& cp, RNG& rng,
                                                            NeutronEnergy ekin ) const
{
  auto xs = this->crossSectionIsotropic( cp, ekin );
  return { xs, this->sampleScatterIsotropic( cp, rng, ekin ) };
}

std::pair<NC::CrossSect,NC::ScatterOutcomeIsotropic>
NCPI::ScatterAnisotropicMat::evalXSAndSampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ekin ) const
{
  NCRYSTAL_THROW(LogicError,"Process::evalXSAndSampleScatterIsotropic can only be called for isotropic materials.");
  return { CrossSect{0.0}, ScatterOutcomeIsotropic::noScat(ekin) };
}

bool NCPI::ScatterIsotropicMat::isPureElasticScatter() const
{
  return false;
}

bool NCPI::ScatterAnisotropicMat::isPureElasticScatter() const
{
  return false;
}

void NCPI::ScatterIsotropicMat::evalManyXS( CachePtr& cp, const double* ekin,
                                            const double*, const double*, const double*,
                                            std::size_t N,
                                            double* out_xs ) const
{
  //Default implementation simply serialises, and ignores direction:
  for ( std::size_t i = 0; i < N; ++i )
    *out_xs++ = this->crossSectionIsotropic( cp, NeutronEnergy{*ekin++} ).dbl();
}

void NCPI::ScatterAnisotropicMat::evalManyXS( CachePtr& cp, const double* ekin,
                                              const double* ux, const double* uy, const double* uz,
                                              std::size_t N,
                                              double* out_xs ) const
{
  //Default implementation simply serialises:
  for ( std::size_t i = 0; i < N; ++i )
    *out_xs++ = this->crossSection( cp,
                                    NeutronEnergy{*ekin++},
                                    NeutronDirection{ *ux++, *uy++, *uz++ } ).dbl();
}

void NCPI::ScatterIsotropicMat::evalManyXSIsotropic( CachePtr& cp,
                                                     const double* ekin,
                                                     std::size_t N,
                                                     double* out_xs ) const
{
  for ( std::size_t i = 0; i < N; ++i )
    *out_xs++ = this->crossSectionIsotropic( cp, NeutronEnergy{*ekin++} ).dbl();
}

void NCPI::ScatterAnisotropicMat::evalManyXSIsotropic( CachePtr&,
                                                       const double*,
                                                       std::size_t,
                                                       double* ) const
{
  NCRYSTAL_THROW(LogicError,"Process::evalManyXSIsotropic can only be called for isotropic materials.");
}

bool NCPI::ProcComposition::isPureElasticScatter() const
{
  if ( m_processType != ProcessType::Scatter )
    return false;
  for ( auto& e : m_components )
    if ( !e.process->isPureElasticScatter() )
      return false;
  return true;
}

std::pair<NC::CrossSect,NC::ScatterOutcome>
NCPI::ProcComposition::evalXSAndSampleScatter( CachePtr& cachePtr,
                                               RNG& rng,
                                               NeutronEnergy ekin,
                                               const NeutronDirection& dir ) const
{
  if (!m_domain.contains(ekin))
    return { CrossSect{0.0}, { ekin, dir } };
  auto& cache = ( m_materialType == MaterialType::Anisotropic
                  ? Impl::updateCacheAnisotropic( this, cachePtr, ekin, dir )
                  : Impl::updateCacheIsotropic( this, cachePtr, ekin ) );
  auto ichoice = pickRandIdxByWeight( rng, cache.componentXSectCommul );
  return { CrossSect{cache.tot_xs},
           m_components[ichoice].process->sampleScatter(cache.componentCache[ichoice].cachePtr,rng,ekin,dir) };
}

std::pair<NC::CrossSect,NC::ScatterOutcomeIsotropic>
NCPI::ProcComposition::evalXSAndSampleScatterIsotropic( CachePtr& cachePtr,
                                                        RNG& rng,
                                                        NeutronEnergy ekin ) const
{
  if (!m_domain.contains(ekin))
    return { CrossSect{0.0}, { ekin, CosineScatAngle{1.0} } };
  nc_assert( m_materialType == MaterialType::Isotropic );
  auto& cache = Impl::updateCacheIsotropic( this, cachePtr, ekin );
  auto ichoice = pickRandIdxByWeight( rng, cache.componentXSectCommul );
  return { CrossSect{cache.tot_xs},
           m_components[ichoice].process->sampleScatterIsotropic(cache.componentCache[ichoice].cachePtr,rng,ekin) };
}

void NCPI::ProcComposition::evalManyXS( CachePtr& cachePtr, const double* ekin,
                                        const double* ux, const double* uy, const double* uz,
                                        std::size_t N, double* out_xs ) const
{
  auto& cache = Impl::initAndAccessCache(this,cachePtr);
  if ( m_components.size() == 1 ) {
    m_components.front().process->evalManyXS( cache.componentCache[0].cachePtr,
                                              ekin, ux, uy, uz, N, out_xs );
    double scale = m_components.front().scale;
    if ( scale != 1.0 )
      for ( auto i : ncrange(N) )
        out_xs[i] *= scale;
    return;
  }

  for ( auto i : ncrange(N) )
    out_xs[i] = 0.0;

  constexpr std::size_t nbuf = 4096;
  double buf[nbuf];

  while ( N > 0 ) {
    std::size_t nstep = std::min<std::size_t>(N,nbuf);
    for ( auto icomp : ncrange(m_components.size()) ) {
      auto& comp = m_components[icomp];
      comp.process->evalManyXS( cache.componentCache[icomp].cachePtr,
                                ekin, ux, uy, uz, nstep, buf );
      double scale = comp.scale;
      for ( std::size_t j = 0; j<nstep; ++j )
        out_xs[j] += scale * buf[j];
    }
    ekin += nstep;
    ux += nstep;
    uy += nstep;
    uz += nstep;
    out_xs += nstep;
    N -= nstep;
  }
}

void NCPI::ProcComposition::evalManyXSIsotropic( CachePtr& cachePtr,
                                                 const double* ekin,
                                                 std::size_t N,
                                                 double* out_xs ) const
{
  auto& cache = Impl::initAndAccessCache(this,cachePtr);
  if ( m_components.size() == 1 ) {
    m_components.front().process->evalManyXSIsotropic( cache.componentCache[0].cachePtr,
                                                       ekin, N, out_xs );
    double scale = m_components.front().scale;
    if ( scale != 1.0 )
      for ( auto i : ncrange(N) )
        out_xs[i] *= scale;
    return;
  }

  for ( auto i : ncrange(N) )
    out_xs[i] = 0.0;

  constexpr std::size_t nbuf = 4096;
  double buf[nbuf];

  while ( N > 0 ) {
    std::size_t nstep = std::min<std::size_t>(N,nbuf);
    for ( auto icomp : ncrange(m_components.size()) ) {
      auto& comp = m_components[icomp];
      comp.process->evalManyXSIsotropic( cache.componentCache[icomp].cachePtr,
                                         ekin, nstep, buf );
      double scale = comp.scale;
      for ( std::size_t j = 0; j<nstep; ++j )
        out_xs[j] += scale * buf[j];
    }
    ekin += nstep;
    out_xs += nstep;
    N -= nstep;
  }
}

void NCPI::NullProcess::evalManyXS( CachePtr&, const double*,
                                    const double*, const double*, const double*,
                                    std::size_t N, double* out_xs ) const
{
  for ( std::size_t i = 0; i < N; ++ i )
    *out_xs++ = 0.0;
}

void NCPI::NullProcess::evalManyXSIsotropic( CachePtr&,
                                             const double*,
                                             std::size_t N,
                                             double* out_xs ) const
{
  for ( std::size_t i = 0; i < N; ++ i )
    *out_xs++ = 0.0;
}


std::pair<NC::CrossSect,NC::ScatterOutcome>
NCPI::NullAbsorption::evalXSAndSampleScatter( CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const
{
  NCRYSTAL_THROW(LogicError,"Process::evalXSAndSampleScatter can not be called for an absorption process.");
  return { CrossSect{0.0}, { NeutronEnergy{0.0}, NeutronDirection{0,0,1} } };
}

std::pair<NC::CrossSect,NC::ScatterOutcomeIsotropic>
NCPI::NullAbsorption::evalXSAndSampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const
{
  NCRYSTAL_THROW(LogicError,"Process::evalXSAndSampleScatterIsotropic can not be called for an absorption process.");
  return { CrossSect{0.0}, { NeutronEnergy{0.0}, CosineScatAngle{0.0} } };
}

std::pair<NC::CrossSect,NC::ScatterOutcome>
NCPI::NullScatter::evalXSAndSampleScatter( CachePtr&, RNG&, NeutronEnergy ekin, const NeutronDirection& dir ) const
{
  return { CrossSect{0.0}, { ekin, dir } };
}

std::pair<NC::CrossSect,NC::ScatterOutcomeIsotropic>
NCPI::NullScatter::evalXSAndSampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ekin ) const
{
  return { CrossSect{0.0}, { ekin, CosineScatAngle{1.0} } };
}

void NCPI::AbsorptionIsotropicMat::evalManyXS( CachePtr& cp, const double* ekin,
                                               const double*, const double*, const double*,
                                               std::size_t N,
                                               double* out_xs ) const
{
  for ( std::size_t i = 0; i < N; ++i )
    *out_xs++ = this->crossSectionIsotropic( cp, NeutronEnergy{*ekin++} ).dbl();
}

void NCPI::AbsorptionIsotropicMat::evalManyXSIsotropic( CachePtr& cp, const double* ekin, std::size_t N,
                                                        double* out_xs ) const
{
  for ( std::size_t i = 0; i < N; ++i )
    *out_xs++ = this->crossSectionIsotropic( cp, NeutronEnergy{*ekin++} ).dbl();
}


#endif

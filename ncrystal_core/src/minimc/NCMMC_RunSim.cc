////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

#include "NCrystal/internal/minimc/NCMMC_RunSim.hh"
#include "NCrystal/internal/minimc/NCMMC_StdEngine.hh"
#include "NCrystal/internal/minimc/NCMMC_SimMgrMT.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/factories/NCFactImpl.hh"

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

NCMMC::SimOutputMetadata NCMMC::runSim_StdEngine( GeometryPtr geom,
                                                  SourcePtr src,
                                                  TallyPtr tally,
                                                  MatDef matdef,
                                                  const EngineOpts& eopts )
{
  using SimClass = NCMMC::StdEngine;
  if ( eopts.includeAbsorption
       == EngineOpts::IncludeAbsorption::NO )
    matdef.absorption = nullptr;

  auto sim_engine = NC::makeSO<SimClass>( std::move( matdef ),
                                          std::move( eopts ) );
  auto tallymgr = makeSO<TallyMgr>( tally->cloneSetup() );
  NCMMC::SimMgrMT<SimClass> mgr(geom,src,eopts,
                                sim_engine,tallymgr);
  auto missCounts = mgr.launchSimulations( eopts.nthreads,
                                           eopts.seed );
  auto tally_result = tallymgr->getFinalResult();
  tally->merge( std::move( *tally_result ) );

  NCMMC::SimOutputMetadata simoutmd;
  simoutmd.miss = missCounts;
  simoutmd.provided = src->particlesProvided();
  return simoutmd;
}

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace {
      MatDef cfg2MatDef( const MatCfg& cfg )
      {
        auto info = FactImpl::createInfo(cfg);
        auto sct = FactImpl::createScatter(cfg);
        auto absn = FactImpl::createAbsorption(cfg);
        return MatDef( std::move( sct ),
                       std::move( absn ),
                       info->getNumberDensity() );
      }
      void pcsToJSON(std::ostream& os, const ParticleCountSum p)
      {
        streamJSONDictEntry( os, "count",p.count, JSONDictPos::FIRST);
        streamJSONDictEntry( os, "weight",p.weight, JSONDictPos::LAST);
      }

    }
  }
}

NCMMC::MatDef::MatDef( OptionalProcPtr scatter_,
                       OptionalProcPtr absorption_,
                       NumberDensity nd )
  : scatter( std::move( scatter_ ) ),
    absorption( std::move( absorption_ ) ),
    numDens( nd )
{
}

NCMMC::MatDef::MatDef( const MatCfg& cfg )
  : MatDef( cfg2MatDef( cfg ) )
{
  matcfg = cfg;
}

void NCMMC::simOutMetaDataToJSON(std::ostream& os,const SimOutputMetadata& md )
{
  os << "{\"provided\":";
  pcsToJSON(os,md.provided);
  os << ",\"miss\":";
  pcsToJSON(os,md.miss);
  os << '}';
}

void NCMMC::resultsToJSON( std::ostream& os,
                           GeometryPtr geometry,
                           SourcePtr source,
                           TallyPtr tally,
                           MatDef matdef,
                           const EngineOpts& engine_opts,
                           const SimOutputMetadata& simoutmd )
{
  streamJSONDictEntry( os, "datatype", "NCrystalMiniMCResults_v1",
                       JSONDictPos::FIRST);

  os << ",\"input\":{";

  streamJSON( os, "material" );
  os << ':';
  if ( matdef.matcfg.has_value() ) {
    streamJSONDictEntry( os, "cfgstr",
                         matdef.matcfg.value().toStrCfg(),
                         JSONDictPos::FIRST );
    streamJSONDictEntry( os, "decoded",
                         matdef.matcfg.value().toJSONCfg(),
                         JSONDictPos::LAST );
  } else {
    streamJSON( os, json_null_t{} );
  }

  os << ',';
  streamJSON( os, "geometry" );
  os << ':';
  geometry->toJSON(os);

  os << ',';
  streamJSON( os, "source" );
  os << ':';
  source->toJSON(os);

  os << ',';
  streamJSON( os, "engine" );
  os << ':';
  engineOptsToJSON( os, engine_opts );

  os << "},\"output\":{";

  streamJSON( os, "tally" );
  os << ":{";
  {
    bool firstitem(true);
    for ( auto& itemName : tally->tallyItemNames() ) {
      if (firstitem)
        firstitem = false;
      else
        os << ',';
      streamJSON( os, itemName );
      os << ':';
      tally->tallyItemToJSON(os,itemName);
    }
  }
  os << "},";
  streamJSON( os, "metadata" );
  os << ':';
  simOutMetaDataToJSON( os, simoutmd );
  os << "}}";
}

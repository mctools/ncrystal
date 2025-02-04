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

#include "NCrystal/internal/minimc/NCMMC_RunSim.hh"
#include "NCrystal/internal/minimc/NCMMC_StdEngine.hh"
#include "NCrystal/internal/minimc/NCMMC_SimMgrMT.hh"
#include "NCrystal/factories/NCFactImpl.hh"

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

void NCMMC::runSim_StdEngine( ThreadCount nthreads,
                              GeometryPtr geom,
                              SourcePtr src,
                              TallyPtr tally,
                              MatDef matdef,
                              StdEngineOptions engine_options )
{
  using SimClass = NCMMC::StdEngine;
  auto sim_engine = NC::makeSO<SimClass>( std::move( matdef ),
                                          std::move( engine_options ) );
  auto tallymgr = makeSO<TallyMgr>( tally->clone() );
  NCMMC::SimMgrMT<SimClass> mgr(geom,src,sim_engine,tallymgr);
  mgr.launchSimulations( nthreads );
  auto tally_result = tallymgr->getFinalResult();
  tally->merge( std::move( *tally_result ) );
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
    }
  }
}

NCMMC::MatDef::MatDef( ProcPtr scatter_,
                       ProcPtr absorption_,
                       NumberDensity nd )
  : scatter( std::move( scatter_ ) ),
    absorption( std::move( absorption_ ) ),
    numDens( nd )
{
}

NCMMC::MatDef::MatDef( const MatCfg& cfg )
  : MatDef( cfg2MatDef( cfg ) )
{
}

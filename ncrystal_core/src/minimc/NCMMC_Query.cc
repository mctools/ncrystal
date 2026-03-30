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

#include "NCrystal/internal/minimc/NCMMC_Query.hh"
#include "NCrystal/internal/minimc/NCMMC_Utils.hh"
#include "NCrystal/internal/minimc/NCMMC_RunSim.hh"
#include "NCrystal/internal/minimc/NCMMC_Geom.hh"
#include "NCrystal/internal/minimc/NCMMC_Source.hh"
#include "NCrystal/internal/minimc/NCMMC_StdTallies.hh"
#include "NCrystal/internal/minimc/NCMMC_SimEngine.hh"
#include "NCrystal/internal/minimc/NCMMC_EngineOpts.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/factories/NCFactImpl.hh"

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace {
      void query_mmcrun( std::ostream& os,
                         const StrView& cfgstr,
                         const StrView& geomcfg,
                         const StrView& srccfg,
                         const StrView& enginecfg,
                         const Optional<CB::CBMgrInput>& cb )
      {
        MatDef matdef( cfgstr.to_string() );
        auto geom = createGeometry( geomcfg );
        auto src = createSource( srccfg );
        auto eopts = parseEngineOpts( enginecfg );
        auto tally = NC::makeSO<TallyStdHists>( eopts,
                                                src->metaData(),
                                                matdef.matTemp );
        auto resmd = runSim_StdEngine( geom, src, tally, matdef, eopts, cb );
        resultsToJSON( os, geom, src, tally, matdef, eopts, resmd );
      }
    }
  }
}

void NCMMC::Query::JSONQuery( std::ostream& os, const Query& query )
{
  JSONQuery_flexmmcrun( os, query, NullOpt );
}

void NCMMC::Query::JSONQuery_flexmmcrun( std::ostream& os,
                                         const Query& query,
                                         const Optional<CB::CBMgrInput>& cb )
{
  auto invalid = [&query](const char * reason){
    std::ostringstream ss;
    ss << "Invalid MiniMC JSON query: ";
    streamJSON( ss, query );//trick: use streamJSON for easy format.
    if ( reason )
      ss<< " ("<<reason<<')';
    NCRYSTAL_THROW( BadInput, ss.str() );
  };
  constexpr auto sv_mmc = StrView::make("mmc");
  if ( query.size() < 2 || query.front() != sv_mmc ) {
    invalid(nullptr);
    return;
  }
  //shift off the "mmc" and key entries:
  const auto& key = query.at(1).trimmed();
  auto arg = [&query]( std::size_t i ) { return query.at(i+2); };
  auto argstr = [&arg]( std::size_t i ) { return arg(i).to_string(); };
  const std::size_t nargs = static_cast<std::size_t>(query.size()-2);

  constexpr auto sv_list = StrView::make("list");
  constexpr auto sv_run = StrView::make("run");
  constexpr auto sv_scenario = StrView::make("scenario");
  constexpr auto sv_inspectcfg = StrView::make("inspectcfg");
  constexpr auto sv_cfgdoc = StrView::make("cfgdoc");

  if ( cb.has_value() && key != sv_run )
    NCRYSTAL_THROW(BadInput,"MiniMC flexmmcrun only works when combined"
                   "with an ['mmc','run',...] JSON query");

  if ( key == sv_run ) {
    if ( nargs != 4 && nargs != 3 )
      invalid("correct usage: [\"mmc\",\"run\","
              "CFGSTR,GEOMCFG,SRCCFG,ENGINECFG]"
              " or [\"mmc\",\"run\",CFGSTR,SCENARIOCFG,ENGINECFG]");
    if ( nargs == 4 ) {
      query_mmcrun(os,arg(0),arg(1),arg(2),arg(3),cb);
    } else {
      auto matcfg = MatCfg(argstr(0));
      auto d = Utils::decodeScenario( matcfg,argstr(1).c_str() );
      query_mmcrun(os,arg(0),d.geomcfg,d.srccfg,arg(2),cb);
    }
  } else if ( key == sv_scenario ) {
    if ( nargs != 2 )
      invalid("correct usage: [\"mmc\",\"scenario\",CFGSTR,SCENARIOSTR]");
    auto matcfg = MatCfg(argstr(0));
    auto d = Utils::decodeScenario( matcfg, argstr(1).c_str() );
    streamJSONDictEntry( os, "geomcfg", d.geomcfg, JSONDictPos::FIRST );
    streamJSONDictEntry( os, "srccfg", d.srccfg, JSONDictPos::LAST );
  } else if ( key == sv_list ) {
    if ( nargs != 0 )
      invalid("no arguments should come after: [\"mmc\",\"list\"]");
    streamJSON( os,
                std::array<StrView,4>{ sv_run,
                                       sv_scenario,
                                       sv_inspectcfg,
                                       sv_cfgdoc } );
  } else if ( key == sv_inspectcfg ) {
    const char * usage
      = "correct usage: [\"mmc\",\"inspectcfg\",\"src|geom|engine\",STRCFG]";
    if ( nargs != 2 )
      invalid(usage);
    if ( arg(0) == "src" ) {
      createSource( argstr(1).c_str() )->toJSON(os);
    } else if ( arg(0) == "geom" ) {
      createGeometry( argstr(1).c_str() )->toJSON(os);
    } else if ( arg(0) == "engine" ) {
      engineOptsToJSON( os, parseEngineOpts( arg(1) ) );
    } else {
      invalid(usage);
    }
  } else if ( key == sv_cfgdoc ) {
    const char * usage
      = "correct usage: [\"mmc\",\"cfgdoc\",\"src|geom|engine|tally\"]";
    if ( nargs != 1 )
      invalid(usage);
    if ( arg(0) == "src" ) {
      sourceOptsDocsToJSON( os );
    } else if ( arg(0) == "geom" ) {
      geometryOptsDocsToJSON( os );
    } else if ( arg(0) == "engine" ) {
      engineOptsDocsToJSON( os );
    } else if ( arg(0) == "tally" ) {
      tallyHistInfoToJSON(os);
    } else {
      invalid(usage);
    }
  } else {
    invalid(nullptr);
  }
}

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

#include "NCrystal/internal/sab/NCSABRefEval.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/internal/extd_utils/NCInfoUtils.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCSABQuery.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {

    namespace {
      void query_impl_SABRefEval( std::ostream& os,
                                  const MatCfg& cfg,
                                  NeutronEnergy eval,
                                  Optional<std::string> atomDisplayLabel,
                                  std::uint64_t nsample )
      {
        auto info = FactImpl::createInfo(cfg);
        const DynamicInfo* di = InfoUtils::findDynInfo( info, atomDisplayLabel );
        auto di_knl = dynamic_cast<const DI_ScatKnl*>(di);
        if ( !di_knl )
          NCRYSTAL_THROW(BadInput,"Indicated DynInfo object does not provide"
                         " S(alpha,beta) kernels.");
        auto sabdata = extractSABDataFromDynInfo( di_knl, cfg.get_vdoslux() );
        SABRefEval<> refeval( sabdata, eval );
        RNG* rngptr = nullptr;
        std::shared_ptr<RNGStream> rngholder;
        if ( nsample > 0 ) {
          rngholder = getRNG();
          rngptr = rngholder.get();
        }
        refeval.toJSON( os, nsample, rngptr );
      }
    }
  }
}

void NC::SABUtils::JSONQuery( std::ostream& os, const Query& query )
{
  auto invalid = [&query](const char * reason){
    std::ostringstream ss;
    ss << "Invalid SAB JSON query: ";
    streamJSON( ss, query );//trick: using streamJSON for easy format.
    if ( reason )
      ss<< " ("<<reason<<')';
    NCRYSTAL_THROW( BadInput, ss.str() );
  };
  constexpr auto sv_sab = StrView::make("sab");
  if ( query.size() < 2 || query.front() != sv_sab ) {
    invalid(nullptr);
    return;
  }
  //shift off the "sab" and key entries:
  const auto& key = query.at(1).trimmed();
  auto arg = [&query]( std::size_t i ) { return query.at(i+2); };
  auto argstr = [&arg]( std::size_t i ) { return arg(i).to_string(); };
  const std::size_t nargs = static_cast<std::size_t>(query.size()-2);

  constexpr auto sv_list = StrView::make("list");
  constexpr auto sv_refeval = StrView::make("refeval");

  if ( key == sv_refeval ) {
    //query like: ncrystal_query sab refeval 0.025 'bla.ncmat' 1000 ['Al']
    if ( nargs != 4 && nargs != 3 )
      invalid("correct usage: [\"sab\",\"refeval\","
              "ENERGY_EV,MATCFGSTR,NSAMPLE,ATOMDISPLAYLABEL]"
              " (the ATOMDISPLAYLABEL can be left out for monoatomic materials)");
    auto opt_eval = arg(0).toDbl();
    if ( !opt_eval.has_value() || !(opt_eval.value()>0.0) )
      invalid("invalid ENERGY_EV value");
    auto opt_nsample = arg(2).toUInt64();
    if ( !opt_nsample.has_value() )
      invalid("invalid NSAMPLE value");
    Optional<std::string> atomDisplayLabel;
    if ( nargs == 4 )
      atomDisplayLabel = argstr(3);
    query_impl_SABRefEval( os,
                           argstr(1),
                           NeutronEnergy{ DoValidate, opt_eval.value() },
                           atomDisplayLabel,
                           opt_nsample.value() );
  } else if ( key == sv_list ) {
    if ( nargs != 0 )
      invalid("no arguments should come after: [\"mmc\",\"list\"]");
    streamJSON( os, std::array<StrView,1>{ sv_refeval } );
  } else {
    invalid(nullptr);
  }
}

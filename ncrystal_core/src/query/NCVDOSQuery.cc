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

//#include "NCrystal/internal/sab/NCSABRefEval.hh"
//#include "NCrystal/internal/sab/NCSABKBCellSmpl.hh"
//#include "NCrystal/internal/sab/NCSABSurveyor.hh"
//#include "NCrystal/internal/sab/NCSABProcessor.hh"
//#include "NCrystal/internal/sab/NCSABCellInteg.hh"
//#include "NCrystal/internal/sab/NCSABCellSample.hh"
//#include "NCrystal/internal/sab/NCSABRefSampler.hh"
//#include "NCrystal/internal/sab/NCSABExtended.hh"
//#include "NCrystal/internal/sab/NCSABFactory.hh"
//#include "NCrystal/internal/sab/NCSABIntegrator.hh"
//#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
//#include "NCrystal/internal/extd_utils/NCInfoUtils.hh"
//#include "NCrystal/factories/NCFactImpl.hh"
//#include "NCrystal/factories/NCMatCfg.hh"
//#include "NCrystal/core/NCSmallVector.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/vdos/NCVDOSGn.hh"
#include "NCrystal/internal/vdos/NCVDOSToScatKnl.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/vdos/NCVDOSKnlGrid.hh"
#include "NCVDOSQuery.hh"
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace VDOSUtils {
    namespace {

      struct FlexVDOS final : private MoveOnly {
        const VDOSData* vdosData = nullptr;
        int vdoslux = -1;
        bool is_vdosdebye = false;
        //Info:
        std::string cfgstr;
        std::string atomlbl;
        //Lifetime guards:
        std::shared_ptr<const Info> guard_info;
        Optional<VDOSData> guard_vdosData;
      };

      FlexVDOS loadVDOS( const std::string& cfgstr,
                         const std::string& requested_lbl = "" )
      {
        MatCfg matcfg( cfgstr );
        FlexVDOS res;
        res.guard_info = FactImpl::createInfo( matcfg );
        const auto& info = *res.guard_info;
        if ( !info.isSinglePhase() )
          NCRYSTAL_THROW(BadInput,"Single phase material required");
        res.cfgstr = matcfg.toStrCfg();
        res.vdoslux = matcfg.get_vdoslux();
        for ( auto& di : info.getDynamicInfoList() ) {
          const std::string& lbl = info.displayLabel(di->atom().index);
          if ( !requested_lbl.empty() && requested_lbl != lbl )
            continue;
          const VDOSData* vdosData = nullptr;
          auto divdos = dynamic_cast<const DI_VDOS*>(di.get());
          if ( divdos ) {
            res.is_vdosdebye = false;
            vdosData = &divdos->vdosData();
         } else {
            res.is_vdosdebye = true;
            auto divdosdebye = dynamic_cast<const DI_VDOSDebye*>(di.get());
            if ( divdosdebye ) {
              res.guard_vdosData
                = createVDOSDebye( divdosdebye->debyeTemperature(),
                                   divdosdebye->temperature(),
                                   divdosdebye->atomData().scatteringXS(),
                                   divdosdebye->atomData().averageMassAMU() );
              vdosData = &res.guard_vdosData.value();
            }
          }
          if (!vdosData) {
            if ( lbl == requested_lbl )
              NCRYSTAL_THROW(BadInput,"Requested label is present but does"
                             " not have VDOS data");
            continue;
          }
          if ( res.vdosData )
            NCRYSTAL_THROW(BadInput,"Multiple dyninfos with VDOS data present"
                           " and label not provided");
          res.atomlbl = lbl;
          res.vdosData = vdosData;
        }
        if ( res.vdosData == nullptr ) {
          if ( !requested_lbl.empty() ) {
            NCRYSTAL_THROW2(BadInput,"Could not find VDOS data with label \""
                            <<requested_lbl
                            <<"\" (wrong label or not a dyninfo with VDOS)");
          } else {
            NCRYSTAL_THROW(BadInput,"Could not find any VDOS data in input");
          }

        }
        nc_assert_always( res.vdosData != nullptr );
        return res;
      }

      void query_impl_expand( std::ostream& os,
                              const std::string& cfgstr,
                              const std::string& requested_lbl,
                              Optional<VDOSGn::Order> gnmax )
      {
        auto vdos = loadVDOS( cfgstr, requested_lbl );
        const Optional<double> targetEmax;//fixme allow as parameter?
        Optional<unsigned> gnmax_arg;
        if ( gnmax.has_value() )
          gnmax_arg = gnmax.value().value();

        auto gnexpn = expandVDOSToGnFcts( *vdos.vdosData,
                                          vdos.vdoslux,
                                          targetEmax.value_or(0.0),
                                          VDOSGn::TruncAndThinningChoices::Default,
                                          gnmax_arg );

        auto refE0ABGrid = VDOS::setupE0ABGrid( gnexpn, 100000 );

        auto combinedGnFct = VDOS::getCombinedGnFct( gnexpn.Gn );

        unsigned nalpha, nbeta;
        std::tie(nalpha, nbeta) = VDOS::gridDimFromLux( vdos.vdoslux );
        auto abgrid = VDOS::determineAlphaBetaGridFromGn( gnexpn, nalpha, nbeta );

        //fixme make complete:
        os << "{\"input\":{\"cfgstr\":";
        streamJSON(os,vdos.cfgstr);
        os << ",\"lbl\":";
        streamJSON(os,vdos.atomlbl);
        os << ",\"vdoslux\":";
        streamJSON(os,vdos.vdoslux);
        os << ",\"gnmax\":";
        streamJSON(os,gnmax_arg);
        os << ",\"emax\":";
        streamJSON(os,targetEmax);
        os << "},\"output\":{\"alpha2x\":";
        streamJSON(os,gnexpn.alpha2x);
        os << ",\"emax\":";
        streamJSON(os,gnexpn.suggestedEmax);
        os << ",\"alphaMax\":";
        streamJSON(os,gnexpn.upper_alpha);
        os << ",\"betaMax\":";
        streamJSON(os,gnexpn.upper_beta);
        os << ",\"kT\":";
        const double kT = gnexpn.Gn.kT();
        streamJSON(os,kT);
        os << ",\"alphaGrid\":";
        streamJSONHugeDblVect(os,std::move(abgrid.first));
        os << ",\"betaGrid\":";
        streamJSONHugeDblVect(os,std::move(abgrid.second));
        os << ",\"refE0ABGrid\":[";
        streamJSONHugeDblVect(os,std::move(refE0ABGrid.first));
        os <<',';
        streamJSONHugeDblVect(os,std::move(refE0ABGrid.second));
        os << "],\"combinedGnFct\":[";
        streamJSONHugeDblVect(os,std::move(combinedGnFct.first));
        os <<',';
        streamJSONHugeDblVect(os,std::move(combinedGnFct.second));
        os << "],\"Gn\":[";
        const auto& Gn = gnexpn.Gn;
        for ( auto nm1 : ncrange( Gn.maxOrder().value() ) ) {
          const auto n = VDOSGn::Order{nm1 + 1};
          if ( nm1 )
            os << ',';
          os << "{\"n\":";
          streamJSON(os,n.value());
          os << ",\"energy_range\":";
          PairDD range( Gn.eRange(n) );
          streamJSON(os,range);
          // range.first /= kT;
          // range.second /= kT;
          // os << ",\"beta_range\":";
          //streamJSON(os,range);
          os << ",\"values\":";
          streamJSONHugeDblVect(os,VectD(Gn.getRawSpectrum(n)));
          os << '}';
        }
        os << "]}}";
      }
    }
  }
}

void NC::VDOSUtils::JSONQuery( std::ostream& os, const Query& query )
{
  auto invalid = [&query](const char * reason){
    std::ostringstream ss;
    ss << "Invalid VDOS JSON query: ";
    streamJSON( ss, query );//trick: using streamJSON for easy format.
    if ( reason )
      ss<< " ("<<reason<<')';
    NCRYSTAL_THROW( BadInput, ss.str() );
  };
  constexpr auto sv_vdos = StrView::make("vdos");
  if ( query.size() < 2 || query.front() != sv_vdos ) {
    invalid(nullptr);
    return;
  }
  //shift off the "vdos" and key entries:
  const auto& key = query.at(1).trimmed();
  auto arg = [&query]( std::size_t i ) { return query.at(i+2); };
  auto argstr = [&arg]( std::size_t i ) { return arg(i).to_string(); };
  const std::size_t nargs = static_cast<std::size_t>(query.size()-2);

  constexpr auto sv_list = StrView::make("list");
  constexpr auto sv_expand = StrView::make("expand");

  if ( key == sv_expand ) {
    //query like: ncrystal_query vdos gn CFGSTR [GNMAXOVERRIDE]
    if ( nargs<1||nargs>3 )
      invalid("correct usage: [\"vdos\",\"expand\","
              "MATCFGSTR,ATOMDISPLAYLABEL,GNMAXOVERRIDE]"
              " (ATOMDISPLAYLABEL can be left out or as an empty string for"
              " monoatomic materials and GNMAXOVERRIDE can be left out)");
    Optional<VDOSGn::Order> gnmax;
    if ( nargs == 3 ) {
      auto iv = arg(2).toInt32();
      nc_assert_always(iv.has_value()&&iv.value()>=1);
      gnmax.emplace(static_cast<unsigned>(iv.value()));
    }
    std::string lbl;
    if ( nargs >= 2 )
      lbl = arg(1).trimmed().to_string();
    query_impl_expand( os, argstr(0), lbl, gnmax );

  } else if ( key == sv_list ) {
    if ( nargs != 0 )
      invalid("no arguments should come after: [\"mmc\",\"list\"]");
    streamJSON( os, std::array<StrView,2>{ sv_expand } );
  } else {
    invalid(nullptr);
  }
}

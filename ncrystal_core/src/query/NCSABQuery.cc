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
#include "NCrystal/internal/sab/NCSABKBCellSmpl.hh"
#include "NCrystal/internal/sab/NCSABSurveyor.hh"
#include "NCrystal/internal/sab/NCSABCellInteg.hh"//fixme: reconsider filename?
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/internal/extd_utils/NCInfoUtils.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/core/NCSmallVector.hh"
#include "NCSABQuery.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {

    namespace {
      void query_impl_sglcell( std::ostream& os, double E_div_kT,
                               double a1, double a2, double b1, double b2,
                               double s11, double s12, double s21, double s22 )
      {
        PairDD alpha(a1,a2), beta(b1,b2);
        VectD alpha_v = {a1,a2};
        VectD beta_v = {b1,b2};
        SABSurveyor surv( alpha_v, beta_v );
        SABCellSurvey cellsurv( a1, a2, b1, b2, E_div_kT );
        double cellinteg_full;
        SmallVector<std::pair<StrView,double>,32> cellinteg_pb_list;
        {
          CellData cell;//fixme: use this in all relevant interfaces?
          cell.a1 = a1;
          cell.a2 = a2;
          cell.b1 = b1;
          cell.b2 = b2;
          cell.S[0] = s11;
          cell.S[1] = s12;
          cell.S[2] = s21;
          cell.S[3] = s22;
          for ( auto i : ncrange(4) )
            cell.logS[i] = ( cell.S[i] > 0.0 ? std::log(cell.S[i]) : 0.0 );
          StableSum sum_full;
          using SCI = StdLogLinCellIntegrator;
          SCI::integrateFullCell(cell,sum_full);
          cellinteg_full = sum_full.sum();

          StrView sv_allschemes(SCI::allIntegSchemesAsStr());
          for ( auto sv_scheme : sv_allschemes.split(';') ) {
            auto scheme = SCI::str2IntegScheme( sv_scheme );
            StableSum sum_pb;
            SCI::integrateWithinKB( cell, E_div_kT, scheme, sum_pb);
            cellinteg_pb_list.emplace_back(sv_scheme,sum_pb.sum());
          }
        }

        nc_assert_always(surv.getTouchList().size()==1);
        nc_assert_always(surv.getCoverList().size()==1);
        os<<"{\"alpha\":";
        streamJSON(os,alpha);
        os<<",\"beta\":";
        streamJSON(os,beta);
        os<<",\"S\":";
        SmallVector<double,4> s_v = {s11, s12, s21, s22};
        streamJSON(os,s_v);
        os<<",\"surveyor\":{\"E_div_kT_touch\":";
        streamJSON(os,surv.getTouchList().front().first);
        os<<",\"E_div_kT_cover\":";
        streamJSON(os,surv.getCoverList().front().first);

        {
          //Fixme: SABCellEval is not trustworthy -> replace eventually with new
          //code!
          double s[4] = {s11, s12, s21, s22};//fixme check order
          using CellEval = SABCellEval<InterpolationScheme::LOGLIN,
                                       SABInterpolationOrder::ALPHA_FIRST>;
          CellEval ce( alpha, beta, s );
          os<<"},\"celleval_OBSOLETE\":{\"full_integral\":";//fixme: remove this obsolete class
          streamJSON(os,ce.integral());
          os<<",\"phasespace_integral\":";
          streamJSON(os,ce.integralWithinKinematicBounds( E_div_kT ));
          os<<",\"phasespace_E_div_kT\":";
          streamJSON(os,E_div_kT);
        }
        os<<"},\"cellintegral\":{\"full_integral\":";
        streamJSON(os,cellinteg_full);
        os<<",\"phasespace_E_div_kT\":";
        streamJSON(os,E_div_kT);
        os <<",\"integration_regions\":";
        cellsurv.toJSON(os);
        os<<",\"phasespace_integral_format\":[\"scheme\",\"value\"]";
        os<<",\"phasespace_integral\":[";
        bool first=true;
        for ( auto& e : cellinteg_pb_list ) {
          if (!first)
            os<<',';
          first = false;
          os<<'[';
          streamJSON(os,e.first);
          os<<',';
          streamJSON(os,e.second);
          os<<']';
        }
        os <<"]}}";
      }

      void query_impl_surveyor( std::ostream& os,
                                const VectD& alpha,
                                const VectD& beta )
      {
        SABSurveyor surv( alpha, beta );
        os<<"{\"alpha\":";
        streamJSON(os,alpha);
        os<<",\"beta\":";
        streamJSON(os,beta);
        auto streamCellList = [&os](const SABSurveyor::CellList& cl)
        {
          nc_assert_always(cl.size()>=1);
          os << '[';
          bool first = true;
          for ( auto& e : cl ) {
            auto cell_idx = SABSurveyor::unpackCellIdx<unsigned>(e.second);
            if ( first )
              first = false;
            else
              os << ',';
            os<<"[";
            streamJSON(os,e.first);
            os<<',';
            streamJSON(os,cell_idx.first);
            os<<',';
            streamJSON(os,cell_idx.second);
            os<<']';
          }
          os << ']';
        };
        os<<",\"list_format\":[\"E_div_kT\",\"cell_ialpha\",\"cell_ibeta\"]";
        os<<",\"touch_list\":";
        streamCellList(surv.getTouchList());
        os<<",\"cover_list\":";
        streamCellList(surv.getCoverList());
        os<<'}';
      }

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

      class RNGSpy final : public NC::RNGStream {
        shared_obj<RNGStream> m_realrng;
        std::uint64_t m_count = 0;
      public:
        RNGSpy(shared_obj<RNGStream> rs) : m_realrng( std::move(rs) ) {}
        std::uint64_t count() const { return m_count; }
      protected:
        double actualGenerate() override {
          ++m_count;
          return m_realrng->generate(); }
      };

      void query_impl_samplexyparabolicband( std::ostream& os,
                                             double x0, double y0,
                                             double x1, double y1,
                                             std::uint64_t nsample )
      {
        RNGSpy rng( getRNG() );

        ParabolicBandBoxSampler sampler( x0, y0, x1, y1 );
        if (!sampler.canSample())
          nsample = 0;

        nc_assert_always(x0>=0);
        nc_assert_always(y0>=0);
        nc_assert_always(x1>x0);
        nc_assert_always(y1>y0);

        os << "{\"x0\":";
        streamJSON(os,x0);
        os << ",\"y0\":";
        streamJSON(os,y0);
        os << ",\"x1\":";
        streamJSON(os,x1);
        os << ",\"y1\":";
        streamJSON(os,y1);
        os << ",\"samples\":[";
        for ( std::uint64_t i = 0; i < nsample; ++i ) {
          auto res = sampler.sample(rng);
          if ( i )
            os << ',';
          streamJSON(os,res);
        }
        os << "],\"rng_per_sample\":";
        if ( nsample )
          streamJSON(os, double(rng.count())/nsample );
        else
          streamJSON(os, json_null_t{} );
        os << ",\"sampler_details\":";
        sampler.toJSON(os);
        os<<'}';
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
  constexpr auto sv_samplepb = StrView::make("samplepb");
  constexpr auto sv_surveyor = StrView::make("surveyor");
  constexpr auto sv_sglcell = StrView::make("sglcell");
  constexpr auto sv_integschemes = StrView::make("integschemes");

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
  } else if ( key == sv_samplepb ) {
    //sab samplepb x0 y0 x1 y1 nsample
    if ( nargs != 5 )
      invalid("correct usage: [\"sab\",\"samplepb\",X0,Y0,X1,Y1,NSAMPLE");
    auto opt_nsample = arg(4).toUInt64();
    if ( !opt_nsample.has_value() )
      invalid("invalid NSAMPLE value");
    auto s_x0 = arg(0).toDbl();
    auto s_y0 = arg(1).toDbl();
    auto s_x1 = arg(2).toDbl();
    auto s_y1 = arg(3).toDbl();
    if ( !s_x0.has_value() || !s_y0.has_value()
         || !s_x1.has_value() || !s_y1.has_value() )
      invalid("bad coordinate value");
    query_impl_samplexyparabolicband( os, s_x0.value(), s_y0.value(),
                                      s_x1.value(), s_y1.value(),
                                      opt_nsample.value() );
  } else if ( key == sv_surveyor ) {
    const char * usage = ( "correct usage: [\"sab\",\"surveyor\",ALPHAGRID,"
                           "BETAGRID], with grids encoded as '@val1@..@valn" );
    if ( nargs != 2 )
      invalid(usage);
    auto decodeGrid = [](StrView v) -> VectD
    {
      if (!v.startswith('@'))
        return {};
      auto vs = v.substr(1).splitTrimmed('@');
      VectD res;
      res.reserve(vs.size());
      for ( auto& e : vs ) {
        auto val = e.toDbl();
        if (!e.has_value())
          return {};
        res.push_back(val.value());
      }
      return res;
    };
    const VectD alpha = decodeGrid( arg(0) );
    const VectD beta = decodeGrid( arg(1) );
    if ( beta.size() < 2 || alpha.size() < 2)
      invalid(usage);
    query_impl_surveyor( os, alpha, beta );
  } else if ( key == sv_sglcell ) {
    const char * usage = ( "correct usage: [\"sab\",\"sglcell\",EDIVKT,ALPHA1,"
                           "ALPHA2,BETA1,BETA2,SA1B1,SA2B1,SA1B2,SA2B2],"
                           " with negative beta values prefixed with '@'." );
    if ( nargs != 9 )
      invalid(usage);
    if (!arg(3).startswith('@')||!arg(4).startswith('@'))
      invalid(usage);
    SmallVector<Optional<double>,9> v;
    for ( auto i : ncrange(9)) {
      v.push_back( arg(i).startswith('@')
                   ? arg(i).substr(1).toDbl()
                   : arg(i).toDbl() );
    }
    for ( auto e : v ) {
      if (!e.has_value()||ncisnan(e.value())||!std::isfinite(e.value()))
        invalid(usage);
    }
    const double E_div_kT = v.at(0).value();
    const double a1 = v.at(1).value();
    const double a2 = v.at(2).value();
    const double b1 = v.at(3).value();
    const double b2 = v.at(4).value();
    const double s11 = v.at(5).value();
    const double s12 = v.at(6).value();
    const double s21 = v.at(7).value();
    const double s22 = v.at(8).value();
    if ( !(E_div_kT>0.0) || !(a1<a2) || !(a1>=0.0) || !(b2>b1)
         || !(s11>=0.0) || !(s12>=0.0)
         || !(s21>=0.0) || !(s22>=0.0) )
      invalid(usage);
    query_impl_sglcell( os, E_div_kT, a1, a2, b1, b2,
                        s11, s12, s21, s22);
  } else if ( key == sv_integschemes ) {
    if ( nargs != 0 )
      invalid("[\"sab\",\"sglcell\",\"integschemes\"] does"
              " not support any arguments");
    streamJSON(os,StrView(StdLogLinCellIntegrator::
                          allIntegSchemesAsStr()).split(';'));
  } else if ( key == sv_list ) {
    if ( nargs != 0 )
      invalid("no arguments should come after: [\"mmc\",\"list\"]");
    streamJSON( os, std::array<StrView,5>{ sv_integschemes,
                                           sv_refeval,
                                           sv_samplepb,
                                           sv_sglcell,
                                           sv_surveyor } );
  } else {
    invalid(nullptr);
  }
}

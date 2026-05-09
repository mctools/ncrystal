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
  } else if ( key == sv_list ) {
    if ( nargs != 0 )
      invalid("no arguments should come after: [\"mmc\",\"list\"]");
    streamJSON( os, std::array<StrView,2>{ sv_refeval, sv_samplepb } );
  } else {
    invalid(nullptr);
  }
}

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
#include "NCrystal/internal/sab/NCSABProcessor.hh"
#include "NCrystal/internal/sab/NCSABCellInteg.hh"
#include "NCrystal/internal/sab/NCSABCellSample.hh"
#include "NCrystal/internal/sab/NCSABRefSampler.hh"
#include "NCrystal/internal/sab/NCSABExtended.hh"
#include "NCrystal/internal/sab/NCSABFactory.hh"
#include "NCrystal/internal/sab/NCSABIntegrator.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/internal/extd_utils/NCInfoUtils.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/core/NCSmallVector.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCSABQuery.hh"
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {

    namespace {

      shared_obj<const SABData>
      query_impl_extractSabData( const std::string& matcfgstr,
                                 Optional<std::string> atomdsplbl )
      {
        MatCfg cfg(matcfgstr);
        auto info = FactImpl::createInfo(cfg);
        const DynamicInfo* di = nullptr;
        if ( !atomdsplbl.has_value() ) {
          if ( info->getDynamicInfoList().size() > 1 )
            NCRYSTAL_THROW(BadInput, "ATOMDISPLAYLABEL required for"
                           " polyatomic materials");
          di = info->getDynamicInfoList().front().get();
        } else {
          di = InfoUtils::findDynInfo( info, atomdsplbl );
        }
        auto di_knl = dynamic_cast<const DI_ScatKnl*>(di);
        if ( !di_knl )
          NCRYSTAL_THROW(BadInput,"Indicated DynInfo object does not provide"
                         " S(alpha,beta) kernels.");
        return extractSABDataFromDynInfo( di_knl, cfg.get_vdoslux() );
      };

      struct SampleResult {
        VectD a;
        VectD b;
        std::uint64_t ntries = 0;
        double prob1 = -1.0;//probability of edge @ beta=b1
      };
      SampleResult testFullCellOverlaySample( double E_div_kT,
                                              RNG& rng,
                                              const CellData& c,
                                              std::uint64_t nsample,
                                              std::uint64_t ntries_max )
      {
        FullCellSampler fcsampler( &c );
        const double foure = 4*E_div_kT;
        SampleResult res;
        res.prob1 = fcsampler.probabilityEdge1();
        res.a.reserve(nsample);
        res.b.reserve(nsample);
        while( res.a.size() < nsample ) {
          if ( res.ntries > ntries_max ) {
            NCRYSTAL_WARN("Early abort of FullCellOverlaySample due to bad AR.");
            break;
          }
          ++res.ntries;
          auto ab = fcsampler.sampleAlphaBeta( rng );
          if ( ncsquare(ab.first-ab.second) <= ab.first*foure ) {
            res.a.push_back( ab.first );
            res.b.push_back( ab.second );
          }
        }
        return res;
      }

      SampleResult testBoundedCellSample( const SABCellSurvey& cellsurv,
                                          double E_div_kT,
                                          RNG& rng,
                                          const CellData& c,
                                          std::uint64_t nsample,
                                          std::uint64_t ntries_max )
      {
        //Use contributions from each individual cell edge to find probability
        //of scattering from that edge.

        //Calculate prob1 (will be cached in actual production usage alg):
        double prob1;
        const auto scheme
          = SABCfg::IntegrationScheme::MaxPrec;//fixme
        {
          CellData c1 = c;
          c1.S[2]=c1.S[3]=c1.logS[2]=c1.logS[3]=0.0;
          CellData c2 = c;
          c2.S[0]=c2.S[1]=c2.logS[0]=c2.logS[1]=0.0;
          StableSumKahan sum;
          StdLogLinCellIntegrator::integrateWithinKB( c1, E_div_kT,
                                                      scheme, sum );
          const double W1 = sum.sum();
          StdLogLinCellIntegrator::integrateWithinKB( c2, E_div_kT,
                                                      scheme, sum );
          //fixme: can also just use the regular dual-edge integration result in
          //place of W1plusW2.
          const double W1plusW2 = sum.sum();
          nc_assert_always(W1plusW2>0.0);
          prob1 = W1 / W1plusW2;
        }

        auto bcsdata
          = BoundedCellSampler::prepareBCSData( prob1, cellsurv, c, E_div_kT );

        {
          //Small serialisation test:
          MixedDataVector dummy_cache;
          dummy_cache.append(123.0);
          auto packidx = BoundedCellSampler::pack(dummy_cache, bcsdata );
          nc_assert_always(packidx>0&&packidx<dummy_cache.size_bytes());
          dummy_cache.append(456.0);
          std::memset( &bcsdata, 17, sizeof(bcsdata) );
          bcsdata = BoundedCellSampler::unpack(dummy_cache, packidx );
        }

        BoundedCellSampler sampler( c, bcsdata, E_div_kT );
        SampleResult res;
        res.prob1 = prob1;
        res.a.reserve(nsample);
        res.b.reserve(nsample);
        while( res.a.size() < nsample ) {
          auto ab = sampler.sampleAlphaBeta(rng);
          res.a.push_back(ab.alpha);
          res.b.push_back(ab.beta);
          res.ntries += ab.ntries;
          if ( res.ntries > ntries_max ) {
            NCRYSTAL_WARN("Early abort of FixedKBCellSample due to bad AR.");
            break;
          }
        }
        return res;
      }

      SampleResult produceRefCellSamples( double E_div_kT,
                                          RNG& rng,
                                          const CellData& c,
                                          std::uint64_t nsample,
                                          std::uint64_t ntries_max )
      {
        RefCellSampler refsampler( c, E_div_kT );
        const double foure = 4*E_div_kT;
        SampleResult res;
        res.prob1 = -1.0;//not available
        res.a.reserve(nsample);
        res.b.reserve(nsample);
        while( res.a.size() < nsample ) {
          if ( res.ntries > ntries_max )
            NCRYSTAL_THROW(CalcError,"RefCellSampler too inefficient.");
          auto ab = refsampler.sampleAlphaBeta( rng );
          res.ntries += ab.ntries;
          if ( !( ncsquare(ab.alpha-ab.beta) <= ab.alpha*foure ) )
            NCRYSTAL_THROW(CalcError,
                           "RefCellSampler produced inaccessible point.");
          res.a.push_back( ab.alpha );
          res.b.push_back( ab.beta );
        }
        return res;
      }

      void query_impl_proc( std::ostream& os,
                            shared_obj<const SABData> sab,
                            VectD egrid,
                            std::uint64_t nsample,
                            std::uint64_t seed,
                            NeutronEnergy sample_ekin,
                            int knllux )
      {
        if ( knllux < 0 )
          knllux = SABCfg::sablux_default_luxury;

        //Create with diagnostics:
        auto sp = makeSO<SABProcessor>
          ( SABCfg::createConfig(knllux),
            sab,
            std::make_shared<VectD>( std::move(egrid) ),
            SABProcessor::SampleSupport::YES,
            SABProcessor::StoreExtraDiagnostics::YES );
        auto spe = SABExtended::createWithFGExtender( sp );

        os << "{\"sabproc\":";
        spe->processor().toJSON(os);
        os << ",\"sample\":";
        if ( !nsample ) {
          streamJSON(os,json_null_t{});
          os << '}';
          return;
        }
        auto rng = createBuiltinRNG( seed );
        VectD a, b;//, ekin, mu;
        a.reserve(nsample);
        b.reserve(nsample);
        ++nsample;
        while (nsample-- > 1) {
          auto evt =spe->sampleScatterAlphaBeta(rng,sample_ekin);
          a.push_back( evt.first );
          b.push_back( evt.second );
        }
        os << "{\"alpha\":";
        streamJSONHugeDblVect(os,std::move(a));
        os << ",\"beta\":";
        streamJSONHugeDblVect(os,std::move(b));
        os << "}}";
      }

      void query_impl_refsample( std::ostream& os,
                                 shared_obj<const SABData> sab,
                                 std::uint64_t nsample,
                                 std::uint64_t seed,
                                 NeutronEnergy sample_ekin )
      {
        auto rng = createBuiltinRNG( seed );
        const double E_div_kT = sample_ekin.dbl() / sab->temperature().kT();
        auto ab = SABRef::refSampleAlphaBeta( rng, sab, E_div_kT, nsample,
                                              SABRef::RefSampleExtension{} );
        os << "{\"refsample\":{\"E\":";
        streamJSON(os,sample_ekin);
        os << ",\"E_div_kT\":";
        streamJSON(os,E_div_kT);
        os << ",\"kT\":";
        streamJSON(os,sab->temperature().kT());
        os << ",\"seed\":";
        streamJSON(os,seed);
        os << ",\"nsample\":";
        streamJSON(os,nsample);
        os << ",\"alpha\":";
        streamJSONHugeDblVect(os,std::move(ab.first));
        os << ",\"beta\":";
        streamJSONHugeDblVect(os,std::move(ab.second));
        os << "}}";
      }

      void query_impl_legacysample( std::ostream& os,
                                    shared_obj<const SABData> sab,
                                    std::uint64_t nsample,
                                    std::uint64_t seed,
                                    NeutronEnergy sample_ekin,
                                    bool legacy_oversample )
      {
        auto rng = createBuiltinRNG( seed );
        auto scathelper = SAB::
          createScatterHelperWithCache( sab, nullptr,
                                        ( legacy_oversample
                                          ? SAB::LegacySABAlgOpts::OVERSAMPLE10
                                          : SAB::LegacySABAlgOpts::DEFAULT ) );
        auto& sampler = scathelper->sampler;
        const double E_div_kT = sample_ekin.dbl() / sab->temperature().kT();
        VectD alpha, beta;
        alpha.reserve(nsample);
        beta.reserve(nsample);
        for ( std::uint64_t i = 0; i < nsample; ++i ) {
          auto ab = sampler.sampleAlphaBeta( sample_ekin, rng );
          alpha.push_back( ab.first );
          beta.push_back( ab.second );
        }
        os << "{\"legacysample\":{\"E\":";
        streamJSON(os,sample_ekin);
        os << ",\"E_div_kT\":";
        streamJSON(os,E_div_kT);
        os << ",\"kT\":";
        streamJSON(os,sab->temperature().kT());
        os << ",\"seed\":";
        streamJSON(os,seed);
        os << ",\"nsample\":";
        streamJSON(os,nsample);
        os << ",\"alpha\":";
        streamJSONHugeDblVect(os,std::move(alpha));
        os << ",\"beta\":";
        streamJSONHugeDblVect(os,std::move(beta));
        os << "}}";
      }


      void query_impl_sglcell( std::ostream& os, double E_div_kT,
                               double a1, double a2, double b1, double b2,
                               double s11, double s12, double s21, double s22,
                               std::uint64_t nsample, std::uint64_t seed )
      {
        const std::uint64_t ntries_max = nsample*100;
        PairDD alpha(a1,a2), beta(b1,b2);
        VectD alpha_v = {a1,a2};
        VectD beta_v = {b1,b2};
        SABSurveyor surv( alpha_v, beta_v );
        SABCellSurvey cellsurv( a1, a2, b1, b2, E_div_kT );
        double cellinteg_full;
        SmallVector<std::pair<StrView,double>,32> cellinteg_pb_list;
        double cellinteg_pb_chosen_for_fcsample = -1.0;
        CellData cell;//fixme: use this in all relevant interfaces?
        {
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
        }
        {
          StableSumKahan sum_full;
          using SCI = StdLogLinCellIntegrator;
          SCI::integrateFullCell(cell,sum_full);
          cellinteg_full = sum_full.sum();

          StrView sv_allschemes(SABCfg::allIntegSchemesAsStr());
          for ( auto sv_scheme : sv_allschemes.split(';') ) {
            auto scheme = SABCfg::str2IntegScheme( sv_scheme );
            StableSumKahan sum_pb;
            SCI::integrateWithinKB( cell, E_div_kT, scheme, sum_pb);
            cellinteg_pb_list.emplace_back(sv_scheme,sum_pb.sum());
            if ( scheme == SABCfg::IntegrationScheme::Flex5 )
              cellinteg_pb_chosen_for_fcsample = sum_pb.sum();
          }
        }

        const double predicted_fc_ar
          = ( cellinteg_full > 0.0
              ? cellinteg_pb_chosen_for_fcsample/cellinteg_full
              : -1.0 );

        const bool useBoundedCellSample = ( cellinteg_pb_chosen_for_fcsample
                                            < 0.1*cellinteg_full );//fixme: thr?

        nc_assert_always(surv.data().size()==1);
        const double E_div_kT_touch = surv.data().front().e_touch;
        SampleResult samples_fc, samples_bc, samples_ref;
        Optional<double> prob1_fc, prob1_bc;
        auto rng = createBuiltinRNG( seed );
        if ( nsample>0 && E_div_kT>E_div_kT_touch ) {
          samples_ref = produceRefCellSamples( E_div_kT, rng, cell,
                                               nsample, ntries_max );
          nc_assert_always( samples_ref.a.size() == nsample );
          nc_assert_always( samples_ref.b.size() == nsample );
          if ( predicted_fc_ar > 0.001 || !useBoundedCellSample )
            samples_fc = testFullCellOverlaySample( E_div_kT, rng, cell,
                                                    nsample, ntries_max );
          samples_bc = testBoundedCellSample( cellsurv, E_div_kT, rng,
                                              cell, nsample, ntries_max );
          prob1_fc = samples_fc.prob1;
          prob1_bc = samples_bc.prob1;
        }
        os<<"{\"alpha\":";
        streamJSON(os,alpha);
        os<<",\"beta\":";
        streamJSON(os,beta);
        os<<",\"S\":";
        SmallVector<double,4> s_v = {s11, s12, s21, s22};
        streamJSON(os,s_v);
        os<<",\"surveyor\":{\"E_div_kT_touch\":";
        streamJSON(os,surv.data().front().e_touch);
        os<<",\"E_div_kT_cover\":";
        streamJSON(os,surv.data().front().e_cover);

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
        os <<"]},\"sampling\":{";
        os <<"\"chosen_sample_method\":";
        streamJSON(os,(useBoundedCellSample?"bc":"fc"));
        os <<",\"fc_predicted_AR\":";
        streamJSON(os,predicted_fc_ar);
        os <<",\"fc_sampled_alpha\":";
        streamJSONHugeDblVect(os,std::move(samples_fc.a));
        os <<",\"fc_sampled_beta\":";
        streamJSONHugeDblVect(os,std::move(samples_fc.b));
        os <<",\"bc_sampled_alpha\":";
        streamJSONHugeDblVect(os,std::move(samples_bc.a));
        os <<",\"bc_sampled_beta\":";
        streamJSONHugeDblVect(os,std::move(samples_bc.b));
        os <<",\"bc_sampled_ntries\":";
        streamJSON(os,samples_bc.ntries);
        os <<",\"fc_sampled_ntries\":";
        streamJSON(os,samples_fc.ntries);
        os <<",\"fc_prob_edge_1\":";
        streamJSON(os,prob1_fc);
        os <<",\"bc_prob_edge_1\":";
        streamJSON(os,prob1_bc);
        os <<",\"ref_sampled_alpha\":";
        streamJSONHugeDblVect(os,std::move(samples_ref.a));
        os <<",\"ref_sampled_beta\":";
        streamJSONHugeDblVect(os,std::move(samples_ref.b));
        os <<",\"ref_sampled_ntries\":";
        streamJSON(os,samples_ref.ntries);
        os <<"}}";
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
        //Each entry in a CellList is (E/kT,cell index):
        using CellList = std::vector<std::pair<double,SABSurveyor::cellidx_t>>;
        auto sortCellList = [](CellList& cl)
        {
          using E = CellList::value_type;
          std::sort( cl.begin(), cl.end(),
                     []( const E& a, const E& b )
                     {
                       if ( a.first != b.first )
                         return a.first < b.first;
                       return a.second.val < b.second.val;
                     });
        };

        auto streamCellList = [&os](const CellList& cl)
        {
          nc_assert_always(cl.size()>=1);
          os << '[';
          bool first = true;
          for ( auto& e : cl ) {
            if ( first )
              first = false;
            else
              os << ',';
            os<<"[";
            streamJSON(os,e.first);
            os<<',';
            streamJSON(os,e.second.unpackAlphaIdx());
            os<<',';
            streamJSON(os,e.second.unpackBetaIdx());
            os<<']';
          }
          os << ']';
        };
        os<<",\"list_format\":[\"E_div_kT\",\"cell_ialpha\",\"cell_ibeta\"]";
        CellList touchList, coverList;
        {
          touchList.reserve( surv.data().size() );
          for ( auto& cc : surv.data() )
            touchList.emplace_back( cc.e_touch, cc.cellidx );
          sortCellList( touchList );

          coverList.reserve( surv.data().size() );
          for ( auto& cc : surv.data() )
            coverList.emplace_back( cc.e_cover, cc.cellidx );
          sortCellList( coverList );
        }
        os<<",\"touch_list\":";
        streamCellList(touchList);
        os<<",\"cover_list\":";
        streamCellList(coverList);
        os<<'}';
      }

      void query_impl_SABRefEval( std::ostream& os,
                                  const std::string& matcfg,
                                  NeutronEnergy eval,
                                  Optional<std::string> atomDisplayLabel,
                                  std::uint64_t nsample )
      {
        //Fixme: is this obsolete? Should we remove ths option + the
        //NCSABRefEval header again?
        auto sabdata = query_impl_extractSabData( matcfg, atomDisplayLabel );
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
  constexpr auto sv_proc = StrView::make("proc");
  constexpr auto sv_refsample = StrView::make("refsample");

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
                           "ALPHA2,BETA1,BETA2,SA1B1,SA2B1,SA1B2,SA2B2,"
                           "NSAMPLE,SEED], with negative beta values prefixed"
                           " with '@' and NSAMPLE+SEED being optional." );
    std::uint64_t nsample = 0;
    std::uint64_t seed = 123456;
    if ( nargs>=10 ) {
      auto optns = arg(9).toUInt64();
      if ( !optns.has_value() )
        invalid(usage);
      nsample = optns.value();
    }
    if ( nargs>=11 ) {
      auto optsd = arg(10).toUInt64();
      if ( !optsd.has_value() )
        invalid(usage);
      seed = optsd.value();
    }
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
                        s11, s12, s21, s22, nsample, seed );
  } else if ( key == sv_proc || key == sv_refsample ) {
    //Almost same usage+parsing of proc/refsample:
    const bool is_proc = key==sv_proc;
    const char * usage
      = ( is_proc
          ? ( "correct usage: [\"sab\",\"proc\",MATCFGSTR,ATOMDISPLAYLABEL,"
              "NSAMPLE,SEED,SAMPLE_EVAL,EGRID0,..,EGRIDN]"
              " (the EGRID pts are optional, and the ATOMDISPLAYLABEL can"
              " be left as an empty string for monoatomic materials)" )
          : ( "correct usage: [\"sab\",\"refsample\",MATCFGSTR,"
              "ATOMDISPLAYLABEL,NSAMPLE,SEED,SAMPLE_EVAL,TYPE] where TYPE is one"
              " of \"legacy\", \"legacy_oversample\", and \"ref\" (the"
              " ATOMDISPLAYLABEL can be left as an empty string for monoatomic"
              " materials)" ) );
    if ( nargs < 5 )
      invalid(usage);

    Optional<std::string> atomdsplbl;
    if ( !arg(1).empty() )
      atomdsplbl = argstr(1);
    auto sabdata = query_impl_extractSabData( argstr(0), atomdsplbl );
    const int knllux = MatCfg(argstr(0)).get_knllux();

    auto optns = arg(2).toUInt64();
    if ( !optns.has_value() )
      invalid(usage);
    std::uint64_t nsample = optns.value();
    auto optsd = arg(3).toUInt64();
    if ( !optsd.has_value() )
      invalid(usage);
    std::uint64_t seed = optsd.value();
    auto sample_ekin_raw = arg(4).toDbl();
    if (!sample_ekin_raw.has_value()||!(sample_ekin_raw.value()>=0.0))
      invalid(usage);
    NeutronEnergy sample_ekin( DoValidate_t{}, sample_ekin_raw.value() );
    VectD egrid;
    bool refsample_type_is_legacy(false);
    bool legacy_oversample(false);
    if ( !is_proc ) {
      if ( nargs != 6 )
        invalid(usage);
      //parse TYPE
      if ( arg(5) == "ref" ) {
        refsample_type_is_legacy = false;
      } else if ( arg(5) == "legacy" ) {
        refsample_type_is_legacy = true;
      } else if ( arg(5) == "legacy_oversample" ) {
        refsample_type_is_legacy = true;
        legacy_oversample = true;
      } else {
        invalid(usage);
      }
    } else if ( nargs > 5 ) {
      if (!is_proc)
        invalid(usage);
      egrid.reserve(static_cast<std::size_t>(nargs-5));
      for ( std::size_t i = 5; i < nargs; ++i ) {
        auto val = arg(i).toDbl();
        if (!val.has_value()||!(val.value()>0.0))
          invalid(usage);
        egrid.push_back(val.value());
      }
    }
    if ( key == sv_proc ) {
      query_impl_proc( os, std::move(sabdata), std::move(egrid),
                       nsample, seed, sample_ekin, knllux );
    } else {
      if ( refsample_type_is_legacy )
        query_impl_legacysample( os, std::move(sabdata), nsample, seed,
                                 sample_ekin, legacy_oversample );
      else
        query_impl_refsample( os, std::move(sabdata), nsample, seed,
                              sample_ekin );
    }
  } else if ( key == sv_list ) {
    if ( nargs != 0 )
      invalid("no arguments should come after: [\"mmc\",\"list\"]");
    streamJSON( os, std::array<StrView,7>{ sv_integschemes,
                                           sv_proc,
                                           sv_refeval,
                                           sv_refsample,
                                           sv_samplepb,
                                           sv_sglcell,
                                           sv_surveyor } );
  } else {
    invalid(nullptr);
  }
}

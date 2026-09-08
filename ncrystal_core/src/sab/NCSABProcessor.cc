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

#include "NCrystal/internal/sab/NCSABProcessor.hh"
#include "NCrystal/internal/sab/NCSABSurveyor.hh"
#include "NCrystal/internal/sab/NCSABCellInteg.hh"
#include "NCrystal/internal/sab/NCSABCellSample.hh"
#include "NCrystal/internal/sab/NCSABIdx.hh"
#include "NCrystal/internal/utils/NCMixedDataVector.hh"
#include "NCrystal/internal/utils/NCTinyVector.hh"
#include "NCrystal/internal/utils/NCFileUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include <fstream>

namespace NC = NCrystal;
namespace NCS = NCrystal::SABUtils;

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {

    namespace {

#ifndef NDEBUG
      static bool dummy_testsabidx = [](){
        namespace SI = SABIdx;
        nc_assert(SI::PackedIndex::fromAlphaIdxBetaIdx(2,6).unpackAlphaIdx()==2);
        nc_assert(SI::PackedIndex::fromAlphaIdxBetaIdx(2,6).unpackBetaIdx()==6);
        nc_assert(SI::PackedIndex::fromAlphaIdxBetaIdx(0,0).unpackAlphaIdx()==0);
        nc_assert(SI::PackedIndex::fromAlphaIdxBetaIdx(0,0).unpackBetaIdx()==0);
        nc_assert(SI::NAlpha{4}.value() == 4);
        nc_assert(SI::NAlphaCells(SI::NAlpha{4}).value() == 3);
        nc_assert(SI::NAlphaCells{4}.value() == 4);
        nc_assert(SI::NAlpha(SI::NAlphaCells{4}).value() == 5);
        return true;
      }();
#endif
      class CellMgr final : private NoCopyMove {

        void initLogSCache()
        {
#ifndef NDEBUG
          (void)dummy_testsabidx;
#endif
          auto& sab = m_sab->sab();
          nc_assert(m_logS==nullptr);
          m_logS = ncmake_unique_array_noinit<double>(sab.size());
          auto itS = sab.begin();
          auto itSE = sab.end();
          double * itLS = m_logS.get();
          for ( ; itS != itSE; ++itS, ++itLS )
            *itLS = ( *itS ? std::log(*itS) : 0.0 );
        }
        void initCellIntegrals()
        {
          auto& sab = m_sab->sab();
          auto& alpha = m_sab->alphaGrid();
          auto& beta = m_sab->betaGrid();
          //Now calculate cell integrals:
          // const auto nalpha = alpha.size();
          const auto nbeta = beta.size();
          SABIdx::NAlpha nalpha(m_sab->alphaGrid());
          SABIdx::NAlphaCells nalphacells(nalpha);
          nc_assert_always(nbeta*nalpha.value()==sab.size());
          nc_assert_always(nbeta>=2&&nalpha.value()>=2);
          nc_assert( nalphacells.value()+1==alpha.size() );
          const auto nbetacells = nbeta-1;
          const auto ncells = nalphacells.value()*nbetacells;
          m_cellIntegral.reserve(ncells);
          VectD alphaint_cache;
          alphaint_cache.resize(nalphacells.value(),0.0);
          auto itS = sab.begin();
          nc_assert( m_logS != nullptr );
          const double * itLogS = m_logS.get();
          //NB: S(alpha[ia],beta[ib]) = sab[ ib*nalpha+ia ]
          {
            //Init alphaint_cache with cell edge integrals at ib=0:
            auto itA = alpha.begin();
            for ( auto& e : alphaint_cache ) {
              const double a_low = *itA++;
              const double s_low = *itS++;
              const double logS_low = *itLogS++;
              e = integrateAlphaInterval_fast(a_low,s_low,*itA,*itS,
                                              logS_low, *itLogS);
            }
          }

          //fixme: check all _always in file
          for ( auto ib : ncrange(nbetacells) ) {
            //skip to beginning of beta slice:
            ++itS;
            ++itLogS;
            nc_assert( (std::size_t)std::distance(sab.begin(),itS) == (ib+1)*nalpha.value() );
            nc_assert( (std::size_t)std::distance(static_cast<const double*>(m_logS.get()),itLogS) == (ib+1)*nalpha.value() );
            nc_assert( itS == std::next(sab.begin(),(ib+1)*nalpha.value()) );
            nc_assert( itLogS == std::next(m_logS.get(),(ib+1)*nalpha.value()) );
            const double b1 = vectAt(beta,ib);//fixme: use iterators instead
            const double b2 = vectAt(beta,ib+1);
            const double half_db = 0.5 * (b2-b1);
            auto itAIntCache = alphaint_cache.begin();
            auto itAIntCacheE = alphaint_cache.end();
            auto itA = alpha.begin();
            nc_assert(alphaint_cache.size()+1==alpha.size());
            for ( ; itAIntCache!=itAIntCacheE; ++itAIntCache ) {
              nc_assert(itA<alpha.end());
#ifndef NDEBUG
              //Verify consistency of indexing:
              auto dbg_ialpha = static_cast<std::size_t>( itA-alpha.begin() );
#endif
              const double a_low = *itA++;
              const double s_low = *itS++;
              const double logS_low = *itLogS++;
              // nc_assert_always( logS_low == *itLogS++);
              double aint = integrateAlphaInterval_fast(a_low,s_low,*itA,*itS,
                                                        logS_low, *itLogS);
              //combine linearly with previous edge integration to get full cell
              //integral:
              nc_assert( m_cellIntegral.size()
                         == SABIdx::SABCellIdx( nalphacells,
                                                dbg_ialpha, ib ).value() );
              m_cellIntegral.push_back( half_db*(aint+*itAIntCache) );
              //store this edge integration for use during next beta line:
              *itAIntCache = aint;
            }
          }
          nc_assert(m_cellIntegral.size()==ncells);
          nc_assert(std::next(itS)==sab.end());
          nc_assert(std::next(itLogS)==m_logS.get()+sab.size());
        }

      public:
        using cellidx_t = SABIdx::PackedIndex;
        const SABData& sabData() const { return *m_sab; }
        shared_obj<const SABData> sabDataPtr() const { return m_sab; }

        CellMgr( shared_obj<const SABData> sab )
          : m_sab(std::move(sab))
        {
          m_nalpham1 = m_sab->alphaGrid().size()-1;
          initLogSCache();
          initCellIntegrals();
        }

        CellData lookupCellInfo( cellidx_t packed_idx ) const
        {
          auto ia = packed_idx.unpackAlphaIdx();
          auto ib = packed_idx.unpackBetaIdx();
          auto& sab = m_sab->sab();
          auto& alpha = m_sab->alphaGrid();
          auto& beta = m_sab->betaGrid();
          CellData res;
          res.a1 = vectAt(alpha,ia);
          res.a2 = vectAt(alpha,ia+1);
          res.b1 = vectAt(beta,ib);
          res.b2 = vectAt(beta,ib+1);
          SABIdx::NAlpha na(m_sab->alphaGrid());//fixme cache?
          SABIdx::SABIdx idx_sab( na, ia, ib );
          //const auto na = alpha.size();
          //auto idx_sab = ib*na + ia;
          nc_assert( idx_sab.value() == ib*na.value() + ia );
          auto itS = std::next(sab.begin(),idx_sab.value());
          const double * itLogS = m_logS.get() + idx_sab.value();
          nc_assert(std::next(itS)<sab.end());
          nc_assert(std::next(itLogS)<m_logS.get()+sab.size());
          res.S[0] = *itS++;
          res.S[1] = *itS;
          res.logS[0] = *itLogS++;
          res.logS[1] = *itLogS;
          //Advance one row, but we already advanced one item, so only add
          //nalpha-1:
          itS += m_nalpham1;
          itLogS += m_nalpham1;
          nc_assert(std::next(itS)<sab.end());
          nc_assert(std::next(itLogS)<m_logS.get()+sab.size());
          res.S[2] = *itS++;
          res.S[3] = *itS;
          res.logS[2] = *itLogS++;
          res.logS[3] = *itLogS;
          return res;
        }

        double getCellIntegral( SABIdx::PackedIndex p_idx ) const
        {
          SABIdx::SABCellIdx sab_idx( SABIdx::NAlphaCells{m_nalpham1}, p_idx );
          return vectAt(m_cellIntegral, sab_idx.value() );
        }

      private:
        std::unique_ptr<double[]> m_logS;
        VectD m_cellIntegral;
        shared_obj<const SABData> m_sab;
        std::size_t m_nalpham1;
      };

      struct SABProcImpl final : private NoCopyMove {
        using SampleSupport = SABProcessor::SampleSupport;
        using StoreExtraDiagnostics = SABProcessor::StoreExtraDiagnostics;
        SABProcImpl( const SABCfg::Cfg&,
                     shared_obj<const SABData>,
                     std::shared_ptr<const VectD> egrid,
                     SampleSupport,
                     StoreExtraDiagnostics );

        double phaseSpaceIntegral( NeutronEnergy ) const;
        using AlphaBetaOutcome = SABProcessor::AlphaBetaOutcome;
        AlphaBetaOutcome sampleAlphaBeta( RNG&,
                                          NeutronEnergy ) const;
        ScatterOutcomeIsotropic
        ab2scatout( RNG&, NeutronEnergy, const AlphaBetaOutcome& ) const;

        void toJSON( std::ostream& ) const;

        double m_kT;
        double m_invkT;

        //Sorted energy grid (E/kT):
        VectD m_eGrid;
        //Equivalent phase-spaced bounded integral of S for each of the energy
        //points:
        VectD m_sIntegral;

        //Extra diagnostics data
        VectD m_diagnostics_fcAR; //Full cell AR at at each m_eGrid pt

        //We know how to sample at all the energy values in the grid, and the
        //method used is governed by a single idx. If it is negative, it represent
        //minus the number of touched cells, and we can simply use full-cell
        //rejection sampling with all the touched cells. If it is non-negative, a
        //more careful method must be used and the value represents the idx into
        //the m_bcEptInfo vector.
        std::vector<std::int32_t> m_sampleIdx;
        bool canSample() const { return !m_sampleIdx.empty(); }
        ////////////////////////////////////////////
        //Data needed for the full-cell sampling:

        //The cumulative full cell integrals, ordered in the order in which they
        //are touched by neutron phasespace as energy increases.
        VectD m_cumulFCInt;
        //The corresponding packed cell indices:
        std::vector<SABIdx::PackedIndex> m_cumulFCInt_cellidx;

        ///////////////////////////////////////////
        //Data needed for bounded-cell sampling:
        //fixme: not the most efficient structure initially:

        //The sampleIdx gives us a way to select the information needed to
        //sample at a particular energy point (fixme: pack all this more
        //efficiently, avoid large amounts of small vectors and just use offsets
        //into a single large vector?):
        struct BCEnergyPoint {
          //fixme: something different than vectors? ~std::vector<>() uses ~1% of
          //the total time in benchmarks.
          VectD cumulCellContrib;
          struct CellInfo {
            SABIdx::PackedIndex sabCellIndex;
            std::uint32_t bcSampleInfoIndex;//== ::max() for full cell sampling
          };
          std::vector<CellInfo> cellInfo;
        };
        CellMgr m_cellmgr;
      private:
        std::vector<BCEnergyPoint> m_bcEptInfo;
        //sample info (indexed by the .bcSampleInfoIndex properties):
        MixedDataVector m_bcSample;

        //Sample a scattering event:
        struct ABRes { double a, b; std::size_t ntries; };
        ABRes sampleAlphaBetaFC( RNG&, double E_div_kT,
                                 std::size_t nCellsTouched ) const;
        ABRes sampleAlphaBetaBC( RNG&, double E_div_kT,
                                 const BCEnergyPoint& ) const;

      };

      struct sIntegralAtE_Result {
        double integralWithinKB;
        double integralTouchedCells;

        //Detailed crossed-cell integrals are also recorded, in case they are
        //needed for BoundedCell sampling later. Note that these are in the same
        //order as on the original SABSurveyor object.
        std::vector<std::pair<double,SABSurveyor::cellidx_t>> crossedIntegrals;

        //fixme: better with : <0 is crossedIntegral, >=0 means full cell
        //integral. And for final cumul-value storage we should sort first so
        //the lowest contrib comes first.

        //Number of cells touched (they will all be located at the beginning of
        //SABSurveyor::data()):
        std::size_t nTouchedCells;

        //Internal cache, speeding up repeated calls if they are done from low
        //to high energy values:
        struct PrevPtInfo {
          double E_div_kT = -1.0;//-1.0 means "n/a"
          StableSumKahan sumFullyCoveredCells;
          std::size_t nLastTouched = 0;
        } prevCache;

      };

      void sIntegralAtE( double E_div_kT,
                         const CellMgr& mgr,
                         const SABSurveyor& surv,
                         StdLogLinCellIntegrator::IntegrationScheme scheme,
                         sIntegralAtE_Result& result,
                         bool do_sample )
      {
        auto& prev = result.prevCache;
        const double prev_E = ( ( prev.E_div_kT >= 0.0
                                  && prev.E_div_kT <= E_div_kT )
                                ? prev.E_div_kT : -1.0 );
        const double prev_contrib( prev_E >= 0.0
                                   ? prev.sumFullyCoveredCells.sum() : 0.0 );


        result.crossedIntegrals.clear();
        nc_assert(E_div_kT>0&&std::isfinite(E_div_kT));
        //Find all touched cells:
        auto survCells = surv.data();
        //std::upper bound returns the first element where the comparision
        //returns true. The list must be sorted so the comparison will return
        //[false,false,...,false,true,true,...,true].

        //we find the first element coming AFTER E_div_kT

        //FIXME: Why do we need the binary search? Why not just proceed over all
        //cells until e < ci.e_touch?

        auto itSearchBegin = survCells.begin();
        if ( prev_E >= 0.0 && prev.nLastTouched > 2 )
          itSearchBegin += prev.nLastTouched-1;//fixme: is this actually helping in a measurable way??

        auto itLastTouchedE
          = std::upper_bound(itSearchBegin,
                             survCells.end(),
                             E_div_kT,
                             [](double e,
                                const SABSurveyor::CellInfo& ci)
                             { return e < ci.e_touch; } );
        //upper_bound returns one past the last touched:
        if( itLastTouchedE == survCells.begin() )
          NCRYSTAL_THROW2(BadInput,"Phase-space of neutron energy E/kT="
                          <<fmt(E_div_kT)
                          <<" does not touch any S(alpha,beta) cells!");
        //fixme: should we require all sab tables grids to intersect the E=0
        //phasespace line?
        StableSumKahan sum;
        StableSumKahan sumFullCells;
        StableSumKahan sumCoveredCells;
        result.nTouchedCells = std::distance(survCells.begin(),itLastTouchedE);
        //cells with relative max contrib below this can never contribute:
        // (fixme: sablux should modify 1e10 in range [1e6,1e17])
        const double threshold = ncclamp( 1.0/(1e10*result.nTouchedCells),
                                          1e-30,1e-6 );

        for ( auto it = survCells.begin(); it!=itLastTouchedE; ++it) {
          nc_assert( E_div_kT >= it->e_touch );
          if ( prev_E >= it->e_cover ) {
            nc_assert( E_div_kT >= it->e_cover );
            //Was already fully covered at previous energy point, so the same
            //will be true at this energy point. We already have the accumulated
            //info for this cell included in the .prevCache, so skip it here.
            nclikely continue;
          }

          const double fullCellIntegral = mgr.getCellIntegral( it->cellidx );
#define SUMADD_NEW//fixme ensure this is correct and then enable
#ifndef SUMADD_NEW
          sumFullCells.add(fullCellIntegral);
#endif
          if ( E_div_kT >= it->e_cover ) {
            //fully covers cells:
            sumCoveredCells.add(fullCellIntegral);
          } else {
            double contrib;
            if ( fullCellIntegral < (sumFullCells.sumUncorrected()+prev_contrib)*threshold ) {
              contrib = 0.0;//don't waste time on irrelevant cell
            } else {
              CellData cellData = mgr.lookupCellInfo( it->cellidx );
              StableSumKahan crossedRes;
              StdLogLinCellIntegrator::integrateWithinKB( cellData, E_div_kT,
                                                          scheme, crossedRes );
              sum.add( crossedRes );
              contrib = crossedRes.sum();
            }
            if ( do_sample ) {
              result.crossedIntegrals.emplace_back( contrib, it->cellidx );
            }
          }
        }
#ifdef SUMADD_NEW
        sumFullCells.add(sumCoveredCells);
        sumFullCells.add(sum);
#endif

        //Handle interactions with prev_E cache:
        if ( prev_E >= 0.0 )
          sumCoveredCells.add(prev.sumFullyCoveredCells);

        //Prepare cache for next:
        prev.sumFullyCoveredCells.set(sumCoveredCells);
        prev.E_div_kT = E_div_kT;
        prev.nLastTouched = result.nTouchedCells;
        sum.add(sumCoveredCells);
        sumFullCells.add(sumCoveredCells);


        result.integralWithinKB = sum.sum();
        result.integralTouchedCells = sumFullCells.sum();
        nc_assert(std::isfinite(result.integralWithinKB));
        nc_assert(std::isfinite(result.integralTouchedCells));
      }


      //Automatically determine EMax by looking at the full cell integrals of
      //cells touched at a given energy:
      NeutronEnergy determineEMax( Temperature temp,
                                   //AtomMass,//fixme: add if needed (for FG model perhaps?)
                                   const VectD& cumulFullCellIntegral,
                                   const VectD& cellTouch_E_div_kT )
      {
        const double kT = temp.kT();
        const VectD& e = cellTouch_E_div_kT;
        const VectD& s = cumulFullCellIntegral;
        const double extra_safety = 0.98;

        //First, discard trailing entries where the integral(S) remains
        //~constant, since that is a sign that increased phase-space does not
        //add new contributions to the integral. Thus, we discard until we have
        //removed eps_emax (in logIntegral space) of the integral(S) difference:

        nc_assert_always(s.size()==e.size());
        nc_assert_always(s.size()>=16);
        const double s1(s.front());
        const double s2(s.back());
        nc_assert_always( s2>s1 && s1>=0.0 && std::isfinite(s2) );
        double tgt_s;
        constexpr double eps_emax = 1e-2;
        if ( s1==0.0 ) {
          tgt_s = s2*(1.0-eps_emax);
          } else {
          tgt_s = s1 * std::pow( s2/s1, 1.0-eps_emax );
        }
        std::size_t i = s.size()-1;
        while ( i && vectAt(s,i)>tgt_s )
          --i;
        nc_assert_always( i > 4 );//fixme nicer error

        auto max_E_div_kT = vectAt(e,i);
        nc_assert( max_E_div_kT > 0.0 );
        return NeutronEnergy{ ncmin( max_E_div_kT*kT*extra_safety,
                                     1000.0,
                                     10000*kT ) };
      }


      //Given f(e) = Sintegral(e)/sqrt(e), and a search range [e1,e2] find in as
      //few f(e) evaluations as possible, the highest possible emin value where
      //f(emin)~=f(e1). It is assumed that [e1,e2] is already narrow enough that
      //f(e) is a monotonic function over the interval.
      double determineEMinDivKT( const SABCfg::Cfg& cfg,
                                 const PairDD& e_search_range,
                                 const std::function<double(double)>& f_of_e )
      {
        //fixme: this function could use some luxury parameters
        const double root_acc = cfg.egrid_emin_accuracy;

        const double e1 = e_search_range.first;
        const double e2 = e_search_range.second;
        const double logf1 = std::log(f_of_e( e1 ));
        const double logf2 = std::log(f_of_e( e2 ));
        const double logS_raise_target = 1e-3*(logf2-logf1);//fixme: lux?
        nc_assert_always(logS_raise_target > 0.0);
        const double logf1pluseps = logf1 + logS_raise_target;
        if ( logf2 <= logf1pluseps )
          return e2;
        const double loge1 = std::log(e1);
        const double loge2 = std::log(e2);
        unsigned long ncalls = 2;
        // NCRYSTAL_MSG("TKTEST "
        //              <<" logS_raise_target="<<fmt(logS_raise_target)
        //              <<" logf1pluseps="<<fmt(logf1pluseps)
        //              <<" logf1="<<fmt(logf1)
        //              <<" logf2="<<fmt(logf2)
        //              <<" loge1="<<fmt(loge1)
        //              <<" loge2="<<fmt(loge2));
        auto froot = [&f_of_e,logf1pluseps,&ncalls,
                      logf1,logf2,loge1,loge2]( double loge )
        {
          ++ncalls;
          nc_assert(std::isfinite(loge));
          double val;
          if ( !(loge>loge1) ) {
            val = logf1;
          } else if ( !(loge<loge2) ) {
            val = logf2;
          } else {
            val = std::log(f_of_e(std::exp(loge)));
          }
          return val - logf1pluseps;
          //          return std::log(f_of_e(std::exp(loge)))-logf1pluseps;
        };
        if ( false ) {//fixme
          std::string fn = "ncrystal_logemin_info.txt";
          NCRYSTAL_WARN("Writing to "<<fn
                        <<" if it does not already exist");
          if (!file_exists(fn)) {
            std::ofstream ofs(fn.c_str(), std::ofstream::out);
            ofs << "#ncrystal_xycurve\n";
            ofs << "#colnames = loge;deltalogf\n";
            for ( auto loge : linspace(loge1,loge2,100000) )
              ofs<<fmt(loge)<<" "<<fmt(froot(loge))<<"\n";
            ofs.close();
          }
        }

        double logemin;
        try {
          logemin = findRoot2(froot,loge1,loge2, root_acc);
        } catch ( NC::Error::CalcError& e ) {
          logemin = -1.0;
          //Write fct to file for debugging:
          std::string fn = "ncrystal_logemin_fail_info.txt";
          NCRYSTAL_WARN("Failed to find emin via root-finding. Will attempt"
                        " to write debugging info to "<<fn
                        <<" if it does not already exist");
          if (!file_exists(fn)) {
            std::ofstream ofs(fn.c_str(), std::ofstream::out);
            ofs << "#ncrystal_xycurve\n";
            ofs << "#colnames = loge;deltalogf\n";
            for ( auto loge : linspace(loge1,loge2,100000) )
              ofs<<fmt(loge)<<" "<<fmt(froot(loge))<<"\n";
            ofs.close();
          }
          throw;
        }

        return ncclamp(std::exp(logemin),e1,e2);
      }

      //Based on full-cell integrals, determine range in which a more refined
      //search for EMin/kT should take place.
      PairDD determineEMinDivKTminRange( const VectD& cumulFullCellIntegral,
                                         const VectD& cellTouch_E_div_kT,
                                         const double EMax_div_kT )
      {
        (void)EMax_div_kT;//fixme??
        const VectD& e = cellTouch_E_div_kT;
        const VectD& s = cumulFullCellIntegral;

        //We check for the point where the S-integral has increased some amount
        //from the initial value.

        nc_assert_always(s.size()==e.size());
        nc_assert_always(s.size()>=16);
        const double s1(s.front());
        const double s2(s.back());
        nc_assert_always( s2>s1 && s1>=0.0 && std::isfinite(s2) );
        const std::size_t ns = s.size();

        std::size_t i1 = 0;
        {
          constexpr double eps_emin1 = 1e-6;
          const double tgt_s1 = (s1?s1*std::pow(s2/s1,eps_emin1):s2*eps_emin1);
          while ( i1<ns && vectAt(s,i1)<tgt_s1 )
            ++i1;
          nc_assert_always( i1+8 < ns );
        }

        std::size_t i2 = i1+1;
        {
          constexpr double eps_emin2 = 0.1;
          const double tgt_s2 = (s1?s1*std::pow(s2/s1,eps_emin2):s2*eps_emin2);
          while ( i2<ns && vectAt(s,i2)<tgt_s2 )
            ++i2;
          nc_assert_always( i2+4 < ns );
        }

        constexpr double safety1 = 0.9;
        PairDD res( vectAt(e,i1)*safety1, vectAt(e,i2) );
        nc_assert_always( std::isfinite(res.first) );
        nc_assert_always( std::isfinite(res.second) );
        nc_assert_always( std::isfinite(res.first) );
        nc_assert_always( res.second >= res.first );
        return res;
      }

      struct SIntByETouch final : private MoveOnly {
        VectD sint, e;
      };

      SIntByETouch determineSIntegralByTouchedCells( const CellMgr& cellmgr,
                                                     const SABSurveyor& surv,
                                                     double kT )
      {
        auto it = surv.data().begin();
        auto itE = surv.data().end();
        nc_assert_always( it != itE );

        SIntByETouch res;
        res.e.reserve(4096);
        res.sint.reserve(4096);
        double eprev(0.0);
        double curr_thr(0.0);
        const double thr_eup = 1e3/kT;//1keV/kT is essentially ~infinity for our
                                      //purposes
        StableSumKahan sum;

        nc_assert_always(it!=itE);
        for ( ;it!=itE;++it ) {
          double ci = cellmgr.getCellIntegral(it->cellidx);
          if ( it->e_touch > thr_eup )
            break;
          if ( !(ci>curr_thr) )
            continue;
          sum.add(ci);
          //fixme: tune?
          curr_thr = sum.sum()*1e-26;//less than this won't ever affect result (~16
          //digits + O(1e9) cells => ~O(1e-25) is enough). (fixme: tighten?)
          if ( it->e_touch > eprev*1.001 ) {
            eprev = it->e_touch;
            res.e.push_back(it->e_touch);
            res.sint.push_back(sum.sum());
          }
        }

        return res;
      }

      //fixme: clarify that this returns in E/kT, but requested_egrid is in eV
      VectD determineEGrid( const SABCfg::Cfg& cfg,
                            const SABData& sab,
                            const CellMgr& cellmgr,
                            const SABSurveyor& surv,
                            std::shared_ptr<const VectD> requested_egrid )
      {
        const auto scheme = cfg.integSchemeDetermineEGrid;
        static_assert( std::is_same<decltype(cfg.egrid_npts),unsigned>
                       ::value, "" );
        auto npts = static_cast<std::size_t>( cfg.egrid_npts );

        //fixme: major cleanup needed
        Optional<NeutronEnergy> suggestedEMax_by_egrid;
        Optional<NeutronEnergy> suggestedEMin_by_egrid;

        if ( requested_egrid && !requested_egrid->empty() ) {
          const auto& re = *requested_egrid;
          //fixme: more sanity checks
          const double invkT = 1.0/sab.temperature().kT();
          if ( re.size() == 1 ) {
            suggestedEMax_by_egrid = NeutronEnergy{ re.front() };
          } else if ( re.size() == 3 ) {
            //fixme: the NCMAT spec says that this should give points
            //distributed "evenly in logarithmic space". Perhaps we can solve
            //this merely by updating the language in the spec, to be a bit more
            //loose?

            // { emin or 0, emax or 0, npts or 0 }
            if ( re.at(0) != 0.0 )
              suggestedEMin_by_egrid = NeutronEnergy{ re.at(0) };
            if ( re.at(1) != 0.0 )
              suggestedEMax_by_egrid = NeutronEnergy{ re.at(1) };
            if ( re.at(2) != 0.0 ) {
              npts = static_cast<std::size_t>(re.at(2));
              if ( (double)npts != re.at(2) || npts < 10  )
                NCRYSTAL_THROW2(BadInput,"Invalid egrid npts (must be"
                                " integral value, at least 10): "<<re.at(2));
            }
          } else {
            //complete egrid given directly by request, so in this case we
            //simply honour it and return it directly (we can only assume that
            //the users know what they are doing):
            if ( !( re.size() > 10 ) )
              NCRYSTAL_THROW2(BadInput,"Invalid requested egrid (must have"
                              " length 1, 3, or >=10)");
            VectD res;
            res.reserve( re.size() );
            for ( auto e : re )
              res.push_back( e * invkT );
            return res;
          }
        }

        {
          nc_assert_always( npts >= 10 );
          if ( suggestedEMax_by_egrid.has_value() )
            suggestedEMax_by_egrid.value().validate();
          if ( suggestedEMin_by_egrid.has_value() )
            suggestedEMin_by_egrid.value().validate();
          if ( suggestedEMax_by_egrid.has_value()
               && suggestedEMin_by_egrid.has_value()
               && !( suggestedEMax_by_egrid.value()
                     > suggestedEMin_by_egrid.value() ) ) {
            NCRYSTAL_THROW2(BadInput,"egrid: emax ("
                            <<suggestedEMax_by_egrid.value()
                            <<") not greater than emin ("
                            <<suggestedEMin_by_egrid.value()<<")");
          }
        }

        auto sIntByEtouch
          = determineSIntegralByTouchedCells( cellmgr, surv,
                                              sab.temperature().kT() );
        const double kT = sab.temperature().kT();

        Optional<NeutronEnergy> suggestedEMin = NullOpt;//fixme
        Optional<NeutronEnergy> suggestedEMax;
        if ( sab.suggestedEmax() > 0.0 )
          suggestedEMax = sab.suggestedEmax();//fixme this should not compile!
        if ( suggestedEMax_by_egrid.has_value() )
          suggestedEMax = suggestedEMax_by_egrid;

        std::pair<NeutronEnergy,NeutronEnergy> E_min_max;
        E_min_max.second = ( suggestedEMax.has_value()
                             ? suggestedEMax.value()
                             : determineEMax( sab.temperature(), sIntByEtouch.sint, sIntByEtouch.e ) );

        if ( suggestedEMin_by_egrid.has_value() )
          suggestedEMin = suggestedEMin_by_egrid;

        if ( suggestedEMin.has_value() ) {
          E_min_max.first = suggestedEMin.value();
        } else {
          auto r = determineEMinDivKTminRange( sIntByEtouch.sint, sIntByEtouch.e,
                                               E_min_max.second.get()/kT );

          sIntegralAtE_Result tmp_result;
          tmp_result.crossedIntegrals.reserve(1024);
          auto f = [&cellmgr, &surv, &tmp_result,scheme]( double E_div_kT ) {
            nc_assert(E_div_kT>0.0);
            const double se = std::sqrt(E_div_kT);
            nc_assert(se>0.0);
            sIntegralAtE( E_div_kT, cellmgr, surv, scheme, tmp_result, false );
            return tmp_result.integralWithinKB / se;
          };
          E_min_max.first = NeutronEnergy{ kT*determineEMinDivKT( cfg, r, f ) };
        }
        if ( E_min_max.first >= E_min_max.second ) {
          //fixme: warning!
          E_min_max.first = NeutronEnergy( 0.1 * E_min_max.second.get() );
        }
        if ( E_min_max.first.get() < 1e-50*E_min_max.second.get() ) {
          //fixme: warning!
          E_min_max.first = NeutronEnergy( 1e-50 * E_min_max.second.get() );
        }

        VectD final_egrid;
        {
          //at lowE int(S)(E) = const*sqrt(E).
          //
          //So to have an acceptance rate of 0.9 between E pts, we need E points to be
          //spaced with a factor of k, given by:
          //
          //  const*sqrt(E) = 0.9*const*sqrt(E*k)
          //  <=> sqrt(k) = 1/0.9 <=> k = 1/0.9**2
          //
          // Meaning that N pts will cover a factor k^N values. So to cover En/E0=1e4,
          // one will need:
          //
          // k^N = 1e4 <=> N = ln(1e4)/lnk = ln(1e4)/ln(1/0.9**2) ~= 43.7 ~= 44
          //
          // If on the other hand we wanted to cover En/E0=1e6 with AR>=0.95, we would
          // need ~134 pts. This would also be OK, but note that the quality of
          // sampled results is the same.
          constexpr std::size_t n_lowe = 44;//fixme: we should consider a
                                            //lin-in-sqrt(E) interpolation in
                                            //the region of these cells, and a
                                            //lin-in-E above. Since we have
                                            //already determined a "lowE" cutoff
                                            //this seems sensible. For
                                            //user-provided egrid, we can for
                                            //simplicity just use lin-in-sqrt(E)
                                            //for the lower 5% of epts (or we
                                            //can detect)
          if ( npts < 80 )
            npts = 80;
          const std::size_t n = npts-n_lowe;
          constexpr double lowE_factor = 1e-4;
          final_egrid.reserve(n+n_lowe);
          for ( auto ee : geomspace(E_min_max.first.get()*lowE_factor,
                                    E_min_max.first.get(), n_lowe) )
            final_egrid.push_back(ee);
          final_egrid.pop_back();
          for ( auto ee : geomspace(E_min_max.first.get(),
                                    E_min_max.second.get(), n) )
            final_egrid.push_back(ee);
        }

        const double invkT = 1.0/kT;
        for ( auto& ee : final_egrid )
          ee *= invkT;

        return final_egrid;
      }

      SABProcImpl::BCEnergyPoint
      prepareBCEPt( double E_div_kT, const SABSurveyor& surv,
                    const CellMgr& mgr, const sIntegralAtE_Result& integAtE,
                    SABCfg::IntegrationScheme scheme,
                    MixedDataVector& commonStorage )
      {
        nc_assert( commonStorage.size_bytes() < static_cast<std::size_t>
                   (std::numeric_limits<std::uint32_t>::max()) );
        const double e = E_div_kT;
        SABProcImpl::BCEnergyPoint res;

        auto itNextCrossed = integAtE.crossedIntegrals.begin();
        using TmpE = std::pair<double,SABProcImpl::BCEnergyPoint::CellInfo>;
        std::vector<TmpE> tmp;
        tmp.reserve(2*integAtE.crossedIntegrals.size());

        for( auto& cell : surv.data() ) {
          if ( e < cell.e_touch )//fixme: it used to be "if ( e > cell.e_touch ) break" ?????
            break;//not touching any more cells
          double fullCellInteg = mgr.getCellIntegral( cell.cellidx );
          double contrib;
          if ( e >= cell.e_cover ) {
            contrib = fullCellInteg;
          } else {
            nc_assert( itNextCrossed != integAtE.crossedIntegrals.end() );
            nc_assert( itNextCrossed->second.val == cell.cellidx.val );
            contrib = itNextCrossed->first;
            ++itNextCrossed;
          }

          if ( !(contrib>0.0) )
            continue;//skip irrelevant cell

          tmp.emplace_back();
          auto& entry = tmp.back();
          entry.first = contrib;
          entry.second.sabCellIndex = cell.cellidx;

          //If full cell sampling acceptance rate is high enough, mark this cell
          //for full-cell sampling .
          constexpr double fc_threshold = 0.2;//fixme: tune
          if ( contrib >= fc_threshold*fullCellInteg ) {
            //Full cell sampling acceptance rate is good enough:
            entry.second.bcSampleInfoIndex
              = std::numeric_limits<std::uint32_t>::max();
          } else {
            //Needs bounded cell sampling:
            using BCS = BoundedCellSampler;

            CellData cellData = mgr.lookupCellInfo( cell.cellidx );
            SABCellSurvey cellSurvey( cellData.a1, cellData.a2,
                                      cellData.b1, cellData.b2,
                                      E_div_kT );
            double prob1;
            {
              CellData c1 = cellData;
              c1.S[2]=c1.S[3]=c1.logS[2]=c1.logS[3]=0.0;
              CellData c2 = cellData;
              c2.S[0]=c2.S[1]=c2.logS[0]=c2.logS[1]=0.0;
              StableSumKahan sum;
              StdLogLinCellIntegrator::integrateWithinKB( c1, E_div_kT,
                                                          scheme, sum );
              const double W1 = sum.sum();
              StdLogLinCellIntegrator::integrateWithinKB( c2, E_div_kT,
                                                          scheme, sum );
              //fixme: can also just use the regular dual-edge integration
              //result in place of W1plusW2.
              const double W1plusW2 = sum.sum();
              nc_assert(W1plusW2>0.0);
              prob1 = W1 / W1plusW2;
              nc_assert( prob1 >= 0.0 );
              nc_assert( prob1 <= 1.0 );
            }

            BCS::BCSData bcsdata = BCS::prepareBCSData( prob1, cellSurvey,
                                                        cellData, E_div_kT );
            std::size_t pack_idx = BCS::pack( commonStorage, bcsdata );
            nc_assert( pack_idx < std::numeric_limits<std::uint32_t>::max() );
            entry.second.bcSampleInfoIndex = static_cast<std::uint32_t>(pack_idx);
          }
        }
        nc_assert( itNextCrossed==integAtE.crossedIntegrals.end() );
        nc_assert( commonStorage.size_bytes()
                   < std::numeric_limits<std::uint32_t>::max() );

        //For best numerical stability, sort so smallest contrib comes first in
        //the cumul vector. We do not need stable sort, since sabCellIndex
        //should be unique.
        std::sort( tmp.begin(), tmp.end(),
                   [](const TmpE& a, const TmpE& b) {
                     if ( a.first != b.first )
                       return a.first < b.first;//by ascending contrib
                     //rare contrib tie (written to support self-comparison):
                     if ( a.second.sabCellIndex.val
                          != b.second.sabCellIndex.val )
                       return ( a.second.sabCellIndex.val
                                < b.second.sabCellIndex.val );
                     return ( a.second.bcSampleInfoIndex
                              < b.second.bcSampleInfoIndex );
                   } );

        //Put in final data structures and return:
        res.cumulCellContrib.reserve(tmp.size());
        res.cellInfo.reserve(tmp.size());
        StableSumKahan contribSum;
        for ( auto& ee : tmp ) {
          contribSum.add( ee.first );
          res.cumulCellContrib.push_back(contribSum.sum());
          res.cellInfo.push_back( ee.second );
        }
        return res;
      }

      SABProcImpl::SABProcImpl( const SABCfg::Cfg& cfg,
                                shared_obj<const SABData> sd,
                                std::shared_ptr<const VectD> requested_egrid,
                                SampleSupport sampleSupport,
                                StoreExtraDiagnostics extraDiag )
        : m_kT( sd->temperature().kT() ),
          m_invkT( 1.0 / m_kT ),
          m_cellmgr(sd)
      {
        const bool do_sample( sampleSupport == SampleSupport::YES );
        const bool do_diag( extraDiag == StoreExtraDiagnostics::YES );
        SABSurveyor surv(m_cellmgr.sabData());

        m_eGrid = determineEGrid( cfg, m_cellmgr.sabData(), m_cellmgr,
                                  surv, requested_egrid );

// #if 1 // fixme
//         m_eGrid.push_back(0.09*(1.0-1e-9));
//         m_eGrid.push_back(0.09);
//         m_eGrid.push_back(0.09*(1.0+1e-9));
//         std::sort(m_eGrid.begin(),m_eGrid.end());
// #endif

        const auto scheme = cfg.integScheme;

        std::size_t nTouchedCellsMax(0);
        {
          static_assert
            ( -static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::min())
              >= static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::max()),
              "" );
          constexpr auto n_cellmax
            = static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max());

          if ( surv.data().size() > n_cellmax )
            NCRYSTAL_THROW2(BadInput,"Too many SAB cells (does not fit in int32)");
          m_sIntegral.reserve( m_eGrid.size() );
          sIntegralAtE_Result integAtE;
          if ( do_sample ) {
            m_sampleIdx.reserve( m_eGrid.size() );
            integAtE.crossedIntegrals.reserve(256);
            m_bcEptInfo.reserve(128);
          }
          if ( do_diag )
            m_diagnostics_fcAR.reserve( m_eGrid.size() );
          for ( auto& E_div_kT : m_eGrid ) {
            sIntegralAtE( E_div_kT, m_cellmgr, surv, scheme, integAtE,
                          do_sample );
            m_sIntegral.push_back( integAtE.integralWithinKB );
            if ( do_diag )
              m_diagnostics_fcAR.push_back( integAtE.integralTouchedCells
                                            ? ( integAtE.integralWithinKB
                                                /integAtE.integralTouchedCells )
                                            : 1.0 );
            if (!do_sample)
              continue;

            const bool needsBoundedCellSampler
              = ( integAtE.integralWithinKB
                  < (cfg.fullCellSamplingARThreshold
                     *integAtE.integralTouchedCells) );
            nc_assert( integAtE.nTouchedCells < n_cellmax );
            nc_assert( integAtE.nTouchedCells > 0 );
            if ( needsBoundedCellSampler ) {
              //Needs bounded cell sampler.
              auto sidx = m_bcEptInfo.size();
              nc_assert(sidx<n_cellmax);
              m_sampleIdx.push_back(static_cast<std::int32_t>(sidx));
              m_bcEptInfo.push_back( prepareBCEPt( E_div_kT, surv, m_cellmgr,
                                                   integAtE,
                                                   cfg.integSchemeBCSample,
                                                   m_bcSample ) );
            } else {
              //Full cell sampling is OK, store -nTouchedCells.
              m_sampleIdx.push_back(-static_cast<std::int32_t>
                                    (integAtE.nTouchedCells));
              nTouchedCellsMax = std::max<std::size_t>(nTouchedCellsMax,
                                                       integAtE.nTouchedCells);

            }
          }
        }

        if ( do_sample) {
          auto itCellInfo = surv.data().begin();
          m_cumulFCInt.reserve(nTouchedCellsMax);
          m_cumulFCInt_cellidx.reserve(nTouchedCellsMax);
          StableSumKahan sum_cumulFCInt;
          for ( std::size_t i = 0; i < nTouchedCellsMax; ++i, ++itCellInfo ) {
            nc_assert(itCellInfo != surv.data().end());
            sum_cumulFCInt.add( m_cellmgr.getCellIntegral( itCellInfo->cellidx ) );
            m_cumulFCInt.push_back(sum_cumulFCInt.sum());
            m_cumulFCInt_cellidx.push_back(itCellInfo->cellidx);
          }
          nc_assert( m_cumulFCInt.size() == nTouchedCellsMax );
          nc_assert( m_cumulFCInt_cellidx.size() == nTouchedCellsMax );
          m_bcEptInfo.shrink_to_fit();
        }

      }

      double SABProcImpl::phaseSpaceIntegral( NeutronEnergy ekin ) const
      {
        //fixme: just a temporary implementation, we can certainly improve and
        //for instance cache some of the 1.0/(x1-x0) factors.
        const double t = ekin.dbl() * m_invkT;
        auto& x = m_eGrid;
        auto& y = m_sIntegral;
        nc_assert( x.size() == y.size() );
        nc_assert( x.size() >= 2 );
        nc_assert( x.front() >= 0.0 );
        nc_assert( t >= 0.0 );
        nc_assert( std::isfinite(t) );

        if (t >= x.back())
          return y.back();
        if (t <= x.front()) {
          const double kkk = y.front() * std::sqrt(1.0/x.front());//fixme: cache
          return kkk * std::sqrt(t);
        }
        auto it = std::lower_bound(x.begin(), x.end(), t);
        std::size_t i = static_cast<std::size_t>(it - x.begin()); // 1..n-1
        nc_assert( i>0 );
        nc_assert( i<x.size() );
        double x0 = x[i-1], x1 = x[i];
        double y0 = y[i-1], y1 = y[i];
        double r = (t - x0) / (x1 - x0);//fixme: cache more!
        return y0 * (1.0-r) + r * y1;
      }

      ScatterOutcomeIsotropic
      SABProcImpl::ab2scatout( RNG& rng,
                               NeutronEnergy ekin,
                               const AlphaBetaOutcome& alphabeta ) const
      {
        double dE, mu;
        if ( muIsotropicAtBeta(alphabeta.beta,ekin.get()*m_invkT) ) {//fixme: reuse E/kT
          //close to kinematical end-point, numerically safe fall-back:
          dE = alphabeta.beta*m_kT;
          mu = rng.generate()*2.0 - 1.0;
        } else {
          auto dEmu = convertAlphaBetaToDeltaEMu( alphabeta.alpha,
                                                  alphabeta.beta, ekin, m_kT );
          dE = dEmu.deltaE;
          mu = dEmu.mu;
        }
        return { NeutronEnergy{ncmax(0.0,ekin.dbl()+dE)}, CosineScatAngle{mu} };
      }



      SABProcessor::AlphaBetaOutcome
      SABProcImpl::sampleAlphaBeta( RNG& rng, NeutronEnergy ekin ) const
      {
        nc_assert( canSample() );
        //Find the index of the first value in the energy grid whose value exceeds
        //ekin (or the last value if ekin>=m_eGrid.back()):
        const double E_div_kT = ekin.dbl()*m_invkT;
        std::size_t idx_E_overlay;
        double E_div_kT_overlay;
        {
          nc_assert(m_eGrid.size()>=2);
          auto it = std::lower_bound(m_eGrid.begin(), m_eGrid.end(), E_div_kT);
          if ( it == m_eGrid.end() )
            it = std::prev(it);
          idx_E_overlay = static_cast<std::size_t>(it - m_eGrid.begin());
          E_div_kT_overlay = *it;
        }
        nc_assert(m_sampleIdx.size()==m_eGrid.size());
        std::int32_t sampleIdx = vectAt(m_sampleIdx,idx_E_overlay);
        ABRes alphabeta;
        const double foure = E_div_kT*4.0;
        auto phaseSpaceNotOK = [foure](const ABRes& ab )
        {
          return ( ncsquare( ab.a - ab.b ) > foure * ab.a );
        };
        std::size_t ntries(0);
        if ( sampleIdx < 0 ) {
          std::size_t nCellsTouched = static_cast<std::size_t>(-sampleIdx);
          alphabeta = sampleAlphaBetaFC( rng, E_div_kT, nCellsTouched );
          ntries += alphabeta.ntries;
          nc_assert( !phaseSpaceNotOK(alphabeta) );
        } else {
          const auto& bcEpt = vectAt(m_bcEptInfo,static_cast<std::size_t>(sampleIdx));
          do {
            alphabeta = sampleAlphaBetaBC( rng, E_div_kT_overlay, bcEpt );
            ntries += alphabeta.ntries;
          } while( phaseSpaceNotOK(alphabeta) );
        }
        SABProcessor::AlphaBetaOutcome res;
        res.alpha = alphabeta.a;
        res.beta = alphabeta.b;
        res.ntries = ntries;
        return res;
      }

      SABProcImpl::ABRes
      SABProcImpl::sampleAlphaBetaFC( RNG& rng,
                                      double E_div_kT,
                                      std::size_t nCellsTouched ) const
      {
        ABRes res;
        res.ntries = 1;
        const double foure = 4.0*E_div_kT;
        Span<const double> cumulContrib( m_cumulFCInt.data(),
                                         m_cumulFCInt.data()+nCellsTouched );
        //fixme: possible optimisation: cache last few (idx,FullCellSampler)
        //objects, in case a few cells are hit often?
        while ( true ) {
          std::size_t randidx = pickRandIdxByWeight( rng, cumulContrib );
          auto cellidx = vectAt(m_cumulFCInt_cellidx,randidx);
          auto cellData = m_cellmgr.lookupCellInfo( cellidx );
          auto ab = FullCellSampler::sampleOneAlphaBeta(cellData,rng);
          if ( ncsquare(ab.first-ab.second)<=foure*ab.first ) {
            res.a = ab.first;
            res.b = ab.second;
            return res;
          }
          ++res.ntries;
          nc_assert(res.ntries<100000);
        }
      }

      SABProcImpl::ABRes
      SABProcImpl::sampleAlphaBetaBC( RNG& rng,
                                      double E_div_kT,
                                      const BCEnergyPoint& bce ) const
      {
        std::size_t idx = pickRandIdxByWeight( rng, bce.cumulCellContrib );
        const auto& bce_ci = vectAt(bce.cellInfo,idx);
        CellData cellData = m_cellmgr.lookupCellInfo( bce_ci.sabCellIndex );
        const double foure = 4.0*E_div_kT;
        ABRes res;
        res.ntries = 0;
        if ( bce_ci.bcSampleInfoIndex
             == std::numeric_limits<std::uint32_t>::max() ) {
          //Sample in cell with full cell overlay:
          FullCellSampler fc(&cellData);
          while ( true ) {
            auto ab = fc.sampleAlphaBeta( rng );
            ++res.ntries;
            if ( ncsquare(ab.first-ab.second)<=foure*ab.first ) {
              res.a = ab.first;
              res.b = ab.second;
              return res;
            }
            nc_assert(res.ntries<100000);
          }
        } else {
          //Sample in cell with bounded cell overlay:
          auto bcsdata = BoundedCellSampler::unpack( m_bcSample,
                                                     bce_ci.bcSampleInfoIndex );
          BoundedCellSampler bc( cellData, bcsdata, E_div_kT );
          while ( true ) {
            auto ab = bc.sampleAlphaBeta( rng );
            res.ntries += ab.ntries;
            if ( ncsquare(ab.alpha-ab.beta)<=foure*ab.alpha ) {
              res.a = ab.alpha;
              res.b = ab.beta;
              return res;
            }
            nc_assert(res.ntries<100000);
          }
        }
      }

      void SABProcImpl::toJSON( std::ostream& os ) const
      {
        //Lots of diagnostics: Things related to egrid search (like how many
        //calls needed for root finding, how emin/emax were determined, how many
        //points added below emin ).
        //
        //Things related to: sampling method at each egrid

        nc_assert_always( m_eGrid.size()
                          < std::numeric_limits<std::uint32_t>::max() );
        os <<"{\"egrid\":";
        streamJSONHugeDblVect(os, VectD( m_eGrid ));
        os <<",\"sintegral\":";
        streamJSONHugeDblVect(os, VectD( m_sIntegral ));
        os <<",\"kT\":";
        streamJSON(os,m_kT);

        const bool can_sample = canSample();

        if (!m_diagnostics_fcAR.empty()) {
          os << ",\"diagnostics\":{\"full_cell_AR\":";
          streamJSON(os,m_diagnostics_fcAR);
          os << '}';
        }

        if ( can_sample ) {
          os << ",\"sample_info\":{";
          std::vector<std::uint32_t> epts_with_bc;
          epts_with_bc.reserve(512);
          nc_assert_always( m_sampleIdx.size() == m_eGrid.size() );
          const std::uint32_t n
            = static_cast<std::uint32_t>(m_sampleIdx.size());
          for ( std::uint32_t i = 0; i < n; ++i ) {
            if ( vectAt(m_sampleIdx,i) >= 0 )
              epts_with_bc.push_back(i);
          }
          os << "\"epts_with_bounded_cell_sampling\":";
          streamJSON(os,epts_with_bc);
          os << '}';
        }
        os << '}';
      }

      inline const SABProcImpl* sp_cimpl( const void* implptr ) noexcept
      {
        assert(implptr!=nullptr);//no exceptions
        return static_cast<const SABProcImpl*>(implptr);
      }

    }//anon namespace
  }
}

NCS::SABProcessor::~SABProcessor()
{
  if ( m_impl ) {
    SABProcImpl* sp = static_cast<SABProcImpl*>(m_impl);
    delete sp;
    m_impl = nullptr;
  }
}


double NCS::SABProcessor::phaseSpaceIntegral( NeutronEnergy ekin ) const
{
  return sp_cimpl(m_impl)->phaseSpaceIntegral(ekin);
}

NC::CrossSect NCS::SABProcessor::crossSectionUnitSigmaBound( NeutronEnergy ekin ) const
{
  //Fixme: optimise (most likely by caching in new vector in which we can
  //directly interpolate results as per the old SABXSProvider?)
  //fixme: keep synchronized with getEMaxInfo()
  auto sp = sp_cimpl(m_impl);
  const double Sint = sp->phaseSpaceIntegral(ekin);
  const double kTdiv4 = sp->m_kT * 0.25;//fixme: cache?
  return CrossSect{ Sint * kTdiv4 / ekin.dbl() };
}

NCS::SABProcessor::EPtInfo NCS::SABProcessor::getEMaxInfo() const
{
  auto sp = sp_cimpl(m_impl);
  EPtInfo res;
  res.E_div_kT = sp->m_eGrid.back();
  res.phaseSpaceIntegral = sp_cimpl(m_impl)->m_sIntegral.back();
  res.ekin = NeutronEnergy{ res.E_div_kT * sp->m_kT };
  res.crossSectionUnitSigmaBound
    = CrossSect{ res.phaseSpaceIntegral * 0.25 / res.E_div_kT };
  return res;
}

NC::ScatterOutcomeIsotropic
NCS::SABProcessor::sampleScatter( RNG& rng, NeutronEnergy ekin ) const
{
  auto alphabeta = sp_cimpl(m_impl)->sampleAlphaBeta(rng,ekin);
  return sp_cimpl(m_impl)->ab2scatout( rng, ekin, alphabeta );
}

NCS::SABProcessor::ScatOutcomeDiag
NCS::SABProcessor::sampleScatterDiag( RNG& rng, NeutronEnergy ekin ) const
{
  ScatOutcomeDiag res;
  res.alphabeta = sp_cimpl(m_impl)->sampleAlphaBeta(rng,ekin);
  res.ekinmu = sp_cimpl(m_impl)->ab2scatout( rng, ekin, res.alphabeta );
  return res;
}

NCS::SABProcessor::AlphaBetaOutcome
NCS::SABProcessor::sampleScatterAlphaBeta( RNG& rng, NeutronEnergy ekin ) const
{
  return sp_cimpl(m_impl)->sampleAlphaBeta(rng,ekin);
}

NCS::SABProcessor::SABProcessor( const SABCfg::Cfg& cfg,
                                 shared_obj<const SABData> sd,
                                 std::shared_ptr<const VectD> req_egrid,
                                 SampleSupport sampleSupport,
                                 StoreExtraDiagnostics extraDiag )
  : m_impl( static_cast<void*>( new SABProcImpl( cfg,
                                                 std::move(sd),
                                                 std::move(req_egrid),
                                                 sampleSupport,
                                                 extraDiag ) ) )
{
}

NCS::SABProcessor::SABProcessor( SABProcessor&& o ) noexcept
{
  std::swap(m_impl,o.m_impl);
}

NCS::SABProcessor& NCS::SABProcessor::operator=( SABProcessor&& o ) noexcept
{
  std::swap(m_impl,o.m_impl);
  return *this;
}

void NCS::SABProcessor::toJSON( std::ostream& os ) const
{
  return sp_cimpl(m_impl)->toJSON(os);
}

bool NCS::SABProcessor::hasSampleSupport() const
{
  return sp_cimpl(m_impl)->canSample();
}

double NCS::SABProcessor::kT() const
{
  return sp_cimpl(m_impl)->m_kT;
}

const NC::VectD& NCS::SABProcessor::getEDivKTGrid() const
{
  return sp_cimpl(m_impl)->m_eGrid;
}

const NC::VectD& NCS::SABProcessor::getPhaseSpaceIntegralAtGrid() const
{
  return sp_cimpl(m_impl)->m_sIntegral;
}

NC::shared_obj<const NC::SABData> NCS::SABProcessor::sabDataPtr() const
{
  return sp_cimpl(m_impl)->m_cellmgr.sabDataPtr();
}

void NCS::SABProcessor::toJSONProcessInfo( std::ostream& os,
                                           Optional<SigmaBound> sigma_scale,
                                           Optional<std::string>
                                           extension_method ) const
{
  const auto sp = sp_cimpl(m_impl);
  const auto& sab = sp->m_cellmgr.sabData();
  const auto emax = getEMaxInfo().ekin;
  const auto negrid = sp->m_eGrid.size();
  const auto emin = NeutronEnergy{ sp->m_eGrid.front()*sp->m_kT };
  {
    std::ostringstream tmp;
    tmp << "nalpha="<<sab.alphaGrid().size()<<";nbeta="<<sab.betaGrid().size();
    tmp << ";Emax="<<emax;
    tmp << ";T="<<sab.temperature();
    tmp << ";M="<<sab.elementMassAMU();
    if ( extension_method.has_value() && extension_method.value() != "freegas")
      tmp << ";extend="<<extension_method.value();
    if ( sigma_scale.has_value() )
      tmp << ";sigma_free="
          <<sigma_scale.value().free(sab.elementMassAMU());

    streamJSONDictEntry( os, "summarystr", tmp.str(), JSONDictPos::FIRST );
  }
  streamJSONDictEntry( os, "Emax", emax.dbl()  );
  streamJSONDictEntry( os, "Emin", emin.dbl()  );
  streamJSONDictEntry( os, "negrid", negrid  );
  streamJSONDictEntry( os, "T", sab.temperature().dbl()  );
  streamJSONDictEntry( os, "M", sab.elementMassAMU().dbl()  );
  if ( sigma_scale.has_value() ) {
    streamJSONDictEntry( os, "sigma_bound=",
                         sigma_scale.value().dbl() );
    streamJSONDictEntry( os, "sigma_free=",
                         sigma_scale.value().free(sab.elementMassAMU()).dbl() );
  }
  streamJSONDictEntry( os, "nalpha", sab.alphaGrid().size()  );
  streamJSONDictEntry( os, "nbeta", sab.betaGrid().size()  );
  if ( extension_method.has_value() )
    streamJSONDictEntry( os, "extension_method", extension_method.value() );
  streamJSONDictEntry( os, "sabprocessor_uid",
                       getUniqueID().value, JSONDictPos::LAST );
}

//fixme: remove:
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/phys_utils/NCFreeGasUtils.hh"
void NCS::SABProcessor::testJSON( shared_obj<const SABData> sd,
                                  std::ostream& os )
{
  //MAKE touchedinteg(E) curve first -> reeeelatively cheap I hope, and we
  //anyway need (most) of the cell integrals.
  SABSurveyor surv(sd);
  CellMgr cellmgr(sd);
  nc_assert_always( !surv.data().empty() );

  auto cfg = SABCfg::createConfig( 3 );

  //fixme: some backwards compat overrides:
  cfg.integScheme = SABCfg::IntegrationScheme::Flex17;
  cfg.integSchemeDetermineEGrid = SABCfg::IntegrationScheme::Flex5;
  cfg.egrid_npts = 300;

  VectD test_egrid = determineEGrid( cfg, sd, cellmgr, surv, {} );
  VectD test_sint;
  VectD test_sint_full;
  std::vector<std::size_t> test_sint_ncrossed;
  const auto scheme = cfg.integScheme;
  sIntegralAtE_Result result;
  result.crossedIntegrals.reserve(256);
  for ( auto& E_div_kT : test_egrid ) {
    sIntegralAtE( E_div_kT, cellmgr, surv, scheme, result, true );
    test_sint.push_back( result.integralWithinKB );
    test_sint_full.push_back( result.integralTouchedCells );
    test_sint_ncrossed.push_back( result.crossedIntegrals.size() );
  }

  auto sIntByEtouch
    = determineSIntegralByTouchedCells( cellmgr, surv,
                                        sd->temperature().kT() );


  os << "{\"E/kT\":";
  streamJSON(os,sIntByEtouch.e);
  os << ",\"sumcellint\":";
  streamJSON(os,sIntByEtouch.sint);
  os << ",\"kT\":";
  streamJSON(os,sd->temperature().kT());
  os << ",\"T\":";
  streamJSON(os,sd->temperature().get());
  os << ",\"E/kT_careful\":";
  streamJSON(os,test_egrid);
  os << ",\"Sintegral_careful\":";
  streamJSON(os,test_sint);
  os << ",\"Sintegral_careful_fullcells\":";
  streamJSON(os,test_sint_full);
  os << ",\"Sintegral_careful_ncrossed\":";
  streamJSON(os,test_sint_ncrossed);
  os << ",\"egrid_range\":";
  std::pair<double,double> egrid_range{ test_egrid.front(),
                                        test_egrid.back() };
  streamJSON(os,egrid_range);
  os << "}";

}

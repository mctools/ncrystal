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

#include "NCrystal/internal/sab/NCSABUCN.hh"
#include "NCrystal/internal/sab/NCSABEval.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/fact_utils/NCFactoryUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace UCN {

    std::pair<std::vector<StableSum>,VectD> UCNHelper::getSIntegralsAndOverlayVals( const SABData& sabData,
                                                                                    NeutronEnergy ucn_threshold,
                                                                                    const VectD& eVals )
    {
      nc_assert_always(ucn_threshold.dbl()>0.0);
      const auto& sab = sabData.sab();
      const auto& agrid = sabData.alphaGrid();
      const auto& bgrid = sabData.betaGrid();
      auto nalpha = agrid.size();
      const double kT = sabData.temperature().kT();
      const double inv_kT = 1.0 / kT;
      const double beta_ucn = ucn_threshold.dbl() * inv_kT;

      nc_assert(nc_is_grid(eVals));
      nc_assert_always(eVals.front()>0.0);

      const std::size_t n = eVals.size();

      auto eVals_div_kT = vectorTrf( eVals, [inv_kT](double x) { return x*inv_kT; } );

      std::vector<StableSum> Sintegrals;
      VectD overlayVals;
      Sintegrals.resize(eVals.size());
      overlayVals.resize(eVals.size(),0.0);

      constexpr double overlay_safety = 1.2;//safety (since egrid has limited granularity)

      auto itB = bgrid.begin();
      auto itBEnd = bgrid.end();
      double beta_prev = *itB;
      ++itB;
      decltype(nalpha) row0offset_next = 0;
      for ( ; itB != itBEnd; ++itB ) {
        if ( beta_prev > beta_ucn )
          break;
        PairDD beta{ beta_prev, *itB };
        beta_prev = beta.second;

        const auto row0offset = row0offset_next;
        row0offset_next += nalpha;
        const auto row1offset = row0offset + nalpha;

        //Given beta range and ucn threshold, we can figure out what ranges of
        //energy and alpha to consider:
        const double emax_div_kT = ncmax( 0.0, beta_ucn - beta.first );
        const double emin_div_kT = ncmax( 0.0, -beta.second );

        double effective_b1 = ncmin( beta.second, beta_ucn - emin_div_kT );
        double effective_b0 = ncmax( beta.first, beta_ucn - emax_div_kT );
        if ( !(effective_b1>effective_b0) )
          continue;
        const double alpha_min = ( effective_b1 < 0.0
                                   ? getAlphaMinus( emax_div_kT, effective_b1 )
                                   : ( effective_b0 > 0
                                       ? getAlphaMinus( emax_div_kT, effective_b0 )
                                       : 0.0 ) );
        const double alpha_max = ( effective_b1 <= -emax_div_kT
                                   ? -999.0
                                   : getAlphaPlus( emax_div_kT, ncmin( beta_ucn - emax_div_kT, effective_b1 ) ) );
        if ( alpha_max <= agrid.front() || alpha_min >= agrid.back() )
          continue;

        if ( eVals_div_kT.back() <= emin_div_kT )
          continue;
        if ( eVals_div_kT.front() >= emax_div_kT )
          continue;

        auto itEkinMin = std::lower_bound(eVals_div_kT.begin(),eVals_div_kT.end(),emin_div_kT);
        if ( itEkinMin != eVals_div_kT.begin() )
          itEkinMin = std::prev(itEkinMin);
        const std::size_t iekin_min = static_cast<std::size_t>( std::distance(eVals_div_kT.begin(),itEkinMin) );

        std::size_t iekin_max = iekin_min;
        while( iekin_max+1 < n && vectAt(eVals_div_kT,iekin_max+1) <= emax_div_kT )
          ++iekin_max;
        ++iekin_max;
        if ( !(iekin_min<n) )
          continue;
        nc_assert_always(iekin_max<=n);
        nc_assert_always(iekin_min<iekin_max);

        for ( auto ialpha_low : ncrange(nalpha-1) ) {
          auto itA0 = std::next(agrid.begin(), ialpha_low );
          if ( *itA0 >= alpha_max )
            break;

          PairDD alpha{ *itA0, *(std::next(itA0)) };
          if ( alpha.second <= alpha_min )
            continue;
          const double svals[4] = { vectAt(sab, row0offset + ialpha_low ),
                                    vectAt(sab, row0offset + ialpha_low + 1 ),
                                    vectAt(sab, row1offset + ialpha_low),
                                    vectAt(sab, row1offset + ialpha_low + 1 ) };
          SABUtils::SABCellEval<> cell{ alpha, beta, svals };//todo: we could in principle reuse two of the log values.
          for ( std::size_t iekin = iekin_min; iekin<iekin_max; ++iekin ) {
            const double ekin_div_kT = vectAt(eVals_div_kT,iekin);
            if ( beta.second <= -ekin_div_kT )
              continue;
            if ( beta.first >= -ekin_div_kT + beta_ucn )
              continue;
            const double b0_m_a1 = beta.first - alpha.second;
            if ( b0_m_a1 >= 0.0 && ncsquare(b0_m_a1) >= 4.0*ekin_div_kT*alpha.second )
              continue;
            if ( alpha.first >= beta.second ) {
              if ( alpha.second <= ekin_div_kT ) {
                if ( ncsquare( alpha.second - beta.second ) >= 4.0*ekin_div_kT*alpha.second )
                  continue;
              } else if ( alpha.first >= ekin_div_kT ) {
                if ( ncsquare( alpha.first - beta.second ) >= 4.0*ekin_div_kT*alpha.first )
                  continue;
              }
            }
            double& soverlay_val = vectAt(overlayVals,iekin);
            soverlay_val = ncmax( soverlay_val,
                                  overlay_safety * cell.sOverlayValueWithinKinematicBoundsBelowBetamax(ekin_div_kT,beta_ucn-ekin_div_kT) );
            cell.addIntegralWithinKinematicBoundsBelowBetamax( vectAt(Sintegrals,iekin),
                                                               ekin_div_kT,
                                                               beta_ucn-ekin_div_kT );
          }
        }
      }
      return { std::move(Sintegrals), std::move(overlayVals) };
    }



    namespace {
      VectD setupEGrid( const SABData& sabData,
                        NeutronEnergy ucn_threshold )
      {

        nc_assert_always(ucn_threshold.dbl()>0.0);
        const auto& agrid = sabData.alphaGrid();
        const auto& bgrid = sabData.betaGrid();
        const double kT = sabData.temperature().kT();
        const double inv_kT = 1.0 / kT;
        PairDD egrid_range( ucn_threshold.dbl(),
                            ncmax( ucn_threshold.dbl()*1e2, 1e2 ) );

        //Figure out emin/emax based on reach of grid
        const double eu = ucn_threshold.dbl()*inv_kT;
        const double amin = agrid.front() / eu;
        if ( amin > 1.0 ) {
          //There is a minimum energy, below which XS will be strictly 0 since
          //we use S=0 for alpha<agrid.front().
          egrid_range.first = ncmax( ucn_threshold.dbl() * ( 1.0 + amin - 2 * std::sqrt( amin ) ),
                                     egrid_range.first );
        }

        if ( ! ( agrid.back() > eu ) )
          NCRYSTAL_THROW2(BadInput,"UCN threshold ("<<ucn_threshold<<") is too high for provided "
                          "kernel (alpha-grid reach imposes a limit of "<<NeutronEnergy{kT*agrid.back()}<<")");
        const double amax = agrid.back() / eu;
        const double emax_a = ucn_threshold.dbl() * ( 1.0 + amax - 2 * std::sqrt( amax ) );
        nc_assert_always( bgrid.front() < 0.0 );
        const double emax_b = - kT * bgrid.front();
        egrid_range.second = ncmin(emax_a,emax_b);
        if ( ! ( egrid_range.second > egrid_range.first ) )
          NCRYSTAL_THROW2(BadInput,"UCNHelper could not determine egrid reach. Perhaps UCN threshold ("
                          <<ucn_threshold<<") is too high for provided kernel?");

        constexpr auto ngeomspace = 3000;
        constexpr auto nwl = 500;

        auto eVals0 = geomspace( egrid_range.first, egrid_range.second,ngeomspace );
        //Add some points linearly in wavelength as well:
        auto addWlPtsLinearly = [&eVals0,egrid_range](double wlmin, double wlmax, std::size_t nnn )
        {
          eVals0.reserve(eVals0.size()+nnn);
          for ( auto wl : linspace(wlmin,wlmax,nnn) ) {
            auto ee = wl2ekin(wl);
            if ( ee > egrid_range.first && ee < egrid_range.second )
              eVals0.push_back(ee);
          }
        };
        addWlPtsLinearly( 0.1, 12, nwl/2 );
        addWlPtsLinearly( 12, 300, nwl/2 );
        //And the pts indicated by the (negative) beta grid:
        for ( auto beta : bgrid ) {
          auto ee = -kT * beta;
          if ( ee > 0.0 && ee > egrid_range.first && ee < egrid_range.second )
            eVals0.push_back(ee);
        }
        std::sort(eVals0.begin(),eVals0.end());
        eVals0.erase( std::unique( eVals0.begin(), eVals0.end() ), eVals0.end() );
        eVals0.shrink_to_fit();
        return eVals0;
      }

      struct XSAndOverlay {
        PiecewiseLinearFct1D xs, overlay;
        double worst_AR = 0.0;
        double average_AR = 0.0;
        XSAndOverlay(PiecewiseLinearFct1D xx,PiecewiseLinearFct1D oo, double war, double aar)
          :xs(std::move(xx)), overlay(std::move(oo)), worst_AR(war), average_AR(aar)
        {}//for C++11 (annoying!)

      };

      XSAndOverlay extractXSAndOverlayCurves( const SABData& sabData,
                                              NeutronEnergy ucn_threshold,
                                              VectD eVals )
      {
        const std::size_t n = eVals.size();
        const auto& sab = sabData.sab();
        const double kT = sabData.temperature().kT();
        const double inv_kT = 1.0 / kT;

        std::vector<StableSum> Sintegrals;
        VectD overlayVals;
        std::tie(Sintegrals,overlayVals) = UCNHelper::getSIntegralsAndOverlayVals( sabData, ucn_threshold, eVals );

        VectD xsVals;
        xsVals.reserve(eVals.size());
        const double scale_fact = 0.25 * sabData.boundXS().get() * kT;
        const double globalOverlay = *std::max_element(sab.begin(), sab.end());
        nc_assert( std::max_element(sab.begin(), sab.end()) != sab.end() );

        if ( !scale_fact ) {
          //special signature for null scatter (likely boundXS was 0)
          return XSAndOverlay{ PiecewiseLinearFct1D{ {0.0,1e-199}, {0.0,0.0}, {0.0,0.0} },
                               PiecewiseLinearFct1D{ {0.0,1e-199}, {globalOverlay,globalOverlay}, {globalOverlay,globalOverlay} },
                   0.0, 0.0 };
        }

        const double eu = ucn_threshold.dbl()*inv_kT;
        const double Sarea_factor = (8.0/3.0)*eu*std::sqrt(eu);
        double worst_AR = 1.0;
        double worst_AR_e = -1.0;
        StableSum tot_AR;

        for ( auto i : ncrange(n) ) {
          const double sint = vectAt(Sintegrals,i).sum();
          const double eval = vectAt(eVals,i);
          xsVals.push_back( scale_fact * ( sint / eval ) );
          //Estimate sampling acceptance rate (because an error at init time is
          //preferable to runtime slow sampling/warnings/errors):
          const double ol = vectAt(overlayVals,i);
          const double Sarea = Sarea_factor * std::sqrt(eval*inv_kT);
          const double ols = ol*Sarea;
          const double acceptance_rate = ( ols ? sint/ols : 1.0);
          tot_AR.add( acceptance_rate );
          if ( acceptance_rate < worst_AR ) {
            worst_AR_e = eval;
            worst_AR = acceptance_rate;
          }
        }

        const double worst_allowed_AR = ncgetenv_dbl("UCNHELPER_WORST_AR",1e-3);
        if ( worst_AR < worst_allowed_AR ) {
          NCRYSTAL_THROW2(BadInput,"UCNHelper worst scatter acceptance rate is "<<worst_AR*100
                          <<"% (at "<<NeutronEnergy{worst_AR_e}<<") which is less than the"
                          " allowed worst case of "<<worst_allowed_AR*100
                          <<"%. Most likely the chosen UCN threshold ("<<ucn_threshold
                          <<") is too high. Alternatively (if you are willing to live with"
                          " very slow sampling) lower the limit via the environment variable"
                          " NCRYSTAL_UCNHELPER_WORST_AR.");
        }

        //Peel off 0 cross sections:
        while ( xsVals.size() > 2 && xsVals.back() == 0.0 ) {
          xsVals.pop_back();
          eVals.pop_back();
          overlayVals.pop_back();
        }

        std::size_t npopfront(0);
        for ( auto xs : xsVals ) {
          if ( xs )
            break;
          ++npopfront;
        }
        if ( npopfront > xsVals.size() - 2 )
          npopfront = xsVals.size() - 2;
        if ( npopfront ) {
          xsVals      = VectD( std::next( xsVals.begin(), npopfront ), xsVals.end() );
          eVals       = VectD( std::next( eVals.begin(), npopfront ), eVals.end() );
          overlayVals = VectD( std::next( overlayVals.begin(), npopfront ), overlayVals.end() );
        }
        xsVals.shrink_to_fit();
        eVals.shrink_to_fit();
        overlayVals.shrink_to_fit();

        return { PiecewiseLinearFct1D( eVals, std::move(xsVals), { 0.0, 0.0} ),
                 PiecewiseLinearFct1D( eVals, std::move(overlayVals), {globalOverlay,globalOverlay} ),
                 worst_AR, tot_AR.sum() / n };
      }

      XSAndOverlay extractXSAndOverlayCurves( const SABData& sabData,
                                              NeutronEnergy ucn_threshold )
      {
        nc_assert_always(ucn_threshold.dbl()>0.0);
        VectD eVals = setupEGrid( sabData, ucn_threshold );
        return extractXSAndOverlayCurves( sabData, ucn_threshold, std::move(eVals) );
      }
    }

    UCNHelper::UCNHelper( shared_obj<const SABData> sabData, NeutronEnergy ucn_threshold )
      : m_data([&sabData,ucn_threshold]()
      {
        auto curves = extractXSAndOverlayCurves( sabData, ucn_threshold );
        return Data{ std::move(curves.xs),
                     std::move(curves.overlay),
                     ucn_threshold,
                     std::move(sabData),
                     curves.worst_AR,
                     curves.average_AR };
      }())
    {
    }

    ScatterOutcomeIsotropic UCNHelper::sampleScatterIsotropic( RNG& rng, NeutronEnergy ekin ) const
    {
      const double emin = m_data.xs.xValues().front();
      if ( isNull() || ekin.dbl() > m_data.xs.xValues().back() || ekin.dbl() < emin )
        return ScatterOutcomeIsotropic::noScat( ekin );//scattering not possible
      const double overlay = m_data.overlay(ekin.dbl());
      if (!(overlay>0.0))
        return ScatterOutcomeIsotropic::noScat( ekin );//scattering not possible (should only happen if xs was vanishing...)

      SABUtils::SABEval<> sabEvaluator(m_data.m_sabData);
      const double kT = m_data.m_sabData->temperature().kT();
      const double inv_kT = 1.0 / kT;
      const double ekin_div_kT = ekin.dbl() * inv_kT;
      const double eucn_div_kT = m_data.ucnthr.dbl() * inv_kT;

      auto genFlatAlphaBeta = [&rng,ekin_div_kT,eucn_div_kT]() -> PairDD
      {
        //Flat in alpha/beta, within kinematic bounds (including efinal<=eucn_div_kT).
        const double beta = -ekin_div_kT + eucn_div_kT * std::cbrt( ncsquare( rng.generate() ) );
        auto alim = getAlphaLimitsWithDiff( ekin_div_kT, beta );
        const double alpha = ncclamp( alim.aminus + alim.adiff * rng.generate(),
                                      alim.aminus, alim.aplus );
        return { alpha, beta };
      };

      auto toOutcome = [ekin,kT,ekin_div_kT,&rng]( const PairDD& alphabeta ) -> ScatterOutcomeIsotropic
      {

        const auto& beta = alphabeta.second;
        if ( beta <= -ekin_div_kT || muIsotropicAtBeta( beta, ekin_div_kT ) ) {
          nc_assert( beta >= -ekin_div_kT * 1.001 );
          //close enough to kinematical end-point to make numerical trouble.
          return { NeutronEnergy{ ncmax( 0.0, ekin.dbl()+beta * kT ) },
                   randIsotropicScatterMu(rng) };
        } else {
          auto demu = convertAlphaBetaToDeltaEMu( alphabeta, ekin, kT );
          return { NeutronEnergy{ ncmax(0.0,ekin.dbl() + demu.deltaE) },
                   CosineScatAngle{demu.mu} };
        }
      };

      unsigned long ii{0};
      while ( true ) {
        PairDD ab = genFlatAlphaBeta();
        const double S_at_ab = sabEvaluator.eval( ab );
        if ( S_at_ab > overlay ) {
          NCRYSTAL_WARN("UCNHelper: Overlay value in sampling too small"
                        " by at least a factor: "<<S_at_ab/overlay
                        <<" (please report to NCrystal developers)");
        }
        if ( S_at_ab >= overlay * rng.generate() )
          return toOutcome(ab);
        if ( ++ii == 1000000000 ) {
          NCRYSTAL_THROW(CalcError,"UCNHelper: Sampling is taking ridiculously long (infinite loop?).");
        }
      }
    }

    UCNScatter::UCNScatter( shared_obj<const SABData> data, NeutronEnergy ucn_threshold )
      : m_helper(std::move(data),ucn_threshold)
    {
    }

    CrossSect UCNScatter::crossSectionIsotropic( CachePtr&, NeutronEnergy ekin ) const
    {
      return m_helper.crossSection( ekin );
    }

    ScatterOutcomeIsotropic UCNScatter::sampleScatterIsotropic( CachePtr&, RNG& rng, NeutronEnergy ekin ) const
    {
      return m_helper.sampleScatterIsotropic( rng, ekin );
    }

    EnergyDomain UCNHelper::domain() const noexcept
    {
      if ( isNull() )
        return EnergyDomain::null();
      const auto& x = m_data.xs.xValues();
      return { NeutronEnergy{ x.front() }, NeutronEnergy{ x.back() } };
    }

    EnergyDomain UCNScatter::domain() const noexcept
    {
      return m_helper.domain();
    }

    Optional<std::string> UCNScatter::specificJSONDescription() const
    {
      CrossSect maxXS {-1.0};
      NeutronEnergy maxXS_E{-1.0};
      auto xscurve = m_helper.accessXSCurve();
      for ( auto iE : ncrange(xscurve.size()) ) {
        const auto probeE = NeutronEnergy{ vectAt(xscurve.xValues(),iE) };
        const auto probeXS = CrossSect{ vectAt(xscurve.yValues(),iE) };
        if ( maxXS < probeXS ) {
          maxXS   = probeXS;
          maxXS_E = probeE;
        }
      }
      const auto ucnthr = m_helper.ucnThreshold();
      nc_assert_always(maxXS.dbl()>=0.0);
      std::ostringstream ss;
      {
        std::ostringstream tmp;
        tmp << "Eucn="<<ucnthr;
        tmp << ";max_xs="<<fmt(maxXS.dbl(),"%.3g")<<maxXS.unit();
        tmp << "@"<<maxXS_E;
        tmp << ";avgAR=" << fmt(m_helper.averageSamplingAcceptanceRate(),  "%.2g");
        tmp << ";worstAR=" << fmt(m_helper.worstSamplingAcceptanceRate(),  "%.2g");
        tmp << ";xs_grid_E=" << fmt(xscurve.xValues().front(),"%.3g")
            << ".." << fmt(xscurve.xValues().back(),"%.3g") << NeutronEnergy::unit();
        streamJSONDictEntry( ss, "summarystr", tmp.str(), JSONDictPos::FIRST );
      }
      streamJSONDictEntry( ss, "ucn_threshold", ucnthr.dbl() );
      streamJSONDictEntry( ss, "max_xs", maxXS.dbl() );
      streamJSONDictEntry( ss, "max_xs_at_E", maxXS_E.dbl() );
      streamJSONDictEntry( ss, "worst_AR",m_helper.worstSamplingAcceptanceRate() );
      streamJSONDictEntry( ss, "average_AR",m_helper.averageSamplingAcceptanceRate() );
      streamJSONDictEntry( ss, "xs_grid_n", xscurve.size() );
      streamJSONDictEntry( ss, "xs_grid_E0", xscurve.xValues().front() );
      streamJSONDictEntry( ss, "xs_grid_Emax", xscurve.xValues().back(), JSONDictPos::LAST );
      return ss.str();
    }

    ExcludeUCNScatter::ExcludeUCNScatter( shared_obj<const ProcImpl::ScatterIsotropicMat> ws,
                                          shared_obj<const UCNScatter> us )
      : m_wrappedScatter( std::move(ws) ),
        m_ucnScatter( std::move(us) ),
        m_ucnDomain(m_ucnScatter->domain())
    {
    }

    Optional<std::string> ExcludeUCNScatter::specificJSONDescription() const
    {
      std::ostringstream ss;
      ss << "{\"components\":[";
      ss << "[1.0," <<m_wrappedScatter->jsonDescription()
         << "],[1.0,"<<m_ucnScatter->jsonDescription()<<"]]}";
      return ss.str();
    }

    CrossSect ExcludeUCNScatter::crossSectionIsotropic( CachePtr& cp, NeutronEnergy ekin ) const
    {
      auto xs = m_wrappedScatter->crossSectionIsotropic(cp,ekin);
      auto xs_ucn = m_ucnScatter->ucnHelper().crossSection( ekin );
      return CrossSect{ ncmax( 0.0, xs.dbl() - xs_ucn.dbl() ) };
    }

    namespace {

      static constexpr int detail_sampleScatterIsotropic_loopmax = 20;
      void detail_sampleScatterIsotropic_emit_loopmax_warning()
      {
        //This should hopefully never happen in practice, so we simply use a
        //static mutex to guard this section:
        static std::mutex mtx;
        NCRYSTAL_LOCK_GUARD(mtx);
        static unsigned nwarnings = 0;
        constexpr unsigned nwarnings_max = 50;
        if ( nwarnings < nwarnings_max ) {
          nwarnings += 1;
          NCRYSTAL_WARN("ExcludeUCNScatter: Wrapped process could not"
                        " sample non-UCN final state in "
                        << detail_sampleScatterIsotropic_loopmax
                        << " attempts!"<<
                        ( nwarnings == nwarnings_max
                          ? " (suppressing further WARNINGS of this type)"
                          : "" ));
        }
      }
    }

    ScatterOutcomeIsotropic ExcludeUCNScatter::sampleScatterIsotropic( CachePtr& cp, RNG& rng, NeutronEnergy ekin ) const
    {
      if ( !m_ucnDomain.contains(ekin) || m_ucnScatter->ucnHelper().crossSection(ekin).dbl() <= 0.0 ) {
        //No cross section of m_ucnScatter here, let m_wrappedScatter do its
        //thing competely unhindered:
        return m_wrappedScatter->sampleScatterIsotropic(cp,rng,ekin);
      }
      //m_ucnScatter supposedly provides all UCN events here, so we resample any
      //of those emitted by m_wrappedScatter (with infinite loop safeguards):
      auto ucnThreshold = m_ucnScatter->ucnHelper().ucnThreshold();
      int i = 0;
      while ( true ) {
        auto outcome = m_wrappedScatter->sampleScatterIsotropic(cp,rng,ekin);
        if  ( outcome.ekin >= ucnThreshold )
          return outcome;//the usual path, most often in first attempt.
        if ( ++i == detail_sampleScatterIsotropic_loopmax ) {
          detail_sampleScatterIsotropic_emit_loopmax_warning();
          return outcome;
        }
      }
    }

    namespace {
      struct UCNScatter_ThinnedKey {
        UniqueIDValue sabuid;
        ShortStrDbl ucnthr_str;
        bool operator<(const UCNScatter_ThinnedKey&o) const noexcept
        {
          if ( sabuid != o.sabuid )
            return sabuid < o.sabuid;
          return ucnthr_str.to_view() < o.ucnthr_str.to_view();
        }
      };
      struct UCNScatter_FullKey {
        UniqueIDValue sabuid;
        ShortStrDbl ucnthr_str;
        shared_obj<const SABData> sabData;
        UCNScatter_ThinnedKey thin() const { return {sabuid,ucnthr_str}; }
      };
      struct UCNScatter_KeyThinner {
        using key_type = UCNScatter_FullKey;
        using thinned_key_type = UCNScatter_ThinnedKey;
        template <class TMap>
        static typename TMap::mapped_type& cacheMapLookup( TMap& map, const key_type& key, Optional<thinned_key_type>& tkey )
        {
          if ( !tkey.has_value() )
            tkey = key.thin();
          return map[tkey.value()];
        }
      };
      using UCNScatPtr = shared_obj<const UCNScatter>;
      constexpr auto ucnscatfact_nstrongrefskept = 20;
      class UCNScatFact final : public CachedFactoryBase<UCNScatter_FullKey, UCNScatter, ucnscatfact_nstrongrefskept, UCNScatter_KeyThinner> {
      public:
        std::string keyToString( const UCNScatter_FullKey& key ) const override
        {
          std::ostringstream ss;
          ss << "UCNScatFactKey{sabuid:"<<key.sabuid.value<<",ucn_threshold:"<<key.ucnthr_str<<"}";
          return ss.str();
        }
        const char* factoryName() const override { return "UCNScatFact"; }
      protected:
        std::shared_ptr<const UCNScatter> actualCreate( const UCNScatter_FullKey& key ) const override
        {
          auto opt_ucnthrval = key.ucnthr_str.to_view().toDbl();
          nc_assert_always(opt_ucnthrval.has_value());
          return makeSO<UCNScatter>( key.sabData,
                                     NeutronEnergy{ opt_ucnthrval.value() } );
        }
      };
    }

    shared_obj<const UCNScatter> UCNScatter::createWithCache( shared_obj<const SABData> sabData,
                                                              NeutronEnergy ucn_threshold)
    {
      //Avoid FP numbers in keys - so (losslessly) format ucn_threshold into a
      //string for the key:
      auto key = UCNScatter_FullKey{ { sabData->getUniqueID().value },
                                     dbl2shortstr(ucn_threshold.dbl()),
                                     std::move(sabData) };
      static UCNScatFact s_db;
      return s_db.create(key);
    }
  }
}

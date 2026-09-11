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

#include "NCrystal/internal/vdos/NCVDOSExpand.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/vdos/NCVDOSUtils.hh"
#include "NCrystal/interfaces/NCSABData.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC=NCrystal;

NC::VDOS::GnExpansion
NC::VDOS::expandVDOSToGnFcts( const VDOSData& vdosdata,
                              VDOSLux vdoslux,
                              Optional<double> targetEmax_requested,
                              Optional<VDOSGn::Order> request_override_maxord )
{
  static bool s_verbose = ncgetenv_bool("DEBUG_PHONON");


  //Hidden unofficial env-vars used for special debugging purposes:
  const double override_alphamax = ncgetenv_dbl("HACK_ALPHAMAX");
  const double override_betamax = ncgetenv_dbl("HACK_BETAMAX");
  const unsigned override_max_order
    = ( request_override_maxord.has_value()
        ? request_override_maxord.value().value()
        : static_cast<unsigned>(ncgetenv_int("HACK_MAXORDER")) );

  //Which Emax should we target (i.e. aim to cover the kinematic reachable area
  //for neutrons of that energy):
  nc_assert_always(targetEmax_requested.value_or(1.0)>0.0);

  double targetEmax;//emax in eV
  if ( targetEmax_requested.has_value() ) {
    targetEmax = targetEmax_requested.value();
  } else {
    if ( vdoslux.isLegacy() ) {
      //used to depend on vdoslux.
      nc_assert_always( vdoslux.lvl() <= 5u );
      double legacy_lux2emax[6] = { 0.5, 1.0, 3.0, 5.0, 8.0, 12.0 };
      targetEmax = legacy_lux2emax[vdoslux.lvl()];
    } else {
      //always exactly 5.0 eV now.
      targetEmax = 5.0;
    }
  }

  if ( s_verbose )
    NCRYSTAL_MSG("VDOS expansion initialising with T="<<vdosdata.temperature()
                 <<", vdoslux="<<vdoslux.raw()
                 <<", aiming for Emax="<<targetEmax<<"eV"
                 <<(targetEmax_requested.has_value()
                    ?" (as requested)":"")<<", ...");

  //Initialise evaluators:
  VDOSEval vdoseval(vdosdata);
  const double kT = vdoseval.kT();
  const double invkT = 1.0/kT;
  const double gamma0 = vdoseval.calcGamma0();
  const double msd = vdoseval.getMSD( gamma0 );
  double targetEmax_div_kT = targetEmax*invkT;
  unsigned max_phonon_order = std::max<unsigned>(override_max_order,1);
  auto ttpars = ( vdoslux.isLegacy()
                  ? VDOSGn::TruncAndThinningChoices::Legacy
                  : VDOSGn::TruncAndThinningChoices::Default );

  constexpr double alpha2x_factor = ( 2.0*const_neutron_mass_evc2
                                      /(constant_hbar*constant_hbar) );
  const double alpha2x = alpha2x_factor*kT*msd;
  nc_assert_always( floateq( msd, alpha2x/(alpha2x_factor*kT) ) );

  GnExpansion res{ VDOSGn(vdoseval,ttpars), 0.0, 0.0, alpha2x };//, msd };
  auto& Gn_asym = res.Gn;

  Gn_asym.growMaxOrder(max_phonon_order);

  //What are the highest phonon order we allow?
  unsigned order_limit;
  if ( vdoslux.isLegacy() ) {
    order_limit = 1000;
    if ( targetEmax_requested.has_value() || vdoslux.lvl() == 5 )
      order_limit *= 10;
    if ( vdoslux.lvl()==0 )
      order_limit /= 10;
  } else {
    //fixme: this seems reasonable but when vdoslux=0 or 1 we might wish to
    //adjust the trunc and thinning parameters.
    order_limit = 1000;
    if ( targetEmax_requested.has_value() || vdoslux.lvl() >= 5 )
      order_limit = 10000;
    else
      order_limit = 1000;
  }

  const double emax_lowest_allowed = ( targetEmax_requested.value_or(1e-15) );

  //Now increase order dynamically until the last order only has contributions
  //to S(alpha,beta) outside the kinematic reach of Emax:
  const double relcontriblvl = [&vdoslux]()
  {
    if ( vdoslux.isLegacy() ) {
      //e.g.: 1e-3 for vdoslux 0, 1e-9 for vdoslux 3, 1e-13 for vdoslux 5
      return std::pow(10.0,-(3.0+2.0*vdoslux.lvl()));
    }
    switch ( vdoslux.lvl() ) {
    case 0: return 1e-6;
    case 1: return 1e-7;
    case 2: return 1e-8;
    default:
    case 3: return 1e-9;
    case 4: return 1e-11;
    case 5: return 1e-12;
    case 6: return 1e-13;
    }
  }();

  const double x2alpha = 1.0 / res.alpha2x;
  auto findAlphaBetaRangeOfOrder
    = [&Gn_asym,x2alpha,invkT,relcontriblvl](unsigned n)
    {
      auto eRange = Gn_asym.eRange(n, relcontriblvl);
      PairDD betaRange( eRange.first * invkT, eRange.second * invkT  );
      auto xRange = rangeXNexpMX( n, relcontriblvl );
      PairDD alphaRange( xRange.first * x2alpha, xRange.second * x2alpha  );
      return std::make_pair(alphaRange,betaRange);
    };

  while (true) {
    if (override_max_order>0)
      break;
    Gn_asym.growMaxOrder(max_phonon_order);
    PairDD alphaRange, betaRange;
    std::tie(alphaRange, betaRange)
      = findAlphaBetaRangeOfOrder(Gn_asym.maxOrder().value());
    if (sabPointWithinAlphaPlusCurve(targetEmax_div_kT,
                                     alphaRange.first,
                                     betaRange.second)) {
      //could consider larger stepsize, but need to carefully check usage in the
      //following.
      ++max_phonon_order;
    } else {
      break;
    }
    if (max_phonon_order>order_limit) {
      //Too slow - unfeasible to fill out S(alpha,beta) all the way out to the
      //kinematic curve for E=targetEmax. In this case it is better to reduce
      //targetEmax, to at least get a consistent table (and hope the free-gas
      //extrapolation mechanisms will be adequate already at this lower
      //threshold).
      double targetEmax_reduced  = targetEmax;
      do {
        targetEmax_reduced *= 0.99;
        if ( targetEmax_reduced < emax_lowest_allowed )
          NCRYSTAL_THROW2(CalcError,"VDOS expansion too slow - can not reach E="
                          <<emax_lowest_allowed<<"eV after "<<order_limit
                          <<" phonon convolutions (likely causes: either"
                          " the target energy value is too high, vdoslux too"
                          " low, the temperature too high, or the VDOS is"
                          " very unusual).");
      } while (sabPointWithinAlphaPlusCurve(targetEmax_reduced*invkT,
                                            alphaRange.first,
                                            betaRange.second));
      if (s_verbose)
        NCRYSTAL_WARN("VDOS expansion could only reach Emax="
                      <<targetEmax_reduced
                      <<"eV and not the requested Emax="<<targetEmax<<"K");
      targetEmax_div_kT = targetEmax_reduced * invkT;
      targetEmax = targetEmax_reduced;
      break;
    }
  }
  nc_assert_always( targetEmax_requested.value_or(targetEmax) == targetEmax );
  Gn_asym.growMaxOrder(max_phonon_order);

  //Ok, we now know how many orders we need to reach targetEmax. Next step is to
  //look at the contribution of each order insided the kinematic reach of
  //targetEmax, and use it to determine alpha/beta limits:

  if ( !(override_max_order>0) ) {
    //If maxorder was overridden, we can't know how far the kernel actually
    //reaches. Otherwise, provide the information:
    res.suggestedEmax = targetEmax;
  }

  double betaMin = 0.0;
  double alphaMax = 0.0;
  for ( unsigned n = 1; n<=max_phonon_order; ++n ) {
    PairDD alphaRange, betaRange;
    std::tie(alphaRange, betaRange) = findAlphaBetaRangeOfOrder(n);
    auto ep = findExtremeSABPointWithinAlphaPlusCurve(targetEmax_div_kT,
                                                      alphaRange, betaRange);
    if ( ep.has_value() ) {
      alphaMax = ncmax(alphaMax,ep.value().first);
      betaMin = ncmin(betaMin,ep.value().second);
    }
  }
  const auto reachmsg = ( "This can happen if temperature is too high. It"
                          " help to increase vdoslux if possible." );
  if ( !( betaMin < 0.0 ) )
    NCRYSTAL_THROW2(CalcError,"Beta range after "<<max_phonon_order
                    <<" phonon orders does not extend to negative beta. "
                    <<reachmsg);
  if ( !( alphaMax > 0.0 ) )
    NCRYSTAL_THROW2(CalcError,"Alpha range after "<<max_phonon_order
                    <<" phonon orders does not extend to positive alpha. "
                    <<reachmsg);
  if ( vdoslux.isLegacy() ) {
    res.betaUpper = -betaMin*1.01;
    res.alphaUpper = alphaMax*1.01;
  } else {
    res.betaUpper = -betaMin;
    res.alphaUpper = alphaMax;
  }
  if (override_alphamax)
    res.alphaUpper = override_alphamax;
  if (override_betamax)
    res.betaUpper = override_betamax;
  nc_assert_always( res.betaUpper>0.0 && res.alphaUpper>0.0 );
  return res;
}

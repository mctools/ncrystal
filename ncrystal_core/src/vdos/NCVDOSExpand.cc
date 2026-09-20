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
#include "NCrystal/internal/utils/NCFileUtils.hh"
#include <fstream>
#include <sstream>

namespace NC=NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace VDOS {
    namespace {
      struct ABRangeInfo {
        Rectangle full;
        Rectangle overlap;
      };

      //FIXME: Turn into query instead, and adopt plots into generic utility
      //script?
      void writeFileWithABRanges( const std::vector<ABRangeInfo>& r,
                                  double targetEmax_div_kT )
      {
        std::string fn = "ncrystal_alphabetaranges.txt";
        NCRYSTAL_WARN("Writing Alpha/Beta ranges to "<<fn
                      <<" if it does not already exist");
        if (file_exists(fn))
          return;
        std::ofstream ofs(fn.c_str(), std::ofstream::out);
        ofs << "#ncrystal_alphabetaranges\n";
        ofs << "#targetEmax_div_kT: "<<fmt(targetEmax_div_kT)<<"\n";
        ofs << "#note: all overlap_xxx=-1"
          " if not within phasespace.\n";
        ofs << "#colnames = n,alow,aup,blow,bup,"
          "overlap_alow,overlap_aup,overlap_blow,overlap_bup\n";
        for ( std::size_t i = 0; i < r.size(); ++i ) {
          auto& a = r.at(i).full.xRange();
          auto& b = r.at(i).full.yRange();
          ofs<<i+1
             <<" "<<fmt(a.first)<<" "<<fmt(a.second)
             <<" "<<fmt(b.first)<<" "<<fmt(b.second);
          if ( !r.at(i).overlap.isEmpty() ) {
            auto ao = r.at(i).overlap.xRange();
            auto bo = r.at(i).overlap.yRange();
            ofs <<" "<<fmt(ao.first)<<" "<<fmt(ao.second)
                <<" "<<fmt(bo.first)<<" "<<fmt(bo.second)
                <<"\n";
          } else {
            ofs << " -1 -1 -1 -1\n";
          }
        }
        ofs.close();
      }
    }
  }
}

NC::VDOS::GnExpansion
NC::VDOS::expandVDOSToGnFcts( const VDOSData& vdosdata,
                              VDOSLux vdoslux,
                              Optional<NeutronEnergy> targetEmax_requested )
{
  static bool s_verbose = ncgetenv_bool("DEBUG_PHONON");

  //Hidden unofficial env-vars used for special debugging purposes:
  static const bool dump_vdosabranges = ncgetenv_bool("HACK_DUMP_VDOSABRANGES");

  //Which Emax should we target (i.e. aim to cover the kinematic reachable area
  //for neutrons of that energy):
  constexpr NeutronEnergy lowestEmaxPossible{ 1e-15 };
  auto targetEMaxValid = [lowestEmaxPossible](NeutronEnergy v)
  {
    return ( v.dbl()>0.0
             && std::isfinite(v.dbl())
             && valueInInterval(lowestEmaxPossible.dbl(),1e4,v.dbl()) );
  };
  if ( targetEmax_requested.has_value()
       && !targetEMaxValid(targetEmax_requested.value()) ) {
    NCRYSTAL_THROW2(BadInput,"Invalid or out-of-range targetEmax value"
                    " requested: "<<targetEmax_requested.value());
  }

  NeutronEnergy targetEmax;
  if ( targetEmax_requested.has_value() ) {
    targetEmax = targetEmax_requested.value();
  } else {
    if ( vdoslux.isLegacy() ) {
      //used to depend on vdoslux.
      nc_assert_always( vdoslux.lvl() <= 5u );
      double legacy_lux2emax[6] = { 0.5, 1.0, 3.0, 5.0, 8.0, 12.0 };
      targetEmax = NeutronEnergy{ legacy_lux2emax[vdoslux.lvl()] };
    } else {
      //always exactly 5.0 eV now.
      targetEmax = NeutronEnergy{ 5.0 };
    }
  }

  if ( s_verbose )
    NCRYSTAL_MSG("VDOS expansion initialising with T="<<vdosdata.temperature()
                 <<", vdoslux="<<vdoslux.raw()
                 <<", aiming for Emax="<<targetEmax
                 <<(targetEmax_requested.has_value()?" (as requested)":""));

  nc_assert_always( targetEMaxValid(targetEmax) );

  //Initialise evaluators:
  VDOSEval vdoseval(vdosdata);
  const double kT = vdoseval.kT();
  const double invkT = 1.0/kT;
  const double gamma0 = vdoseval.calcGamma0();
  const double msd = vdoseval.getMSD( gamma0 );
  double targetEmax_div_kT = targetEmax.dbl()*invkT;
  unsigned max_phonon_order = 1;
  auto vdosgn_cfg = ( vdoslux.isLegacy()
                      ? VDOSGn::Cfg::Legacy
                      : VDOSGn::Cfg::Default );

  constexpr double alpha2x_factor = ( 2.0*const_neutron_mass_evc2
                                      /(constant_hbar*constant_hbar) );
  const double alpha2x = alpha2x_factor*kT*msd;
  nc_assert_always( floateq( msd, alpha2x/(alpha2x_factor*kT) ) );

  GnExpansion res{ VDOSGn(vdoseval,vdosgn_cfg), alpha2x,
                   NeutronEnergy{0.0}, Rectangle() };
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

  const NeutronEnergy emax_lowest_allowed
    = targetEmax_requested.value_or(lowestEmaxPossible);

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
      return Rectangle(alphaRange,betaRange);
    };

  while (true) {
    Gn_asym.growMaxOrder(max_phonon_order);
    auto abRange = findAlphaBetaRangeOfOrder(Gn_asym.maxOrder().value());
    if (!findABExtentWithinKB(abRange,targetEmax_div_kT).isEmpty()) {
      //Could consider larger stepsize, but need to carefully check usage in the
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
      NeutronEnergy targetEmax_reduced  = targetEmax;
      do {
        targetEmax_reduced.dbl() *= 0.99;
        if ( targetEmax_reduced < emax_lowest_allowed )
          NCRYSTAL_THROW2(CalcError,"VDOS expansion too slow - can not reach E="
                          <<emax_lowest_allowed<<" after "<<order_limit
                          <<" phonon convolutions (likely causes: either"
                          " the target energy value is too high, vdoslux too"
                          " low, the temperature too high, or the VDOS is"
                          " very unusual).");
      } while (!findABExtentWithinKB(abRange,
                                     targetEmax_reduced.dbl()*invkT).isEmpty());

      if (s_verbose)
        NCRYSTAL_WARN("VDOS expansion could only reach Emax="
                      <<targetEmax_reduced
                      <<" and not the requested Emax="<<targetEmax<<"K");
      targetEmax = targetEmax_reduced;
      targetEmax_div_kT = targetEmax.dbl() * invkT;
      break;
    }
  }
  Gn_asym.growMaxOrder(max_phonon_order);

  //Record the actual Emax:
  res.suggestedEmax = targetEmax;

  //Ok, we now know how many orders we need to reach targetEmax. Next step is to
  //look at the contribution of each order insided the kinematic reach of
  //targetEmax, and use it to determine alpha/beta limits:

  Optional<std::vector<ABRangeInfo>> abRangesForWrite;
  if( dump_vdosabranges ) {
    abRangesForWrite.emplace();
    abRangesForWrite.value().reserve(max_phonon_order);
  }

  for ( unsigned n = 1; n<=max_phonon_order; ++n ) {
    auto abRange = findAlphaBetaRangeOfOrder(n);
    auto abOverlap = findABExtentWithinKB( abRange, targetEmax_div_kT );
    //FIXME: Also store abOverlap for each order n, so the combinedGnFunction
    //knows which parts of each Gn (and alpha) function to consider. Also
    //fillSABFromVDOS could perhaps take advantage.
    res.sabRange = res.sabRange.getUnion( abOverlap );
    if ( abRangesForWrite.has_value() )
      abRangesForWrite.value().push_back( { abRange, abOverlap } );
  }
  if ( abRangesForWrite.has_value() )
    writeFileWithABRanges(abRangesForWrite.value(),targetEmax_div_kT);

  return res;
}

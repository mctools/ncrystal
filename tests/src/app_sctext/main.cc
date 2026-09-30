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

//Validation of SABSCTExtender, the short-collision-time (SCT) based
//alternative to the free-gas SABFGExtender: checks the exact Teff==T
//free-gas limit, cross sections and sampled moments against independent
//2-D numerical integration of the SCT law, and compares both extenders
//against the exact (converged phonon expansion) kernel of H in
//polyethylene at the kernel boundary Emax. Timing is available via
//NCRYSTAL_SCTEXT_TIMING=1 (not part of the reference log).

#include "NCrystal/internal/phys_utils/NCSCTUtils.hh"
#include "NCrystal/internal/sab/NCSABExtender.hh"
#include "NCrystal/internal/sab/NCSABExtended.hh"
#include "NCrystal/internal/sab/NCSABCfg.hh"
#include "NCrystal/internal/sab/NCSABProcessor.hh"
#include "NCrystal/internal/extd_utils/NCSABAnalyser.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/phys_utils/NCKinUtils.hh"
#include "NCrystal/internal/phys_utils/NCFreeGasUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/interfaces/NCRNG.hh"
#include <iostream>
#include <chrono>
#include <cstdio>

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)
#define REQUIREFLTEQ(x,y,tol) nc_assert_always(NC::floateq((x),(y),(tol),0.0))

namespace {

  //Reference implementation of the model laws, in NCrystal alpha units
  //and units of the actual kT, with r=Teff/T (common constant factors
  //dropped throughout, all usage is via ratios):
  //
  //Free gas at Teff: S(alpha,beta) = exp(-(am+beta)^2/(4*am*r))/sqrt(am*r),
  //am=alpha/A (obeys detailed balance at Teff). The SCT law is identical
  //for beta<=0, and additionally suppressed by exp(-beta*(1-1/r)) for
  //beta>0 (moving its detailed balance to the actual T):

  double refShapeFGTeff( double alpha, double beta, double A, double r )
  {
    const double am = alpha / A;
    nc_assert( am > 0.0 );
    return std::exp( -NC::ncsquare(am+beta)/(4.0*am*r) ) / std::sqrt(am*r);
  }

  double upSuppression( double beta, double r )
  {
    return beta <= 0.0 ? 1.0 : std::exp( -beta*(1.0-1.0/r) );
  }

  //2-D integral of weight(alpha,beta)*S over the kinematic domain at
  //c=E/kT, for the SCT (sct=true) or free-gas-at-Teff (sct=false)
  //shape. Inner alpha-integral uses u=sqrt(alpha) to regularise the
  //integrable 1/sqrt(alpha) endpoint singularity at beta~0:
  template<class WFct>
  double integrateShape( double c, double A, double r, bool sct,
                         const WFct& weight )
  {
    auto betaIntegrand = [&]( double beta )
    {
      auto al = NC::getAlphaLimits( c, beta );
      if ( !( al.second > al.first ) )
        return 0.0;
      auto innerIntegrand = [&]( double u )
      {
        const double alpha = u*u;
        if ( !(alpha>0.0) )
          return 0.0;
        return 2.0 * u * refShapeFGTeff(alpha,beta,A,r)
          * weight(alpha,beta);
      };
      double v = NC::integrateRombergFlex( innerIntegrand,
                                           std::sqrt(al.first),
                                           std::sqrt(al.second),
                                           1e-9, 3, 10 );
      return sct ? v * upSuppression(beta,r) : v;
    };
    //NB: the upscatter tail of the SCT shape decays like exp(-beta) in
    //these actual-T units, but the free-gas-at-Teff shape only decays
    //like exp(-beta/r), so its integration limits must scale with r:
    const double s = ( sct ? 1.0 : r );
    const double edges[7] = { -c, -0.5*c, -0.125*c, 0.0,
                              2.0*s, 10.0*s, 80.0*s };
    double sum = 0.0;
    for ( auto i : NC::ncrange(6) )
      sum += NC::integrateRombergFlex( betaIntegrand, edges[i],
                                       edges[i+1], 1e-8, 3, 10 );
    return sum;
  }

  struct SampledStats {
    double mean_beta, var_beta, upfrac, wstat;
    std::size_t n;
  };

  //Sample n times, collecting mean/variance of beta, upscatter
  //fraction, and the "ridge width statistic" W = <(beta+am)^2/(2am)>
  //which equals 1 for a free gas at T and Teff/T in the impulse/SCT
  //limit:
  template<class TSampler>
  SampledStats sampleStats( const TSampler& sampler, NC::RNG& rng,
                            NC::NeutronEnergy ekin, double A,
                            std::size_t n )
  {
    NC::StableSum sum_b, sum_b2, sum_w;
    std::size_t nup(0);
    for ( auto i : NC::ncrange(n) ) {
      (void)i;
      auto ab = sampler( rng, ekin );
      const double am = ab.first / A;
      sum_b.add( ab.second );
      sum_b2.add( NC::ncsquare(ab.second) );
      if ( am > 0.0 )
        sum_w.add( NC::ncsquare( ab.second + am ) / ( 2.0 * am ) );
      if ( ab.second > 0.0 )
        ++nup;
    }
    SampledStats res;
    res.n = n;
    res.mean_beta = sum_b.sum() / n;
    res.var_beta = sum_b2.sum() / n - NC::ncsquare(res.mean_beta);
    res.upfrac = double(nup) / n;
    res.wstat = sum_w.sum() / n;
    return res;
  }

  void test_kernel_closedform()
  {
    //The closed erf-form free-gas beta-kernel versus direct numerical
    //integration of the asymmetric free-gas S over the kinematically
    //allowed alpha range (u=sqrt(alpha) substitution):
    std::printf("test_kernel_closedform:\n");
    unsigned nchecked(0);
    for ( double A : { 1.008/1.00866491588, 26.98/1.00866491588 } ) {
      const double invA = 1.0 / A;
      for ( double c : { 0.5, 5.0, 50.0 } ) {
        for ( double beta : { 0.0, 0.3, 1.0, 3.0, 10.0, 30.0 } ) {
          auto alim = NC::getAlphaLimits( c, beta );
          if ( !( alim.second > alim.first ) )
            continue;
          auto integrand = [beta]( double u )
          {
            const double t = ( u > 0.0
                               ? 0.5*u + beta/(2.0*u)
                               : 0.5*u );
            return NC::kInvSqrtPi * std::exp( -t*t );
          };
          const double ref
            = NC::integrateRombergFlex( integrand,
                                        std::sqrt(alim.first*invA),
                                        std::sqrt(alim.second*invA),
                                        1e-12, 3, 12 );
          const double val
            = NC::SCTXSProvider::evalFGBetaKernel( c, beta, invA );
          REQUIREFLTEQ( val, ref, 1e-7 );
          ++nchecked;
        }
      }
    }
    std::printf("  closed form matches quadrature in %u cases: OK\n",
                nchecked );
  }

  void test_corr_spline()
  {
    //The production upscatter-correction evaluation (splined lookup
    //table + saturated closed form) versus direct quadrature:
    std::printf("test_corr_spline:\n");
    const NC::Temperature T{293.15};
    struct Case { double amu, r; };
    for ( auto cs : { Case{1.008,4.12}, Case{1.008,20.0},
                      Case{1.008,59.3}, Case{207.2,1.05} } ) {
      const NC::AtomMass mass{cs.amu};
      const NC::Temperature Teff{ T.dbl() * cs.r };
      NC::SCTXSProvider p( T, Teff, mass, NC::SigmaBound{1.0} );
      const double A = mass.relativeToNeutronMass();
      const double csat = 35.0 * NC::ncmax( A, 1.0/A );
      double worst = 0.0;
      for ( auto i : NC::ncrange(20) ) {
        const double c = 2e-4 * std::pow( 0.999*csat/2e-4,
                                          i / 19.0 );
        const double v1 = p.evalUpscatterCorr( c );
        const double v2 = p.evalUpscatterCorrIntegral( c );
        if ( v2 > 1e-10 )
          worst = NC::ncmax( worst, NC::ncabs( v1/v2 - 1.0 ) );
      }
      REQUIRE( worst < 1e-3 );
      //Continuity at the saturation point:
      const double rr = ( cs.r - 1.0 ) / cs.r;
      REQUIREFLTEQ( p.evalUpscatterCorrIntegral( csat ), rr, 1e-10 );
      REQUIRE( p.evalUpscatterCorr( 1.001*csat ) == rr );
      std::printf("  A=%.3f r=%.2f : spline vs quadrature worst"
                  " reldiff %.2g: OK\n", A, cs.r, worst );
    }
  }

  void test_r1_limit()
  {
    std::printf("test_r1_limit:\n");
    const NC::Temperature T{293.15};
    const NC::AtomMass mass{1.008};
    NC::SAB::SABFGExtender fg( T, mass, NC::SigmaBound{1.0} );
    NC::SAB::SABSCTExtender sct( T, T, mass, NC::SigmaBound{1.0} );
    for ( double e : { 0.6, 5.0, 50.0 } ) {
      NC::NeutronEnergy ekin{e};
      REQUIRE( fg.crossSection(ekin).dbl()
               == sct.crossSection(ekin).dbl() );
      auto rng1 = NC::createBuiltinRNG( 12345 );
      auto rng2 = NC::createBuiltinRNG( 12345 );
      for ( auto i : NC::ncrange(1000) ) {
        (void)i;
        auto ab1 = fg.sampleAlphaBeta( *rng1, ekin );
        auto ab2 = sct.sampleAlphaBeta( *rng2, ekin );
        REQUIRE( ab1.first == ab2.first && ab1.second == ab2.second );
      }
    }
    std::printf("  Teff==T reproduces free-gas extender exactly: OK\n");
  }

  void test_xs_vs_bruteforce()
  {
    std::printf("test_xs_vs_bruteforce:\n");
    const NC::Temperature T{293.15};
    const double kT = T.kT();
    struct Case { double amu, r; };
    for ( auto cs : { Case{1.008,1.5}, Case{1.008,4.12},
                      Case{1.008,20.0}, Case{26.98,1.5} } ) {
      const NC::AtomMass mass{cs.amu};
      const double A = mass.relativeToNeutronMass();
      const NC::Temperature Teff{ T.dbl() * cs.r };
      NC::SAB::SABSCTExtender sct( T, Teff, mass, NC::SigmaBound{1.0} );
      NC::FreeGasXSProvider fgteff( Teff, mass, NC::SigmaBound{1.0} );
      for ( double e : { 0.7, 3.0, 5.0, 20.0 } ) {
        NC::NeutronEnergy ekin{e};
        const double c = e / kT;
        auto w1 = []( double, double ) { return 1.0; };
        const double ratio_ref
          = ( integrateShape( c, A, cs.r, true, w1 )
              / integrateShape( c, A, cs.r, false, w1 ) );
        const double ratio_impl = ( sct.crossSection(ekin).dbl()
                                    / fgteff.crossSection(ekin).dbl() );
        std::printf("  A=%.3f r=%5.2f E=%5.2feV : sigma_SCT/sigma_FG(Teff)"
                    " = %.6f (ref %.6f)\n", A, cs.r, e,
                    ratio_impl, ratio_ref );
        REQUIREFLTEQ( ratio_impl, ratio_ref, 1e-4 );
        REQUIRE( ratio_impl <= 1.0 && ratio_impl > 0.5 );
      }
    }
  }

  void test_sampling_vs_bruteforce()
  {
    std::printf("test_sampling_vs_bruteforce:\n");
    const NC::Temperature T{293.15};
    const double kT = T.kT();
    const NC::AtomMass mass{1.008};
    const double A = mass.relativeToNeutronMass();
    const double r = 4.12;
    const NC::Temperature Teff{ T.dbl() * r };
    NC::SAB::SABSCTExtender sct( T, Teff, mass, NC::SigmaBound{1.0} );
    auto rng = NC::createBuiltinRNG( 999 );
    for ( double e : { 0.7, 6.0 } ) {
      NC::NeutronEnergy ekin{e};
      const double c = e / kT;
      auto w1 = []( double, double ) { return 1.0; };
      auto wb = []( double, double beta ) { return beta; };
      auto wb2 = []( double, double beta ) { return NC::ncsquare(beta); };
      auto wup = []( double, double beta )
      { return beta > 0.0 ? 1.0 : 0.0; };
      const double norm = integrateShape( c, A, r, true, w1 );
      const double ref_mb = integrateShape( c, A, r, true, wb ) / norm;
      const double ref_vb = ( integrateShape( c, A, r, true, wb2 ) / norm
                              - NC::ncsquare(ref_mb) );
      const double ref_up = integrateShape( c, A, r, true, wup ) / norm;
      const std::size_t n = 200000;
      auto stats = sampleStats( [&sct]( NC::RNG& r_, NC::NeutronEnergy e_ )
      { return sct.sampleAlphaBeta(r_,e_); }, *rng, ekin, A, n );
      std::printf("  E=%.1feV : <beta> ref %.4g, Var(beta) ref %.4g,"
                  " upfrac ref %.4g\n", e, ref_mb, ref_vb, ref_up );
      //Statistical tolerances (samples are correlated with nothing, so
      //plain sqrt(n) errors, taken with generous 6 sigma margins):
      const double sd_b = std::sqrt( ref_vb / n );
      REQUIRE( NC::ncabs( stats.mean_beta - ref_mb ) < 6.0 * sd_b );
      REQUIRE( NC::ncabs( stats.var_beta - ref_vb ) < 0.02 * ref_vb );
      const double sd_up = std::sqrt( ref_up * (1.0-ref_up) / n );
      REQUIRE( NC::ncabs( stats.upfrac - ref_up ) < 6.0 * sd_up + 1e-6 );
      std::printf("    sampled values consistent (n=%u): OK\n",
                  unsigned(n) );
    }
  }

  void test_sampling_cdf()
  {
    //Distribution-level validation of SCT sampling, beyond the moment
    //checks above: node-based Kolmogorov-Smirnov comparison of sampled
    //beta values against a reference CDF built by direct numerical
    //integration of the model definition (a statistic bounded by the
    //true KS D, so the KS critical value is valid). This is the
    //independent check of the sampler's actual distribution; note that
    //std-vs-ref sampler comparisons in the sabsample tests share this
    //very sampling code for the extension region, so cannot provide it.
    std::printf("test_sampling_cdf:\n");
    const NC::Temperature T{293.15};
    const double kT = T.kT();
    const NC::AtomMass mass{1.008};
    const double A = mass.relativeToNeutronMass();
    const double r = 4.12;
    const NC::Temperature Teff{ T.dbl() * r };
    NC::SAB::SABSCTExtender sct( T, Teff, mass, NC::SigmaBound{1.0} );
    auto rng = NC::createBuiltinRNG( 4242 );
    for ( double e : { 0.7, 6.0 } ) {
      NC::NeutronEnergy ekin{e};
      const double c = e / kT;
      //Beta-marginal of the SCT law (alpha integrated out), same
      //independent integration approach as integrateShape:
      auto marginal = [c,A,r]( double beta )
      {
        auto al = NC::getAlphaLimits( c, beta );
        if ( !( al.second > al.first ) )
          return 0.0;
        auto innerIntegrand = [beta,A,r]( double u )
        {
          const double alpha = u*u;
          return alpha > 0.0
            ? 2.0 * u * refShapeFGTeff(alpha,beta,A,r) : 0.0;
        };
        const double v
          = NC::integrateRombergFlex( innerIntegrand,
                                      std::sqrt(al.first),
                                      std::sqrt(al.second),
                                      1e-9, 3, 10 );
        return v * upSuppression(beta,r);
      };
      //Node grid mirroring integrateShape's segmentation (upscatter
      //tail beyond beta=30 carries ~exp(-30) mass, negligible):
      const double edges[7] = { -c, -0.5*c, -0.125*c, 0.0,
                                2.0, 10.0, 30.0 };
      constexpr unsigned npseg = 48;
      NC::VectD nodes, cdf;
      nodes.push_back( edges[0] );
      cdf.push_back( 0.0 );
      NC::StableSum cum;
      for ( auto iseg : NC::ncrange(6) ) {
        const double a0 = edges[iseg], a1 = edges[iseg+1];
        for ( auto j : NC::ncrange( 1u, npseg + 1 ) ) {
          const double b1 = a0 + (a1-a0) * ( double(j) / npseg );
          cum.add( NC::integrateRombergFlex( marginal, nodes.back(),
                                             b1, 1e-8, 3, 10 ) );
          nodes.push_back( b1 );
          cdf.push_back( cum.sum() );
        }
      }
      const double tot = cum.sum();
      REQUIRE( tot > 0.0 );
      //Sample and compute the node-based KS statistic:
      constexpr std::size_t n = 200000;
      NC::VectD betas;
      betas.reserve( n );
      for ( auto i : NC::ncrange(n) ) {
        (void)i;
        betas.push_back( sct.sampleAlphaBeta( *rng, ekin ).second );
      }
      std::sort( betas.begin(), betas.end() );
      double ks_d = 0.0;
      for ( auto i : NC::ncrange( nodes.size() ) ) {
        const auto it = std::upper_bound( betas.begin(), betas.end(),
                                          NC::vectAt( nodes, i ) );
        const double f_emp
          = double( std::distance( betas.begin(), it ) ) / n;
        ks_d = NC::ncmax( ks_d,
                          NC::ncabs( f_emp - NC::vectAt( cdf, i ) / tot ) );
      }
      //KS 1%-level critical value is 1.63/sqrt(n); generous margin:
      const double bound = 2.5 / std::sqrt( double(n) );
      std::printf("  E=%.1feV : node-KS D below %.4f (n=%u,"
                  " %u nodes): %s\n", e, bound, unsigned(n),
                  unsigned(nodes.size()), ks_d < bound ? "OK" : "FAIL" );
      REQUIRE( ks_d < bound );
    }
  }

  void test_kernel_boundary()
  {
    std::printf("test_kernel_boundary (H in polyethylene, 293.15K):\n");
    auto info = NC::FactImpl::createInfo(
      NC::MatCfg("Polyethylene_CH2.ncmat") );
    const NC::DI_VDOS* di_h = nullptr;
    for ( auto& di : info->getDynamicInfoList() ) {
      auto p = dynamic_cast<const NC::DI_VDOS*>( di.get() );
      if ( p && p->atomData().isElement()
           && p->atomData().Z() == 1 )
        di_h = p;
    }
    REQUIRE( di_h != nullptr );
    auto sab = NC::extractSABDataFromDynInfo( di_h );
    NC::VDOSEval ve( di_h->vdosData() );
    const NC::Temperature teff{ ve.calcEffectiveTemperature() };
    const NC::Temperature T = sab->temperature();
    const double r = teff.dbl() / T.dbl();
    std::printf("  Teff = %.1fK (Teff/T = %.3f)\n", teff.dbl(), r );
    auto cfg = NC::SABCfg::createConfig(
      NC::SABCfg::sablux_default_luxury );
    auto processor = NC::makeSO<NC::SABUtils::SABProcessor>(
      cfg, sab, nullptr );
    const auto emaxinfo = processor->getEMaxInfo();
    const NC::NeutronEnergy emax = emaxinfo.ekin;
    const double A = sab->elementMassAMU().relativeToNeutronMass();
    NC::SAB::SABFGExtender fg( T, sab->elementMassAMU(),
                               NC::SigmaBound{1.0} );
    NC::SAB::SABSCTExtender sct( T, teff, sab->elementMassAMU(),
                                 NC::SigmaBound{1.0} );
    const double xs_table = emaxinfo.crossSectionUnitSigmaBound.dbl();
    const double xs_fg = fg.crossSection(emax).dbl();
    const double xs_sct = sct.crossSection(emax).dbl();
    std::printf("  At Emax=%.4geV: xs(table)=%.4gbarn,"
                " xs(freegas)=%.4gbarn, xs(sct)=%.4gbarn\n",
                emax.dbl(), xs_table, xs_fg, xs_sct );
    std::printf("  Rel. mismatch to table: freegas %.2f%%, sct %.2f%%\n",
                1e2*NC::ncabs(xs_fg/xs_table-1.0),
                1e2*NC::ncabs(xs_sct/xs_table-1.0) );
    //The ridge-width statistic is the sharp discriminator (equals 1.0
    //for a free gas at T, Teff/T in the SCT/impulse limit, and the
    //converged kernel should sit close to the latter):
    const std::size_t n = 100000;
    auto rng = NC::createBuiltinRNG( 777 );
    auto st_tab = sampleStats(
      [&processor]( NC::RNG& r_, NC::NeutronEnergy e_ )
      { auto ab = processor->sampleScatterAlphaBeta(r_,e_);
        return NC::PairDD{ab.alpha,ab.beta}; }, *rng, emax, A, n );
    auto st_fg = sampleStats(
      [&fg]( NC::RNG& r_, NC::NeutronEnergy e_ )
      { return fg.sampleAlphaBeta(r_,e_); }, *rng, emax, A, n );
    auto st_sct = sampleStats(
      [&sct]( NC::RNG& r_, NC::NeutronEnergy e_ )
      { return sct.sampleAlphaBeta(r_,e_); }, *rng, emax, A, n );
    std::printf("  Ridge width statistic W at Emax:"
                " table %.3g, freegas %.3g, sct %.3g (Teff/T=%.3g)\n",
                st_tab.wstat, st_fg.wstat, st_sct.wstat, r );
    REQUIRE( NC::ncabs( st_sct.wstat - st_tab.wstat )
             < NC::ncabs( st_fg.wstat - st_tab.wstat ) );
    //NB: at the level of the *total* cross section all three agree to
    //better than 1% here (the total is a very forgiving integral, and
    //the free-gas model at T is not actually improved upon by SCT for
    //it) -- the large SCT improvement seen in the W statistic above
    //concerns the secondary energy/angle distributions:
    REQUIRE( NC::ncabs( xs_sct / xs_table - 1.0 ) < 0.02 );
    REQUIRE( NC::ncabs( xs_fg / xs_table - 1.0 ) < 0.02 );

    if ( NC::ncgetenv_bool("SCTEXT_TIMING") ) {
      for ( double amu : { 1.008, 207.2 } ) {
        auto t0 = std::chrono::steady_clock::now();
        constexpr unsigned nctor = 20;
        NC::StableSum dump;
        for ( auto i : NC::ncrange(nctor) ) {
          (void)i;
          NC::SCTXSProvider p( T, NC::Temperature{ T.dbl()*4.12 },
                               NC::AtomMass{amu}, NC::SigmaBound{1.0} );
          dump.add( p.crossSection(NC::NeutronEnergy{1.0}).dbl() );
        }
        auto dt = std::chrono::duration<double,std::milli>(
          std::chrono::steady_clock::now() - t0 ).count();
        std::printf("  TIMING ctor A=%.0f: %.2f ms/ctor (checksum %g)\n",
                    amu, dt/nctor, dump.sum() );
      }
      auto rng2 = NC::createBuiltinRNG( 5 );
      const std::size_t nt = 100000;
      NC::NeutronEnergy et{ emax.dbl() * 1.2 };
      for ( int mode = 0; mode < 4; ++mode ) {
        auto t0 = std::chrono::steady_clock::now();
        NC::StableSum dump;
        for ( auto i : NC::ncrange(nt) ) {
          (void)i;
          if ( mode == 0 )
            dump.add( fg.crossSection(et).dbl() );
          else if ( mode == 1 )
            dump.add( sct.crossSection(et).dbl() );
          else if ( mode == 2 )
            dump.add( fg.sampleAlphaBeta(*rng2,et).second );
          else
            dump.add( sct.sampleAlphaBeta(*rng2,et).second );
        }
        auto dt = std::chrono::duration<double,std::nano>(
          std::chrono::steady_clock::now() - t0 ).count();
        const char* names[4] = { "fg_xs", "sct_xs",
                                 "fg_sample", "sct_sample" };
        std::printf("  TIMING %s: %.0f ns/call (checksum %g)\n",
                    names[mode], dt / nt, dump.sum() );
      }
    }
  }

  void test_cryogenic()
  {
    std::printf("test_cryogenic (H in polyethylene, 20K):\n");
    auto info = NC::FactImpl::createInfo(
      NC::MatCfg("Polyethylene_CH2.ncmat;temp=20K") );
    const NC::DI_VDOS* di_h = nullptr;
    for ( auto& di : info->getDynamicInfoList() ) {
      auto p = dynamic_cast<const NC::DI_VDOS*>( di.get() );
      if ( p && p->atomData().isElement() && p->atomData().Z() == 1 )
        di_h = p;
    }
    REQUIRE( di_h != nullptr );
    NC::VDOSEval ve( di_h->vdosData() );
    const NC::Temperature T{20.0};
    const NC::Temperature teff{ ve.calcEffectiveTemperature() };
    const double r = teff.dbl() / T.dbl();
    const NC::AtomMass mass = di_h->atomData().averageMassAMU();
    const double A = mass.relativeToNeutronMass();
    NC::SAB::SABFGExtender fg( T, mass, NC::SigmaBound{1.0} );
    NC::SAB::SABSCTExtender sct( T, teff, mass, NC::SigmaBound{1.0} );
    const NC::NeutronEnergy ekin{6.0};
    auto rng = NC::createBuiltinRNG( 4242 );
    const std::size_t n = 100000;
    auto st_fg = sampleStats(
      [&fg]( NC::RNG& r_, NC::NeutronEnergy e_ )
      { return fg.sampleAlphaBeta(r_,e_); }, *rng, ekin, A, n );
    auto st_sct = sampleStats(
      [&sct]( NC::RNG& r_, NC::NeutronEnergy e_ )
      { return sct.sampleAlphaBeta(r_,e_); }, *rng, ekin, A, n );
    std::printf("  Teff/T = %.1f ; sampled W at 6eV:"
                " freegas %.3g, sct %.3g\n", r,
                st_fg.wstat, st_sct.wstat );
    REQUIREFLTEQ( st_fg.wstat, 1.0, 0.1 );
    REQUIREFLTEQ( st_sct.wstat, r, 0.1 );
    REQUIRE( sct.crossSection(ekin).dbl() > 0.0 );
  }

}

namespace {
  void test_estimateTeffMSD()
  {
    //Turn-key estimates on a VDOS-expanded kernel must agree with the
    //VDOSEval-derived truth, and repeated (cached) calls must be
    //identical:
    std::printf("test_estimateTeffMSD:\n");
    auto info = NC::FactImpl::createInfo(
      NC::MatCfg("Polyethylene_CH2.ncmat") );
    const NC::DI_VDOS* di_h = nullptr;
    for ( auto& di : info->getDynamicInfoList() ) {
      auto p = dynamic_cast<const NC::DI_VDOS*>( di.get() );
      if ( p && p->atomData().isElement()
           && p->atomData().Z() == 1 )
        di_h = p;
    }
    REQUIRE( di_h != nullptr );
    auto sab = NC::extractSABDataFromDynInfo( di_h );
    NC::VDOSEval ve( di_h->vdosData() );
    const double teff_true = ve.calcEffectiveTemperature();
    const double msd_true = ve.getMSD( ve.calcGamma0() );
    auto e = NC::SABAnalyser::estimateTeffMSD( *sab );
    REQUIRE( e.effectiveTemperature.has_value() );
    REQUIRE( e.msd.has_value() );
    const double teff_est = e.effectiveTemperature.value().dbl();
    REQUIRE( NC::ncabs( teff_est/teff_true - 1.0 ) < 0.02 );
    REQUIRE( NC::ncabs( e.msd.value()/msd_true - 1.0 ) < 0.03 );
    auto e2 = NC::SABAnalyser::estimateTeffMSD( *sab );
    REQUIRE( e2.effectiveTemperature.has_value() && e2.msd.has_value() );
    REQUIRE( e2.effectiveTemperature.value().dbl() == teff_est );
    REQUIRE( e2.msd.value() == e.msd.value() );
    std::printf("  estimates match VDOS truth, cache consistent: OK\n");
  }
}

int main()
{
  test_kernel_closedform();
  test_corr_spline();
  test_r1_limit();
  test_xs_vs_bruteforce();
  test_sampling_vs_bruteforce();
  test_sampling_cdf();
  test_kernel_boundary();
  test_cryogenic();
  test_estimateTeffMSD();
  std::printf("All tests passed.\n");
  return 0;
}

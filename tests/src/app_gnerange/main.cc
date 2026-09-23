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

////////////////////////////////////////////////////////////////////////////////
// Testbed for NC::VDOS::estimateGnErange -- a standalone utility, extracted
// from (but NOT YET wired into) VDOSGn::eRange, for estimating the energy
// interval outside of which a Gn spectrum is everywhere below
// relcontriblvl*max(spec). See docs/claude_session_vdos_fma_reprod.md for the
// background: production VDOSGn::eRange snaps the crossing to whichever grid
// point first clears the threshold. This (a) discretises a continuous
// root-finding problem to bin-width resolution, amplifying a tiny (last-ULP,
// cross-platform-varying) input difference into a whole-bin-sized output
// difference whenever the true crossing sits close to a bin edge, and (b) is
// "one-sided": the single point that happens to first clear the threshold is
// trusted as ground truth, with no check on whether that very point's own
// value might be the one affected by such noise.
//
// This app builds synthetic curves (a clean Gaussian-shaped bump, the same
// bump with noise injected right at the naive crossing points, and a hard
// step) plus real Gn spectra generated LIVE from the actual production
// VDOSGn/VDOSEval machinery -- not hand-transcribed, so a broad sweep of
// "nuisance" cases can be tried cheaply and without risk of transcription
// error. The live sweep covers: three Debye models spanning very different
// characteristic scales (narrow/steep, Al-like, wide), each swept across
// phonon orders from 1 (where a Debye model's hard cutoff at the Debye
// energy is barely smoothed by convolution) up to 100 (where the
// central-limit-theorem-driven Gaussian shape dominates) and across
// relcontriblvl spanning the full vdoslux 0-6 range (1e-3 down to 1e-13);
// plus a real multi-element file-based material (Li2O, both elements),
// loaded exactly as NCrystal itself would load it, swept the same way.
// This deliberately includes the two specific cases (a Li2O/Li spectrum,
// and a Debye-model spectrum) that showed up as cross-platform-sensitive on
// actual CI runs earlier this session (see
// docs/claude_session_vdos_fma_reprod.md), now generalised into a much
// broader sweep around them. Compares four crossing estimators against
// each other (and, for the synthetic cases, against the true analytic
// answer):
//
//  * snapToGrid:       today's production VDOSGn::eRange behaviour.
//  * linear2pt:        linear interpolation between the immediate bracket.
//  * logLinear2pt:      log-space linear interpolation between the bracket.
//  * estimateGnErange: the new windowed quadratic-log-fit utility
//                      (NCVDOSUtils.cc): high-order Gn spectra approach a
//                      Gaussian shape via the central limit theorem, i.e.
//                      ln(spec) is close to quadratic (not linear) in the
//                      tails, so a windowed quadratic fit was tried after an
//                      initial windowed LINEAR fit turned out to be *worse*
//                      than the plain 2-point logLinear2pt on the synthetic
//                      Gaussian case (curvature bias outweighing the
//                      noise-averaging benefit of widening the window) --
//                      the quadratic fit gets both: essentially exact on the
//                      noise-free Gaussian case (errors at the level of
//                      floating-point round-off), and far less sensitive
//                      than logLinear2pt to noise injected at the bracket
//                      points (a 1e-3 relative perturbation there moves the
//                      quadratic fit's answer by only ~1e-6 relative,
//                      instead of directly showing up at the ~1e-3 level).
//
// Deliberately verbose/exploratory rather than a narrow pass/fail check --
// this is meant to be read, and rerun with variations, while the algorithm in
// estimateGnErange is refined.
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/vdos/NCVDOSUtils.hh"
#include "NCrystal/internal/vdos/NCVDOSGn.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <iostream>
#include <sstream>

namespace NC = NCrystal;
namespace NCV = NCrystal::VDOS;

namespace {

  //////////////////////////////////////////////////////////
  // Reference implementations of the other three methods //
  //////////////////////////////////////////////////////////

  NC::PairDD snapToGrid( NC::Span<const double> egrid,
                         NC::Span<const double> spec, double threshold )
  {
    NC::PairDD r( egrid.front(), egrid.back() );
    for ( std::size_t i = 0; i < spec.size(); ++i ) {
      if ( spec[i] >= threshold ) {
        r.first = egrid[i];
        break;
      }
    }
    for ( std::size_t i = spec.size(); i > 0; --i ) {
      if ( spec[i-1] >= threshold ) {
        r.second = NC::ncmin( r.second, egrid[i-1] );
        break;
      }
    }
    return r;
  }

  double bracketInterp( NC::Span<const double> egrid, NC::Span<const double> spec,
                        std::size_t idxLo, std::size_t idxHi, double threshold,
                        bool logspace )
  {
    const double v0 = spec[idxLo];
    const double v1 = spec[idxHi];
    double t;
    if ( logspace && v0 > 0.0 && v1 > 0.0 )
      t = ( std::log(threshold)-std::log(v0) ) / ( std::log(v1)-std::log(v0) );
    else
      t = (threshold-v0)/(v1-v0);
    return NC::nclerp( egrid[idxLo], egrid[idxHi], NC::ncclamp(t,0.0,1.0) );
  }

  NC::PairDD twoPoint( NC::Span<const double> egrid, NC::Span<const double> spec,
                       double threshold, bool logspace )
  {
    NC::PairDD r( egrid.front(), egrid.back() );
    for ( std::size_t i = 0; i < spec.size(); ++i ) {
      if ( spec[i] >= threshold ) {
        r.first = ( i==0 ? egrid.front()
                   : bracketInterp(egrid,spec,i-1,i,threshold,logspace) );
        break;
      }
    }
    for ( std::size_t i = spec.size(); i > 0; --i ) {
      if ( spec[i-1] >= threshold ) {
        const double x = ( i==spec.size() ? egrid[i-1]
                          : bracketInterp(egrid,spec,i-1,i,threshold,logspace) );
        r.second = NC::ncmin( r.second, x );
        break;
      }
    }
    return r;
  }

  //////////////////////
  // Reporting helper //
  //////////////////////

  void report( const char* label,
              NC::Span<const double> egrid, NC::Span<const double> spec,
              double relcontriblvl,
              const NC::PairDD* truth = nullptr )
  {
    const double spec_max = *std::max_element(spec.begin(),spec.end());
    const double threshold = relcontriblvl*spec_max;
    std::cout << "=== " << label << " (npts=" << spec.size()
              << ", relcontriblvl=" << NC::fmtg(relcontriblvl) << ") ===" << std::endl;
    auto pr = [&]( const char* name, NC::PairDD r )
    {
      std::cout << "  " << name << ": [" << NC::fmt(r.first,"%.10g")
                << ", " << NC::fmt(r.second,"%.10g") << "]";
      if ( truth ) {
        const double elo = r.first - truth->first;
        const double ehi = r.second - truth->second;
        std::cout << "  (err_lo=" << NC::fmt(elo,"%.3g")
                  << ", err_hi=" << NC::fmt(ehi,"%.3g") << ")";
      }
      std::cout << std::endl;
    };
    pr( "snapToGrid      ", snapToGrid(egrid,spec,threshold) );
    pr( "linear2pt       ", twoPoint(egrid,spec,threshold,false) );
    pr( "logLinear2pt    ", twoPoint(egrid,spec,threshold,true) );
    pr( "estimateGnErange", NC::VDOS::estimateGnErange(egrid,spec,relcontriblvl) );
    if ( truth )
      std::cout << "  truth           : [" << NC::fmt(truth->first,"%.10g")
                << ", " << NC::fmt(truth->second,"%.10g") << "]" << std::endl;
    std::cout << std::endl;
  }

  ///////////////////////////////////////////////
  // Synthetic case 1+2: clean/noisy Gaussian.  //
  ///////////////////////////////////////////////

  //Gaussian bump spec(x) = exp(-k*(x-mu)^2), sampled on an equidistant grid
  //of npts points spanning [mu-halfwidth,mu+halfwidth]. Truth (continuous
  //crossing of relcontriblvl*spec_max=relcontriblvl, since spec_max=1 at
  //x=mu) is exactly mu +- sqrt(ln(1/relcontriblvl)/k), modelling the
  //Gaussian-like shape real high-order Gn spectra approach (central limit
  //theorem acting on repeated self-convolution):
  void gaussianCase( const char* label, double k, double mu, double halfwidth,
                     std::size_t npts, double relcontriblvl,
                     double noiseAmplRel = 0.0 )
  {
    NC::VectD egrid(npts), spec(npts);
    const double binwidth = 2.0*halfwidth/(npts-1);
    for ( auto i : NC::ncrange(npts) ) {
      const double x = mu - halfwidth + i*binwidth;
      egrid[i] = x;
      spec[i] = std::exp( -k*(x-mu)*(x-mu) );
    }
    if ( noiseAmplRel != 0.0 ) {
      //Perturb specifically the points nearest the naive crossing on both
      //sides, mimicking a last-ULP-level fluctuation landing exactly on the
      //single point the old snap-to-grid/2pt methods would trust as their
      //sole anchor -- this is the "one-sidedness" concern:
      const double threshold = relcontriblvl;//spec_max=1
      for ( auto i : NC::ncrange(std::size_t(1),npts) ) {
        const bool crossesUp = spec[i-1] < threshold && spec[i] >= threshold;
        const bool crossesDown = spec[i-1] >= threshold && spec[i] < threshold;
        if ( crossesUp || crossesDown ) {
          spec[i-1] *= ( 1.0 + noiseAmplRel );
          spec[i] *= ( 1.0 - noiseAmplRel );
        }
      }
    }
    //Truth must be relative to the actual (discretely sampled) spec_max, not
    //the continuous peak value of 1.0 -- these differ whenever the grid does
    //not happen to include a point exactly at mu, which would otherwise
    //contaminate the error comparison below with a discretisation artifact
    //unrelated to which crossing-estimator is more accurate:
    const double spec_max = *std::max_element(spec.begin(),spec.end());
    const double d = std::sqrt( std::log( 1.0/(relcontriblvl*spec_max) )/k );
    NC::PairDD truth( mu-d, mu+d );
    report( label, egrid, spec, relcontriblvl, &truth );
  }

  ///////////////////////////////////
  // Synthetic case 3: a hard step //
  ///////////////////////////////////

  //spec(x) = 1 for x in [mu-w,mu+w], 0 outside: no interpolation scheme can
  //improve on the true crossing here (it coincides with a grid point by
  //construction), but this checks that none of the methods misbehave
  //(overshoot outside the window, nan, etc) when confronted with a genuine
  //discontinuity rather than a smooth tail:
  void stepCase( double mu, double w, std::size_t npts, double relcontriblvl )
  {
    NC::VectD egrid(npts), spec(npts);
    const double halfwidth = 2.0*w;
    const double binwidth = 2.0*halfwidth/(npts-1);
    for ( auto i : NC::ncrange(npts) ) {
      const double x = mu - halfwidth + i*binwidth;
      egrid[i] = x;
      spec[i] = ( NC::ncabs(x-mu) <= w ) ? 1.0 : 0.0;
    }
    report( "hard step (discontinuous, robustness check)",
            egrid, spec, relcontriblvl );
  }

  ///////////////////////////////////////////////////////////////////
  // Real data, generated LIVE from the actual production VDOSGn/    //
  // VDOSEval machinery (not hand-transcribed), so a broad sweep of  //
  // "nuisance" cases can be tried cheaply and without any risk of   //
  // transcription error: many phonon orders (from n=1, where the    //
  // underlying VDOS's own shape -- including a Debye model's hard   //
  // cutoff -- is barely smoothed by convolution, up to high orders  //
  // where the central-limit-theorem-driven Gaussian shape dominates),//
  // many relcontriblvl values (spanning the full vdoslux 0-6 range), //
  // several very different Debye temperatures (narrow vs wide        //
  // characteristic spectra), and a real multi-element file-based     //
  // material (Li2O, both elements) loaded exactly as NCrystal itself //
  // would load it.                                                   //
  ///////////////////////////////////////////////////////////////////

  void sweepSpectrum( const char* label, const NC::VDOSData& vdosdata,
                      NC::Span<const unsigned> orders,
                      NC::Span<const double> relcontriblvls )
  {
    NC::VDOSEval ve(vdosdata);
    NCV::VDOSGn gn( ve, NCV::VDOSGn::Cfg::Default );
    unsigned maxorder = 0;
    for ( auto n : orders ) maxorder = std::max(maxorder,n);
    gn.growMaxOrder( maxorder );
    for ( auto n : orders ) {
      const NCV::VDOSGn::Order order(n);
      const auto& spec = gn.getRawSpectrum( order );
      const double lower = gn.eRange( order ).first;//exact grid edge
      const double binwidth = gn.binWidth( order );
      NC::VectD egrid( spec.size() );
      for ( auto i : NC::ncrange(spec.size()) )
        egrid[i] = NCV::equidistantGridPoint( lower, binwidth, i );
      for ( auto relcontriblvl : relcontriblvls ) {
        std::ostringstream oss;
        oss << label << ", order=" << n;
        report( oss.str().c_str(), egrid, spec, relcontriblvl );
      }
    }
  }

  NC::VDOSData makeDebyeVDOS( double debyeTempK, double tempK )
  {
    return NC::createVDOSDebye( NC::DebyeTemperature{debyeTempK},
                                NC::Temperature{tempK},
                                NC::SigmaBound{1.5},
                                NC::AtomMass{26.98} );
  }

  void debyeSweep()
  {
    //vdoslux 0..6 (non-legacy) span relcontriblvl 1e-6..1e-13; legacy spans
    //1e-3..1e-13. Orders from 1 (barely-convolved: for a Debye model this is
    //where the input's own hard cutoff at the Debye energy is most directly
    //visible) up through orders high enough that the central-limit-theorem
    //Gaussian shape dominates:
    const unsigned ordersArr[] = { 1, 2, 3, 5, 10, 30, 100 };
    const double lvlsArr[] = { 1e-3, 1e-6, 1e-9, 1e-12, 1e-13 };
    //Three very different characteristic scales (Debye temp vs material
    //temp ratio controls how "steep"/narrow the resulting spectrum is):
    sweepSpectrum( "Debye T=80K,temp=293.15K (narrow/steep)",
                   makeDebyeVDOS(80.0,293.15), ordersArr, lvlsArr );
    sweepSpectrum( "Debye T=410K,temp=293.15K (Al-like)",
                   makeDebyeVDOS(410.0,293.15), ordersArr, lvlsArr );
    sweepSpectrum( "Debye T=2000K,temp=293.15K (wide)",
                   makeDebyeVDOS(2000.0,293.15), ordersArr, lvlsArr );
  }

  void li2oSweep()
  {
    NC::MatCfg cfg( "Li2O_sg225_LithiumOxide.ncmat;vdoslux=3;temp=293.15K" );
    auto info = NC::FactImpl::createInfo( cfg );
    const unsigned ordersArr[] = { 1, 2, 3, 5, 20, 60 };
    const double lvlsArr[] = { 1e-3, 1e-6, 1e-9, 1e-12, 1e-13 };
    for ( auto& di : info->getDynamicInfoList() ) {
      auto di_vdos = dynamic_cast<const NC::DI_VDOS*>( di.get() );
      if ( !di_vdos )
        continue;
      std::ostringstream oss;
      oss << "Li2O (real material), element=" << di->atomData().elementName();
      sweepSpectrum( oss.str().c_str(), di_vdos->vdosData(), ordersArr, lvlsArr );
    }
  }
}

int main() {
  std::cout << "--- Synthetic: clean Gaussian bump, varying resolution ---" << std::endl;
  gaussianCase( "gaussian, coarse (24 pts)", 5.0, 0.0, 3.0, 24, 1e-6 );
  gaussianCase( "gaussian, medium (60 pts)", 5.0, 0.0, 3.0, 60, 1e-6 );
  gaussianCase( "gaussian, fine (400 pts)",  5.0, 0.0, 3.0, 400, 1e-6 );

  std::cout << "--- Synthetic: same coarse Gaussian, with noise injected"
               " right at the naive crossing points (the \"one-sidedness\""
               " concern) ---" << std::endl;
  gaussianCase( "gaussian, coarse, +-1e-6 rel noise at crossing", 5.0, 0.0, 3.0,
               24, 1e-6, 1e-6 );
  gaussianCase( "gaussian, coarse, +-1e-3 rel noise at crossing", 5.0, 0.0, 3.0,
               24, 1e-6, 1e-3 );

  std::cout << "--- Synthetic: hard step (discontinuity robustness) ---" << std::endl;
  stepCase( 0.0, 1.0, 40, 1e-6 );

  std::cout << "--- Real data: live Debye-model sweep"
               " (many orders/temperatures/relcontriblvl) ---" << std::endl;
  debyeSweep();

  std::cout << "--- Real data: live Li2O (Li,O) sweep"
               " (many orders/relcontriblvl) ---" << std::endl;
  li2oSweep();

  return 0;
}

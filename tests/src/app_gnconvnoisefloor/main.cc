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
// Testbed for NC::VDOS::estimateFFTConvolutionNoiseFloor -- a standalone
// utility, NOT YET wired into VDOSGn::Impl::produceNewOrderByConvolutionImpl,
// for estimating how large FastConvolve's own FFT round-off can get. See
// docs/claude_session_vdos_fma_reprod.md for the background: production
// truncates a self-convolved phonon order's spectrum at the first/last point
// exceeding truncationThreshold*spec_max (1e-13), with no margin. For high
// enough orders (large FFT size N, deep tails), that bare 1e-13*peak cutoff
// can itself sit at or below the convolution's own inherent round-off floor,
// so whichever point happens to land a hair above it is decided by chance
// round-off rather than the true (still physically meaningful, smoothly
// decaying) signal -- confirmed cross-platform on real vdoslux=2004 data
// (orders 2 and 3 of Li_from_Li2O.ncmat;temp=10, relative divergence
// ~3e-4-1e-3 right at the truncation edge, vs the ordinary ~1e-16 floor
// everywhere else).
//
// This app calibrates and validates estimateFFTConvolutionNoiseFloor's
// peak*eps*sqrt(n) scaling and safetyFactor against a ground truth: an O(n^2)
// direct (Neumaier-summed) reference convolution, which is accurate to
// essentially 1ULP regardless of n, so the difference between it and
// FastConvolve's FFT result is (unlike a plain-vs-mfma build comparison) a
// direct measurement of the FFT algorithm's true numerical error -- no second
// build required. Checked on:
//
//  * Synthetic smooth bump/decaying spectra self-convolved at a range of
//    sizes (n) and peak magnitudes (scale-invariance), to establish the
//    eps*sqrt(n) scaling itself.
//  * Real order-1 (G1) spectra generated LIVE from the actual production
//    VDOSEval/VDOSGn machinery (several Debye models and the real Li2O
//    material, matching app_gnerange's approach), self-convolved via
//    FastConvolve directly here (order 1 is never truncated in production,
//    so gn.getRawSpectrum(Order{1}) already gives the exact untruncated
//    input FastConvolve itself sees) to reproduce orders 2 and 3's raw,
//    pre-truncation FFT output exactly as produceNewOrderByConvolutionImpl
//    computes it, but with access to the ground truth besides.
//
// Deliberately verbose/exploratory, in the same spirit as app_gnerange.
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/utils/NCFastConvolve.hh"
#include "NCrystal/internal/vdos/NCVDOSUtils.hh"
#include "NCrystal/internal/vdos/NCVDOSGn.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <iostream>
#include <sstream>

namespace NC = NCrystal;
namespace NCV = NCrystal::VDOS;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  //Exact, FMA-sensitive per-point error values are only printed under this
  //(not part of the routinely cross-platform-compared output -- see
  //checkNoiseFloor below):
  const bool s_verbose = NC::ncgetenv_bool("GNCONVNOISEFLOOR_VERBOSE");

  // O(n^2) direct linear convolution (Neumaier-summed), same {a1,a2,dt}
  // convention as FastConvolve::convolve. Accurate to ~1ULP regardless of n,
  // so serves as ground truth for measuring FastConvolve's own FFT error.
  NC::VectD referenceConvolve( const NC::VectD& a1, const NC::VectD& a2, double dt )
  {
    const std::size_t n1 = a1.size();
    const std::size_t n2 = a2.size();
    nc_assert_always( n1 > 0 && n2 > 0 );
    NC::VectD y( n1 + n2 - 1, 0.0 );
    for ( auto k : NC::ncrange(y.size()) ) {
      const std::size_t ilo = ( k + 1 >= n2 ) ? k + 1 - n2 : std::size_t(0);
      const std::size_t ihi = NC::ncmin( k, n1 - 1 );
      NC::StableSum sum;
      for ( auto i : NC::ncrange(ilo,ihi+1) )
        sum.add( a1.at(i) * a2.at(k-i) );
      y[k] = dt * sum.sum();
    }
    return y;
  }

  // Convolves a1,a2 both via FastConvolve and via referenceConvolve, and
  // *asserts* (REQUIRE, not just a printed value) that FastConvolve's true
  // max absolute error stays within estimateFFTConvolutionNoiseFloor's
  // bound at the given safetyFactor -- this is the actual property the
  // whole utility exists to guarantee, so it must be a real, hard check
  // that fails the test outright on a genuine regression, not something
  // that could be silently papered over by a maintainer running --update
  // on a log-diff mismatch.
  //
  // maxErr and the local error right at the production truncation edge
  // (the first/last point exceeding prodThreshold*peak=1e-13*peak, i.e.
  // exactly the point that decides where produceNewOrderByConvolutionImpl
  // trims the spectrum today) are themselves FastConvolve's own round-off
  // -- by construction, exactly the quantity that legitimately varies with
  // FMA contraction (build flags/platform/compiler). Printing their exact
  // value in the routinely-diffed (committed, cross-platform-compared)
  // output would make the test fail on every platform that contracts
  // differently, even though the REQUIRE above already correctly confirms
  // nothing is actually wrong -- a false positive, not a real one. So the
  // routine output only contains quantities that are themselves
  // deterministic/reproducible (label/n/peak/noiseFloor/edgeVal and the
  // qualitative verdict); the exact noisy values are still available for
  // interactive debugging via GNCONVNOISEFLOOR_VERBOSE=1. See
  // docs/claude_session_vdos_fma_reprod.md.
  double checkNoiseFloor( const std::string& label, const NC::VectD& a1,
                          const NC::VectD& a2, double dt, double safetyFactor )
  {
    NC::FastConvolve fc;
    NC::VectD y_fft;
    fc.convolve( a1, a2, y_fft, dt );
    const NC::VectD y_ref = referenceConvolve( a1, a2, dt );
    nc_assert_always( y_fft.size() == y_ref.size() );
    const double peak = *std::max_element( y_ref.begin(), y_ref.end() );
    double maxErr = 0.0;
    for ( auto i : NC::ncrange(y_fft.size()) )
      maxErr = NC::ncmax( maxErr, NC::ncabs( y_fft.at(i) - y_ref.at(i) ) );
    const double noiseFloor = NCV::estimateFFTConvolutionNoiseFloor(
      peak, y_fft.size(), safetyFactor );
    const double ratio = maxErr / NC::ncmax( noiseFloor, 1e-300 );
    REQUIRE( ratio <= 1.0 );//the actual property under test

    constexpr double prodThreshold = 1e-13;
    const double prodCutoff = prodThreshold * peak;
    std::size_t edgeIdx = y_ref.size();
    for ( auto i : NC::ncrange(y_ref.size()) ) {
      if ( y_ref.at(i) > prodCutoff ) {
        edgeIdx = i;
        break;
      }
    }
    double edgeErr = -1.0, edgeVal = -1.0, edgeRelErr = -1.0;
    if ( edgeIdx < y_ref.size() ) {
      edgeVal = y_ref.at(edgeIdx);
      edgeErr = NC::ncabs( y_fft.at(edgeIdx) - edgeVal );
      edgeRelErr = edgeErr / NC::ncmax( edgeVal, 1e-300 );
    }

    std::cout << "  " << label << ": n=" << y_fft.size()
              << " peak=" << NC::fmt(peak,"%.6g")
              << " noiseFloor=" << NC::fmt(noiseFloor,"%.3e")
              << "  OK" << std::endl;
    std::cout << "      at prod. truncation edge (val=" << NC::fmt(edgeVal,"%.3e")
              << ", =" << NC::fmt(edgeVal/NC::ncmax(peak,1e-300),"%.3g") << "*peak): "
              << "noiseFloor/edgeVal=" << NC::fmt(noiseFloor/NC::ncmax(edgeVal,1e-300),"%.3g")
              << std::endl;
    if ( s_verbose ) {
      std::cout << "      [verbose, not cross-platform-stable] maxErr="
                << NC::fmt(maxErr,"%.3e")
                << " maxErr/noiseFloor=" << NC::fmt(ratio,"%.3g")
                << " edgeErr=" << NC::fmt(edgeErr,"%.3e")
                << " edgeRelErr=" << NC::fmt(edgeRelErr,"%.3g")
                << std::endl;
    }
    return ratio;
  }

  //////////////////////////////////////////////////////////
  // Synthetic case: establishing/checking the sqrt(n) scaling, //
  // and scale-invariance in the overall peak magnitude.         //
  //////////////////////////////////////////////////////////

  NC::VectD gaussianBump( std::size_t n, double k, double peak )
  {
    NC::VectD v(n);
    const double mid = 0.5*(n-1);
    for ( auto i : NC::ncrange(n) ) {
      const double x = (static_cast<double>(i)-mid)/mid;
      v[i] = peak*std::exp(-k*x*x);
    }
    return v;
  }

  void syntheticSweep( double safetyFactor )
  {
    std::cout << "--- Synthetic self-convolved Gaussian bumps:"
                 " varying n ---" << std::endl;
    double worst = 0.0;
    for ( std::size_t n : { std::size_t(200), std::size_t(500),
                            std::size_t(1000), std::size_t(2000),
                            std::size_t(4000) } ) {
      auto a = gaussianBump( n, 8.0, 1.0 );
      std::ostringstream oss;
      oss << "gaussian(n=" << n << ",peak=1)";
      worst = NC::ncmax( worst, checkNoiseFloor( oss.str(), a, a, 1.0, safetyFactor ) );
    }
    std::cout << "--- Synthetic self-convolved Gaussian bumps:"
                 " varying overall scale (n=2000 fixed) ---" << std::endl;
    for ( double peak : { 1e-8, 1e-3, 1.0, 1e5, 1e12 } ) {
      auto a = gaussianBump( 2000, 8.0, peak );
      std::ostringstream oss;
      oss << "gaussian(n=2000,peak=" << NC::fmtg(peak) << ")";
      worst = NC::ncmax( worst, checkNoiseFloor( oss.str(), a, a, 1.0, safetyFactor ) );
    }
    if ( s_verbose )
      std::cout << "      [verbose, not cross-platform-stable] worst"
                   " maxErr/noiseFloor over synthetic sweep: "
                << NC::fmt(worst,"%.3g") << std::endl;
    std::cout << std::endl;
  }

  //////////////////////////////////////////////////////////////////////
  // Real data: order 1 (G1) generated LIVE from production VDOSEval/  //
  // VDOSGn, self-convolved here via FastConvolve directly to reproduce //
  // orders 2 and 3's exact raw, pre-truncation FFT output the same way //
  // produceNewOrderByConvolutionImpl does internally (order 1 itself is //
  // never truncated in production, so this is the exact untruncated    //
  // input FastConvolve sees there too).                                 //
  //////////////////////////////////////////////////////////////////////

  void realDataCase( const char* label, const NC::VDOSData& vdosdata,
                     double safetyFactor, double& worst )
  {
    NC::VDOSEval ve(vdosdata);
    NCV::VDOSGn gn( ve, NCV::VDOSGn::Cfg::Default );
    const auto& g1 = gn.getRawSpectrum( NCV::VDOSGn::Order(1) );
    const double dt = gn.binWidth( NCV::VDOSGn::Order(1) );
    std::ostringstream oss2;
    oss2 << label << ", order2=G1(*)G1";
    worst = NC::ncmax( worst, checkNoiseFloor( oss2.str(), g1, g1, dt, safetyFactor ) );
    NC::FastConvolve fc;
    NC::VectD order2raw;
    fc.convolve( g1, g1, order2raw, dt );
    std::ostringstream oss3;
    oss3 << label << ", order3=G1(*)order2raw";
    worst = NC::ncmax( worst, checkNoiseFloor( oss3.str(), g1, order2raw, dt, safetyFactor ) );
  }

  NC::VDOSData makeDebyeVDOS( double debyeTempK, double tempK )
  {
    return NC::createVDOSDebye( NC::DebyeTemperature{debyeTempK},
                                NC::Temperature{tempK},
                                NC::SigmaBound{1.5},
                                NC::AtomMass{26.98} );
  }

  void realDataSweep( double safetyFactor )
  {
    std::cout << "--- Real data: live Debye-model G1, self-convolved"
                 " (orders 2,3) ---" << std::endl;
    double worst = 0.0;
    realDataCase( "Debye T=80K,temp=10K (narrow/steep, cold)",
                 makeDebyeVDOS(80.0,10.0), safetyFactor, worst );
    realDataCase( "Debye T=410K,temp=293.15K (Al-like)",
                 makeDebyeVDOS(410.0,293.15), safetyFactor, worst );
    realDataCase( "Debye T=2000K,temp=293.15K (wide)",
                 makeDebyeVDOS(2000.0,293.15), safetyFactor, worst );

    std::cout << "--- Real data: live Li2O (Li,O) material G1,"
                 " self-convolved (orders 2,3) ---" << std::endl;
    NC::MatCfg cfg( "Li2O_sg225_LithiumOxide.ncmat;vdoslux=4;temp=10K" );
    auto info = NC::FactImpl::createInfo( cfg );
    for ( auto& di : info->getDynamicInfoList() ) {
      auto di_vdos = dynamic_cast<const NC::DI_VDOS*>( di.get() );
      if ( !di_vdos )
        continue;
      std::ostringstream oss;
      oss << "Li2O (real material, temp=10K), element="
          << di->atomData().elementName();
      realDataCase( oss.str().c_str(), di_vdos->vdosData(), safetyFactor, worst );
    }
    if ( s_verbose )
      std::cout << "      [verbose, not cross-platform-stable] worst"
                   " maxErr/noiseFloor over real-data sweep: "
                << NC::fmt(worst,"%.3g") << std::endl;
    std::cout << std::endl;
  }
}

int main() {
  //Candidate safety factor: chosen generously above 1.0 so the bound holds
  //with headroom (not just barely) across every case checked below; see the
  //printed margins for how conservative it ends up being in practice.
  const double safetyFactor = 8.0;

  syntheticSweep( safetyFactor );
  realDataSweep( safetyFactor );

  return 0;
}

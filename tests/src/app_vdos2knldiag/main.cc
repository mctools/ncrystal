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
// Regression-test dump of NC::VDOS::expandVDOSToGnFcts, NC::VDOS::           //
// createScatteringKernel and NC::SABUtils::transformKernelToStdFormat --    //
// the pipeline the C API's ncrystal_raw_vdos2kernel (and hence              //
// NCrystal.vdos.extractKnl, and hence ncmat2endf) uses to expand a VDOS     //
// into a full S(alpha,beta) kernel -- for the ThO2 material's Thorium and   //
// Oxygen elements (vdoslux=1, requested emax=1eV, T=293.15K).              //
//
// Originally written much more verbosely (full alpha/beta/sab grids, every //
// phonon order's raw Gn diagnostics, at full %.17g precision) to bisect a   //
// cross-platform/compiler reproducibility divergence -- see                //
// docs/claude_session_vdos_fma_reprod.md for that investigation, since     //
// resolved. Trimmed back down to summary statistics (size/front/back/min/   //
// max/sum rather than full array dumps) at a more modest precision, once   //
// no longer needed for active bisection.                                   //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/vdos/NCVDOSExpand.hh"
#include "NCrystal/internal/vdos/NCVDOSToScatKnl.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <iostream>
#include <sstream>
#include <string>

namespace NC = NCrystal;

#define REQUIRE(x) nc_assert_always(x)

namespace {

  //Full ("%.10g") precision is only printed under this: a decaying curve's
  //tail (e.g. sab far from its peak) can sit many orders of magnitude below
  //the array's own max, where ordinary FMA-contraction noise -- utterly
  //negligible in absolute terms -- shows up as an enormous RELATIVE
  //difference (confirmed locally against -mfma: up to ~36% for this file's
  //Th "knl sab" min). So even a coarsened relative precision cannot make the
  //routinely cross-platform-compared default output stable; the REQUIRE
  //checks below are what actually guard against a real regression. See
  //docs/claude_session_vdos_fma_reprod.md.
  const bool s_verbose = NC::ncgetenv_bool("VDOS2KNLDIAG_VERBOSE");
  const char* fmtprec() { return s_verbose ? "%.10g" : "%.4g"; }

  //Summary statistics rather than a full element-by-element dump, so the
  //printed size does not scale with the (potentially huge) array size:
  void dumpVectD( const char* label, const NC::Span<const double> v )
  {
    std::cout << label << " [" << v.size() << "]";
    if ( v.empty() ) {
      std::cout << std::endl;
      return;
    }
    NC::StableSum sum;
    double mn = v.front(), mx = v.front();
    for ( auto x : v ) {
      sum.add(x);
      mn = NC::ncmin(mn,x);
      mx = NC::ncmax(mx,x);
    }
    const double s = sum.sum();
    REQUIRE( std::isfinite(v.front()) && std::isfinite(v.back()) );
    REQUIRE( std::isfinite(mn) && std::isfinite(mx) && std::isfinite(s) );
    REQUIRE( mn <= mx );
    REQUIRE( v.front() >= mn && v.front() <= mx );
    REQUIRE( v.back() >= mn && v.back() <= mx );
    //Print deep-tail-noise-prone values (front/back/min, which for a curve
    //like sab legitimately approach zero) as a stable "~0" once negligible
    //relative to the array's own peak magnitude, rather than at a fixed
    //relative precision that a noise-floor value cannot honour:
    const double magref = NC::ncmax( NC::ncabs(mn), NC::ncabs(mx) );
    auto fmtTail = [&]( double x ) -> std::string
    {
      if ( magref > 0.0 && NC::ncabs(x) < 1e-8*magref )
        return "~0";
      std::ostringstream oss;
      oss << NC::fmt(x,fmtprec());
      return oss.str();
    };
    std::cout << ": front=" << fmtTail(v.front())
              << " back=" << fmtTail(v.back())
              << " min=" << fmtTail(mn)
              << " max=" << NC::fmt(mx,fmtprec())
              << " sum=" << NC::fmt(s,fmtprec())
              << std::endl;
  }

  //Real invariants (not just "the same as last time") for a Rectangle
  //diagnostic value, catching e.g. a NaN or empty-range regression while
  //saying nothing about the FMA-sensitive exact coordinates:
  void checkRect( const NC::Rectangle& r )
  {
    REQUIRE( !r.isEmpty() );
    REQUIRE( std::isfinite(r.x0()) && std::isfinite(r.x1()) );
    REQUIRE( std::isfinite(r.y0()) && std::isfinite(r.y1()) );
  }

  void dumpElement( const NC::DI_VDOS& di_vdos, const NC::DynamicInfo& di )
  {
    std::cout << "=== Element: " << di.atomData().elementName()
              << " (mass=" << NC::fmt(di.atomData().averageMassAMU().dbl(),fmtprec())
              << " amu, temperature=" << NC::fmt(di.temperature().dbl(),fmtprec())
              << "K) ===" << std::endl;

    //Build the VDOSData the same way ncc::createVDOSDataFromRaw does in
    //ncrystal.cc (regularising the ORIGINAL, pre-regularisation VDOS curve
    //fresh), rather than reusing DI_VDOS::vdosData() (regularised once
    //already, at .ncmat load time): this is the exact pipeline
    //NCrystal.vdos.extractKnl / ncrystal_raw_vdos2kernel / ncmat2endf use:
    NC::VectD regEgrid, regDensity;
    std::tie( regEgrid, regDensity )
      = NC::regulariseVDOSGrid( di_vdos.vdosOrigEgrid(), di_vdos.vdosOrigDensity() );
    nc_assert_always( regEgrid.size() == 2 );
    NC::VDOSData vdosdata( NC::PairDD( regEgrid.front(), regEgrid.back() ),
                           std::move(regDensity), di.temperature(),
                           di.atomData().scatteringXS(),
                           di.atomData().averageMassAMU() );
    std::cout << "vdos_egrid: " << NC::fmt(vdosdata.vdos_egrid().first,fmtprec())
              << ' ' << NC::fmt(vdosdata.vdos_egrid().second,fmtprec()) << std::endl;
    dumpVectD( "vdos_density", vdosdata.vdos_density() );

    const NC::VDOS::VDOSLux vdoslux( 1 );
    NC::Optional<NC::NeutronEnergy> targetEmax;
    targetEmax = NC::NeutronEnergy{ 1.0 };

    auto gnexpn = NC::VDOS::expandVDOSToGnFcts( vdosdata, vdoslux, targetEmax );
    REQUIRE( gnexpn.Gn.maxOrder().value() >= 1 );
    REQUIRE( std::isfinite(gnexpn.alpha2x) && gnexpn.alpha2x > 0.0 );
    REQUIRE( std::isfinite(gnexpn.suggestedEmax.dbl()) && gnexpn.suggestedEmax.dbl() > 0.0 );
    checkRect( gnexpn.sabRange );
    std::cout << "gnexpn maxOrder: " << gnexpn.Gn.maxOrder().value() << std::endl;
    std::cout << "gnexpn alpha2x: " << NC::fmt(gnexpn.alpha2x,fmtprec()) << std::endl;
    std::cout << "gnexpn suggestedEmax: "
              << NC::fmt(gnexpn.suggestedEmax.dbl(),fmtprec()) << std::endl;
    std::cout << "gnexpn sabRange: "
              << NC::fmt(gnexpn.sabRange.x0(),fmtprec()) << ' '
              << NC::fmt(gnexpn.sabRange.x1(),fmtprec()) << ' '
              << NC::fmt(gnexpn.sabRange.y0(),fmtprec()) << ' '
              << NC::fmt(gnexpn.sabRange.y1(),fmtprec()) << std::endl;
    //Final order's (alpha,beta) range (before intersection with the
    //kinematically accessible region) -- the per-order table this used to
    //print in full was only needed while actively bisecting which specific
    //order first diverged between platforms; the final order is what
    //determines the downstream result:
    const auto& lastAB = gnexpn.abRanges.back();
    checkRect( lastAB );
    std::cout << "gnexpn last abRange: n=" << gnexpn.abRanges.size()
              << ' ' << NC::fmt(lastAB.x0(),fmtprec())
              << ' ' << NC::fmt(lastAB.x1(),fmtprec())
              << ' ' << NC::fmt(lastAB.y0(),fmtprec())
              << ' ' << NC::fmt(lastAB.y1(),fmtprec()) << std::endl;

    //Raw (unfiltered, pre-relcontriblvl-threshold) G1 spectrum (the base
    //single-phonon spectrum, before any convolution):
    dumpVectD( "gnexpn G1 rawspec", gnexpn.Gn.getRawSpectrum( NC::VDOS::VDOSGn::Order(1u) ) );

    auto knldata = NC::VDOS::createScatteringKernel( vdosdata, vdoslux, targetEmax );
    dumpVectD( "knl alphaGrid", knldata.alphaGrid );
    dumpVectD( "knl betaGrid", knldata.betaGrid );
    dumpVectD( "knl sab", knldata.sab );
    REQUIRE( std::isfinite(knldata.suggestedEmax) && knldata.suggestedEmax > 0.0 );
    std::cout << "knl suggestedEmax: "
              << NC::fmt(knldata.suggestedEmax,fmtprec()) << std::endl;

    auto sabdata = NC::SABUtils::transformKernelToStdFormat( std::move(knldata) );
    dumpVectD( "std alphaGrid", sabdata.alphaGrid() );
    dumpVectD( "std betaGrid", sabdata.betaGrid() );
    dumpVectD( "std sab", sabdata.sab() );
    REQUIRE( std::isfinite(sabdata.suggestedEmax()) && sabdata.suggestedEmax() > 0.0 );
    std::cout << "std suggestedEmax: "
              << NC::fmt(sabdata.suggestedEmax(),fmtprec()) << std::endl;
    REQUIRE( sabdata.alphaGrid().front() <= sabdata.alphaGrid().back() );
    REQUIRE( sabdata.betaGrid().front() <= sabdata.betaGrid().back() );
    std::cout << "alpha range: " << NC::fmt(sabdata.alphaGrid().front(),fmtprec())
              << ' ' << NC::fmt(sabdata.alphaGrid().back(),fmtprec()) << std::endl;
    std::cout << "beta range: " << NC::fmt(sabdata.betaGrid().front(),fmtprec())
              << ' ' << NC::fmt(sabdata.betaGrid().back(),fmtprec()) << std::endl;
  }
}

int main()
{
  NC::MatCfg cfg( "ThO2_sg225_ThoriumDioxide.ncmat;vdoslux=1;temp=293.15K" );
  auto info = NC::FactImpl::createInfo( cfg );

  for ( auto& di : info->getDynamicInfoList() ) {
    auto di_vdos = dynamic_cast<const NC::DI_VDOS*>( di.get() );
    if ( di_vdos )
      dumpElement( *di_vdos, *di );
  }
  return 0;
}

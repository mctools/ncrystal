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
// Diagnostic dump of NC::VDOS::expandVDOSToGnFcts, NC::VDOS::                //
// createScatteringKernel and NC::SABUtils::transformKernelToStdFormat --    //
// the pipeline the C API's ncrystal_raw_vdos2kernel (and hence              //
// NCrystal.vdos.extractKnl, and hence ncmat2endf) uses to expand a VDOS     //
// into a full S(alpha,beta) kernel. Written to help bisect a                //
// cross-platform/compiler reproducibility divergence seen in                //
// tests/scripts/n2endf_bad.py, for the ThO2 material's Thorium element      //
// (vdoslux=1, requested emax=1eV, T=293.15K) -- see that test for the       //
// original symptom, and doc/claude_session_vdos_fma_reprod.md for the       //
// investigation so far. Also dumps Oxygen (same material) as a same-run,   //
// known-good comparison point.                                             //
//                                                                            //
// In particular dumps GnExpansion::abRanges -- the per-phonon-order         //
// (alpha,beta) range computed by expandVDOSToGnFcts's order-growth loop,    //
// before intersection with the kinematically accessible region -- so a     //
// divergence in the *number* of orders grown (seen as a shifted final      //
// alpha/beta range) can be bisected down to which order, and whether its   //
// alpha or beta bound, first differs between platforms/compilers.          //
//                                                                            //
// Deliberately very verbose (full alpha/beta/sab grids, plus every order's  //
// abRange) for bisecting a divergence report from CI without needing the    //
// failing compiler locally; trim once the underlying bug is found and      //
// fixed.                                                                    //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/vdos/NCVDOSExpand.hh"
#include "NCrystal/internal/vdos/NCVDOSToScatKnl.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"
#include <iostream>

namespace NC = NCrystal;

namespace {
  void dumpVectD( const char* label, const NC::VectD& v )
  {
    std::cout << label << " [" << v.size() << "]:";
    for ( auto x : v )
      std::cout << ' ' << NC::fmt(x,"%.17g");
    std::cout << std::endl;
  }

  void dumpElement( const NC::DI_VDOS& di_vdos, const NC::DynamicInfo& di )
  {
    std::cout << "=== Element: " << di.atomData().elementName()
              << " (mass=" << NC::fmt(di.atomData().averageMassAMU().dbl(),"%.17g")
              << " amu, temperature=" << NC::fmt(di.temperature().dbl(),"%.17g")
              << "K) ===" << std::endl;

    //Build the VDOSData the same way ncc::createVDOSDataFromRaw does in
    //ncrystal.cc (regularising the ORIGINAL, pre-regularisation VDOS
    //curve fresh), rather than reusing DI_VDOS::vdosData() (regularised
    //once already, at .ncmat load time): this is the exact pipeline
    //NCrystal.vdos.extractKnl / ncrystal_raw_vdos2kernel / ncmat2endf /
    //tests/scripts/n2endf_bad.py use, and it turned out NOT to be
    //equivalent to vdosData() for reproducibility purposes -- see
    //docs/claude_session_vdos_fma_reprod.md for how this was found (a
    //CI round showed this app's own alpha/beta range matching exactly
    //while n2endf_bad.py still diverged, using the vdosData()-based
    //version of this app):
    NC::VectD regEgrid, regDensity;
    std::tie( regEgrid, regDensity )
      = NC::regulariseVDOSGrid( di_vdos.vdosOrigEgrid(), di_vdos.vdosOrigDensity() );
    nc_assert_always( regEgrid.size() == 2 );
    NC::VDOSData vdosdata( NC::PairDD( regEgrid.front(), regEgrid.back() ),
                           std::move(regDensity), di.temperature(),
                           di.atomData().scatteringXS(),
                           di.atomData().averageMassAMU() );
    std::cout << "vdos_egrid: " << NC::fmt(vdosdata.vdos_egrid().first,"%.17g")
              << ' ' << NC::fmt(vdosdata.vdos_egrid().second,"%.17g") << std::endl;
    dumpVectD( "vdos_density", vdosdata.vdos_density() );

    const NC::VDOS::VDOSLux vdoslux( 1 );
    NC::Optional<NC::NeutronEnergy> targetEmax;
    targetEmax = NC::NeutronEnergy{ 1.0 };

    auto gnexpn = NC::VDOS::expandVDOSToGnFcts( vdosdata, vdoslux, targetEmax );
    std::cout << "gnexpn maxOrder: " << gnexpn.Gn.maxOrder().value() << std::endl;
    std::cout << "gnexpn alpha2x: " << NC::fmt(gnexpn.alpha2x,"%.17g") << std::endl;
    std::cout << "gnexpn suggestedEmax: "
              << NC::fmt(gnexpn.suggestedEmax.dbl(),"%.17g") << std::endl;
    std::cout << "gnexpn sabRange: "
              << NC::fmt(gnexpn.sabRange.x0(),"%.17g") << ' '
              << NC::fmt(gnexpn.sabRange.x1(),"%.17g") << ' '
              << NC::fmt(gnexpn.sabRange.y0(),"%.17g") << ' '
              << NC::fmt(gnexpn.sabRange.y1(),"%.17g") << std::endl;
    std::cout << "gnexpn abRanges [" << gnexpn.abRanges.size() << "] (n alpha0 alpha1 beta0 beta1):" << std::endl;
    for ( auto i : NC::ncrange(gnexpn.abRanges.size()) ) {
      const auto& r = gnexpn.abRanges.at(i);
      std::cout << "  " << (i+1) << ' '
                << NC::fmt(r.x0(),"%.17g") << ' '
                << NC::fmt(r.x1(),"%.17g") << ' '
                << NC::fmt(r.y0(),"%.17g") << ' '
                << NC::fmt(r.y1(),"%.17g") << std::endl;
    }

    //Raw (unfiltered, pre-relcontriblvl-threshold) per-order Gn diagnostics,
    //to see exactly which order's *construction* (not just its thresholded
    //abRange) first diverges -- in particular binWidth reveals when the
    //on-demand thinning path in produceNewOrderByConvolutionImpl kicks in
    //(see NCVDOSGn.cc), and specSize reveals truncation-boundary shifts:
    std::cout << "gnexpn raw Gn [" << gnexpn.Gn.maxOrder().value()
              << "] (n binWidth eRangeLo eRangeHi specSize spec.front spec.back):" << std::endl;
    for ( auto n : NC::ncrange(1u,gnexpn.Gn.maxOrder().value()+1) ) {
      NC::VDOS::VDOSGn::Order order( n );
      auto erange = gnexpn.Gn.eRange( order );
      const auto& spec = gnexpn.Gn.getRawSpectrum( order );
      std::cout << "  " << n << ' '
                << NC::fmt(gnexpn.Gn.binWidth(order),"%.17g") << ' '
                << NC::fmt(erange.first,"%.17g") << ' '
                << NC::fmt(erange.second,"%.17g") << ' '
                << spec.size() << ' '
                << NC::fmt(spec.front(),"%.17g") << ' '
                << NC::fmt(spec.back(),"%.17g") << std::endl;
    }

    //Full raw G1 spectrum (the base single-phonon spectrum, before any
    //convolution at all): if this is already bit-different between
    //platforms, the noise enters in VDOSEval's spectrum construction
    //(std::exp/tanh/sinh calls, inherently libm-version-sensitive); if it
    //is bit-identical, the noise must enter during the G1(x)G1 -> G2
    //self-convolution (order 2) itself:
    dumpVectD( "gnexpn G1 rawspec", gnexpn.Gn.getRawSpectrum( NC::VDOS::VDOSGn::Order(1u) ) );
    //Full raw G2 spectrum (G1(x)G1, the very first actual convolution):
    dumpVectD( "gnexpn G2 rawspec", gnexpn.Gn.getRawSpectrum( NC::VDOS::VDOSGn::Order(2u) ) );

    auto knldata = NC::VDOS::createScatteringKernel( vdosdata, vdoslux, targetEmax );
    dumpVectD( "knl alphaGrid", knldata.alphaGrid );
    dumpVectD( "knl betaGrid", knldata.betaGrid );
    dumpVectD( "knl sab", knldata.sab );
    std::cout << "knl suggestedEmax: "
              << NC::fmt(knldata.suggestedEmax,"%.17g") << std::endl;

    auto sabdata = NC::SABUtils::transformKernelToStdFormat( std::move(knldata) );
    dumpVectD( "std alphaGrid", sabdata.alphaGrid() );
    dumpVectD( "std betaGrid", sabdata.betaGrid() );
    dumpVectD( "std sab", sabdata.sab() );
    std::cout << "std suggestedEmax: "
              << NC::fmt(sabdata.suggestedEmax(),"%.17g") << std::endl;
    std::cout << "alpha range: " << NC::fmt(sabdata.alphaGrid().front(),"%.17g")
              << ' ' << NC::fmt(sabdata.alphaGrid().back(),"%.17g") << std::endl;
    std::cout << "beta range: " << NC::fmt(sabdata.betaGrid().front(),"%.17g")
              << ' ' << NC::fmt(sabdata.betaGrid().back(),"%.17g") << std::endl;
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

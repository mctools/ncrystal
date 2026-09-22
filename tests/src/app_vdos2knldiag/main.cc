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
// Diagnostic dump of NC::VDOS::createScatteringKernel +                     //
// NC::SABUtils::transformKernelToStdFormat -- the exact same two calls the  //
// C API's ncrystal_raw_vdos2kernel (and hence NCrystal.vdos.extractKnl, and //
// hence ncmat2endf) uses to expand a VDOS into a full S(alpha,beta) kernel. //
// Written to help bisect a gcc-10-only reproducibility divergence seen in   //
// tests/scripts/n2endf_bad.py, for the ThO2 material's Thorium element      //
// (vdoslux=1, requested emax=1eV, T=293.15K) -- see that test for the       //
// original symptom. Also dumps Oxygen (same material) as a same-run,       //
// known-good comparison point.                                             //
//                                                                            //
// Deliberately very verbose (full alpha/beta/sab grids) for bisecting a     //
// divergence report from CI without needing the failing compiler locally;   //
// trim once the underlying bug is found and fixed.                         //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/interfaces/NCInfo.hh"
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

    const auto& vdosdata = di_vdos.vdosData();
    std::cout << "vdos_egrid: " << NC::fmt(vdosdata.vdos_egrid().first,"%.17g")
              << ' ' << NC::fmt(vdosdata.vdos_egrid().second,"%.17g") << std::endl;
    dumpVectD( "vdos_density", vdosdata.vdos_density() );

    const NC::VDOS::VDOSLux vdoslux( 1 );
    NC::Optional<NC::NeutronEnergy> targetEmax;
    targetEmax = NC::NeutronEnergy{ 1.0 };
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

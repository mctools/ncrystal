#ifndef NCrystal_SABCfg_hh
#define NCrystal_SABCfg_hh

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

#include "NCrystal/internal/utils/NCStrView.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SABCfg {

    //////////////////////////////////////////////////////////////////////
    // Utilities for configuring the various settings related to SAB    //
    // processing. In particular this includes initialisation of values //
    // based on the user-specified luxury level.                        //
    //////////////////////////////////////////////////////////////////////

    /////////////////////////////////////////////
    // Supported numerical integration schemes //
    /////////////////////////////////////////////

    enum class IntegrationScheme : std::uint_fast32_t {
      //Basic fixed-order schemes:
      Trapez2   = 0x020002,
      Trapez3   = 0x020003,
      Trapez5   = 0x020005,
      Trapez9   = 0x020009,
      Trapez17  = 0x020011,
      Trapez33  = 0x020021,
      Simpson3  = 0x040003,
      Simpson5  = 0x040005,
      Simpson9  = 0x040009,
      Simpson17 = 0x040011,
      Simpson33 = 0x040021,
      Romberg5  = 0x010005,
      Romberg9  = 0x010009,
      Romberg17 = 0x010011,
      Romberg33 = 0x010021,
      //Flex schemes, based on potentially adaptive Romberg integration:
      Flex5  = 0x100005,//mostly Romberg5, occasionally adaptive Romberg17+
      //with target prec=1e-3
      Flex9  = 0x100009,//mostly Romberg9, occasionally adaptive Romberg17+
      //with target prec=1e-4
      Flex17 = 0x100011,//mostly Romberg17, occasionally adaptive Romberg17+
      //with target prec=1e-6
      Flex33 = 0x100021,//Adaptive Romberg33+, target prec=1e-9
      Flex65 = 0x100041,//Adaptive Romberg65+, target prec=1e-12
      MaxPrec = Flex65,
      Default = Flex17//fixme revisit (or remove?)
    };

    //Utilities for encoding to/from strings
    const char * integSchemeToStr( IntegrationScheme );
    IntegrationScheme str2IntegScheme( StrView );
    const char * allIntegSchemesAsStr();//';' separated list

    ////////////////
    // Cfg object //
    ////////////////

    struct Cfg {
      //Cfg object. Note that the default values correspond to sablux=3.
      //fixme: consistent snake case vs. camel case here
      IntegrationScheme integScheme = IntegrationScheme::Flex9;
      IntegrationScheme integSchemeBCSample = IntegrationScheme::Flex9;
      unsigned egrid_npts = 300;
      double egrid_emin_accuracy = 0.01;
      double fullCellSamplingARThreshold = 0.15;
      double bcSamplingLargeSRatioThreshold = 1e-6;
    };

    ////////////////////////////////////////////
    // Factory function based on luxury level //
    ////////////////////////////////////////////

    Cfg createConfig( int sablux );

  }
}

#endif

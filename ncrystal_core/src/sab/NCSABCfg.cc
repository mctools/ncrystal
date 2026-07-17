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

#include "NCrystal/internal/sab/NCSABCfg.hh"

namespace NC = NCrystal;

NC::SABCfg::IntegrationScheme NC::SABCfg::str2IntegScheme( StrView v )
{
  using IS = IntegrationScheme;
  if ( v.startswith("Flex") ) {
    if ( v=="Flex9" )
      return IS::Flex9;
    if ( v=="Flex5" )
      return IS::Flex5;
    if ( v=="Flex17" )
      return IS::Flex17;
    if ( v=="Flex33" )
      return IS::Flex33;
    if ( v=="Flex65" )
      return IS::Flex65;
  } else if ( v.startswith("Romberg") ) {
    if ( v=="Romberg9" )
      return IS::Romberg9;
    if ( v=="Romberg5" )
      return IS::Romberg5;
    if ( v=="Romberg17" )
      return IS::Romberg17;
    if ( v=="Romberg33" )
      return IS::Romberg33;
  } else if ( v.startswith("Trapez") ) {
    if ( v=="Trapez2" )
      return IS::Trapez2;
    if ( v=="Trapez3" )
      return IS::Trapez3;
    if ( v=="Trapez5" )
      return IS::Trapez5;
    if ( v=="Trapez9" )
      return IS::Trapez9;
    if ( v=="Trapez17" )
      return IS::Trapez17;
    if ( v=="Trapez33" )
      return IS::Trapez33;
  } else if ( v.startswith("Simpson") ) {
    if ( v=="Simpson3" )
      return IS::Simpson3;
    if ( v=="Simpson5" )
      return IS::Simpson5;
    if ( v=="Simpson9" )
      return IS::Simpson9;
    if ( v=="Simpson17" )
      return IS::Simpson17;
    if ( v=="Simpson33" )
      return IS::Simpson33;
  }
  NCRYSTAL_THROW2(BadInput,"Invalid integration scheme: \""<<v
                  <<"\" (should be one of \""
                  << allIntegSchemesAsStr() <<"\")");
  return IS::Default;
}

const char * NC::SABCfg::allIntegSchemesAsStr()
{
  return "Flex5;Flex9;Flex17;Flex33;Flex65;"
    "Romberg5;Romberg9;Romberg17;Romberg33;"
    "Trapez2;Trapez3;Trapez5;Trapez9;Trapez17;Trapez33;"
    "Simpson3;Simpson5;Simpson9;Simpson17;Simpson33";
}

const char * NC::SABCfg::integSchemeToStr( IntegrationScheme v )
{
  using IS = IntegrationScheme;
  switch ( v ) {
    //fixme: something shorter, so might be used in cfg strings? r33?
  case IS::Flex5:  return "Flex5";
  case IS::Flex9:  return "Flex9";
  case IS::Flex17: return "Flex17";
  case IS::Flex33: return "Flex33";
  case IS::Flex65: return "Flex65";
  case IS::Romberg5:  return "Romberg5";
  case IS::Romberg9:  return "Romberg9";
  case IS::Romberg17: return "Romberg17";
  case IS::Romberg33: return "Romberg33";
  case IS::Trapez2:   return "Trapez2";
  case IS::Trapez3:   return "Trapez3";
  case IS::Trapez5:   return "Trapez5";
  case IS::Trapez9:   return "Trapez9";
  case IS::Trapez17:  return "Trapez17";
  case IS::Trapez33:  return "Trapez33";
  case IS::Simpson3:  return "Simpson3";
  case IS::Simpson5:  return "Simpson5";
  case IS::Simpson9:  return "Simpson9";
  case IS::Simpson17: return "Simpson17";
  case IS::Simpson33: return "Simpson33";
  default:
    nc_assert_always(false&&"invalid IntegrationScheme");
    return "";
  };
}
NC::SABCfg::Cfg NC::SABCfg::createConfig( int sablux )
{
  Cfg c;
  //fixme: revisit all of these
  switch ( sablux ) {
  case 0:
    //Note, this will be the default scheme for VDOSDebye kernels, so it should
    //be "crude but workable":
    c.integScheme = IntegrationScheme::Simpson3;
    c.integSchemeBCSample = IntegrationScheme::Flex5;
    c.egrid_npts = 100;//fixme: too low for sampling speed?
    c.egrid_emin_accuracy = 0.05;
    c.fullCellSamplingARThreshold = 0.15;
    c.bcSamplingLargeSRatioThreshold = 1e-2;
    break;
  case 1:
    c.integScheme = IntegrationScheme::Romberg5;
    c.integSchemeBCSample = IntegrationScheme::Flex5;
    c.egrid_npts = 140;
    c.egrid_emin_accuracy = 0.02;
    c.fullCellSamplingARThreshold = 0.15;
    c.bcSamplingLargeSRatioThreshold = 1e-3;
    break;
  case 2:
    c.integScheme = IntegrationScheme::Flex5;
    c.integSchemeBCSample = IntegrationScheme::Flex5;
    c.egrid_npts = 200;
    c.egrid_emin_accuracy = 0.01;
    c.fullCellSamplingARThreshold = 0.15;
    c.bcSamplingLargeSRatioThreshold = 1e-4;
    break;
  case 3:
    //The default value of VDOS and direct kernels.
    //
    //Default values of Cfg objects are already right for sablux==3, so we just
    //repeat them here to make it easier to get an overview.
    nc_assert( c.integScheme == IntegrationScheme::Flex9);
    nc_assert( c.integSchemeBCSample == IntegrationScheme::Flex9 );
    nc_assert( c.egrid_npts == 300 );
    nc_assert( c.egrid_emin_accuracy == 0.01 );
    nc_assert( c.fullCellSamplingARThreshold == 0.15 );
    nc_assert( c.bcSamplingLargeSRatioThreshold == 1e-5);
    break;
  case 4:
    c.integScheme = IntegrationScheme::Flex17;
    c.integSchemeBCSample = IntegrationScheme::Flex17;
    c.egrid_npts = 450;
    c.egrid_emin_accuracy = 1e-3;
    c.fullCellSamplingARThreshold = 0.15;
    c.bcSamplingLargeSRatioThreshold = 1e-6;
    break;
  case 5:
    c.integScheme = IntegrationScheme::Flex33;
    c.integSchemeBCSample = IntegrationScheme::Flex33;
    c.egrid_npts = 600;
    c.egrid_emin_accuracy = 1e-4;
    c.fullCellSamplingARThreshold = 0.15;
    c.bcSamplingLargeSRatioThreshold = 1e-7;
    break;
  case 6:
    static_assert( IntegrationScheme::MaxPrec
                   == IntegrationScheme::Flex65, "" );
    c.integScheme = IntegrationScheme::MaxPrec;
    c.integSchemeBCSample = IntegrationScheme::MaxPrec;
    c.egrid_npts = 10000;
    c.egrid_emin_accuracy = 1e-7;
    c.fullCellSamplingARThreshold = 0.15;
    c.bcSamplingLargeSRatioThreshold = 1e-8;
    break;
  default:
    NCRYSTAL_THROW2(BadInput,"SABCfg::createConfig sablux="<<sablux
                    <<" outside the valid range (0 to 6)");
  }

  return c;
}

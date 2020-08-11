#ifndef NCrystal_DynInfoUtils_hh
#define NCrystal_DynInfoUtils_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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

#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCSABData.hh"

namespace NCrystal {

  //Access SABData from DynInfo object (or just expand from VDOSDebye
  //model). This might result in some calculations carried out in order to
  //either convert scattering kernels to the correct format or create it by
  //expansion from a VDOS. The vdoslux parameter is used as described in
  //NCMatCfg.hh, and only affects VDOS-based scattering kernels. For the special
  //case of VDOSDebye based kernels, the vdoslux parameter will be reduced by 3
  //(but not less than 0), to prevent wasting resources on what is anyway rather
  //crude input.
  //
  //Unless useCache=false is set, a MT-safe caching mechanism will be employed
  //behind the scene in order to prevent duplication of work in case of repeated
  //calls. The cache can obviously be cleared with the
  //clearSABDataFromDynInfoCaches function (automatically invoked by the global
  //clearCaches function):
  std::shared_ptr<const SABData> extractSABDataFromDynInfo( const DI_ScatKnl*, unsigned vdoslux = 3, bool useCache = true );
  std::shared_ptr<const SABData> extractSABDataFromVDOSDebyeModel( double debyeTemperature,
                                                                   double temperature, SigmaBound boundXS, double elementMassAMU,
                                                                   unsigned vdoslux = 3, bool useCache = true );
  void clearSABDataFromDynInfoCaches();

  //Idealised VDOS based only on Debye temperature:
  VDOSData createVDOSDebye(double debyeTemperature, double temperature, SigmaBound boundXS, double elementMassAMU);
}




#endif


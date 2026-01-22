#ifndef NCrystal_ExtnBC2025_hh
#define NCrystal_ExtnBC2025_hh

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

#include "NCrystal/core/NCDefs.hh"

////////////////////////////////////////////////////////////////////////////////
// Updated recipes evaluating the models described in:                        //
//                                                                            //
// P. J. Becker and P. Coppens Acta Cryst. (1974). A30, 129-147               //
// https://doi.org/10.1107/S0567739474000337                                  //
//                                                                            //
// However, the original 1974 recipes were flawed, so we use updated "BC2025" //
// recipes from:                                                              //
//                                                                            //
// "Revisiting Becker-Coppens (1974): Updated Recipes for Estimating          //
// Extinction Factors in Spherical Crystallites",                             //
// T. Kittelmann et al., 2025 (in preparation)                                //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace Extn {
    namespace BC2025 {

      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      ///                      BC2025 Standard Recipes                       ///
      ///  Precision guarantee for x<1000: Error less than 1e-3*min(y,1-y)   ///
      ///     Reference: T. Kittelmann 2025 (publication in preparation)     ///
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////

      double y_primary( double x, double sintheta );
      double y_scndgauss( double x, double sintheta );
      double y_scndlorentz( double x, double sintheta );
      double y_scndfresnel( double x, double sintheta );

      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////
      ///                       BC2025 Luxury Recipes                        ///
      ///  Precision guarantee for x<1000: Error less than 1e-6*min(y,1-y)   ///
      ///     Reference: T. Kittelmann 2025 (publication in preparation)     ///
      //////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////

      double y_primary_lux( double x, double sintheta );
      double y_scndgauss_lux( double x, double sintheta );
      double y_scndlorentz_lux( double x, double sintheta );
      double y_scndfresnel_lux( double x, double sintheta );

    }
  }
}

#endif

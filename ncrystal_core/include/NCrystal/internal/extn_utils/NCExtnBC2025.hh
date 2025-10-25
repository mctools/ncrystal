#ifndef NCrystal_ExtnBC2025_hh
#define NCrystal_ExtnBC2025_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
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

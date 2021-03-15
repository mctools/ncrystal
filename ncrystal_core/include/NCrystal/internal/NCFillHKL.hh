#ifndef NCrystal_FillHKL_hh
#define NCrystal_FillHKL_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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

#include "NCrystal/NCMatInfo.hh"

namespace NCrystal {

  //Helper function which finds all (h,k,l) planes and calculates their
  //d-spacings and structure factors (fsquared). The first input must be an Info
  //object which already has StructureInfo and AtomInfo (with
  //mean-squared-displacement and atomic coordinates) added. This helper function
  //will then calculate and add the HKL info, so HKL must not yet have been
  //enabled on the passed Info object.
  //
  //Several parameters can be used to fine-tune the behaviour:

  struct FillHKLCfg {

    double dcutoff = 0.5; // Angstrom. Same meaning as in NCMatCfg.hh, but must be
                          // specified as a finite (non-zero) value, since it
                          // affects the hkl range searched.

    double dcutoffup = kInfinity; //Angstrom. Same meaning as in NCMatCfg.hh.

    bool expandhkl = false;// Request that lists of equivalent HKL planes be
                           // created in Info objects.

    double fsquarecut = 1e-5;// Barn. A cutoff value in barn. HKL reflections
                             // with contribution below this will be skipped
                             // (used to skip weak and impossible
                             // reflections). NB: The value of 1e-5 is also used
                             // (hardcoded) in the .nxs factory.

    double merge_tolerance = 1e-6;// Relative tolerance for Fsquare & dspacing
                                  // comparisons when composing hkl families.

    //For specialised expert usage only, all Debye Waller factors can be forced
    //to be unity. If not set, the default is false unless overriden by an
    //environment variables:
    Optional<bool> use_unit_debye_waller_factor = NullOpt;
  };

  void fillHKL( MatInfo &info, FillHKLCfg = {} );

}

#endif

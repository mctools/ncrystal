#ifndef NCrystal_FillHKL_hh
#define NCrystal_FillHKL_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2022 NCrystal developers                                   //
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

namespace NCrystal {

  // The helper function calculateHKLPlanes finds all (h,k,l) planes and
  // calculates their d-spacings and structure factors (fsquared), based on unit
  // cell info and a few configuration parameters (most importantly the dspacing
  // cutoff value). The function needs StructureInfo, but only for the unit cell
  // parameters (the space group number is not used). It also needs the list of
  // atoms in the unit cell, which is required to contain
  // mean-squared-displacement information so Debye-Waller factors can be
  // calculated internally.
  //
  // The parameters which can be used to tune the behaviour are:

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

  HKLList calculateHKLPlanes( const StructureInfo&,
                              const AtomInfoList&,
                              FillHKLCfg = {} );

}

#endif

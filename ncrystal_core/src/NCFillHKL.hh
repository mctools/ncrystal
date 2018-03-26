#ifndef NCrystal_FillHKL_hh
#define NCrystal_FillHKL_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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
#include <limits>

namespace NCrystal {

  //Helper function which finds all (h,k,l) planes and calculates their
  //d-spacings and structure factors (fsquared). The first input must be an Info
  //object which already has StructureInfo and AtomInfo (with
  //mean-squared-displacement and wyckoff positions) added. This helper function
  //will then calculate and add the HKL info, so HKL must not yet have been
  //enabled on the passed Info object.
  //
  //Several parameters can be used to fine-tune the behaviour:
  //
  // dcutoff   : Same meaning as in NCMatCfg.hh, but must be specified as a finite
  //             (non-zero) value, since it affects the hkl range searched.
  // dcutoffup : Same meaning as in NCMatCfg.hh.
  // expandhkl : Request that lists of equivalent HKL planes be created in Info
  //             objects.
  // fsquarecut : A cutoff value in barn. HKL reflections with contribution
  //              below this will be skipped (used to skip weak and impossible
  //              reflections).
  // merge_tolerance : Tolerance for Fsquare & dspacing comparisons when
  //                   composing hkl families.
  void fillHKL( Info &info,
                double dcutoff = 0.5,//angstrom
                double dcutoffup = std::numeric_limits<double>::infinity(),//angstrom
                bool expandhkl = false,
                double fsquarecut = 1e-5,//barn
                double merge_tolerance = 1e-6 );

}

#endif

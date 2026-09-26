#ifndef NCrystal_FillHKL_hh
#define NCrystal_FillHKL_hh

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

#include "NCrystal/interfaces/NCInfo.hh"

namespace NCRYSTAL_NAMESPACE {

  // The helper function calculateHKLPlanes finds all (h,k,l) planes and
  // calculates their d-spacings and structure factors (fsquared), based on unit
  // cell info and a few configuration parameters (most importantly the dspacing
  // cutoff value). The function always needs StructureInfo for the unit cell
  // parameters, including atom positions. If the space group number is also
  // available, the constructed HKL groups will be exactly the
  // symmetry-equivalent groups (i.e. of type HKLInfoType::SymEqvGroup).
  // Otherwise the entries will be grouped by calculated (dspacing,fsquared)
  // values (i.e. of type HKLInfoType::ExplicitHKLs). The list of atoms in the
  // unit cell is required to contain mean-squared-displacement information so
  // Debye-Waller factors can be calculated internally.
  //
  // The parameters which can be used to tune the behaviour are:

  struct FillHKLCfg {

    double dcutoff = 0.5; // Angstrom. Same meaning as in NCMatCfg.hh, but must be
                          // specified as a finite (non-zero) value, since it
                          // affects the hkl range searched.

    double dcutoffup = kInfinity; //Angstrom. Same meaning as in NCMatCfg.hh.

    double fsquarecut = 1e-5;// Barn. A cutoff value in barn. HKL reflections
                             // with contribution below this will be skipped
                             // (used to skip weak and impossible
                             // reflections). NB: The value of 1e-5 is also used
                             // (hardcoded) in the .nxs factory.

    //For usage only when a spacegroup number is NOT available:

    double merge_tolerance = 1e-6;// Relative tolerance for Fsquare & dspacing
                                  // comparisons when composing hkl families
                                  // (this is only used when a space group
                                  // number is NOT available).
    // Why 1e-6: Rounding noise between symmetry-equivalent planes reaches
    // ~4e-11 (relative) in F2 for weak reflections near fsquarecut (d agrees
    // to ~1e-16), so 1e-12 splits true families in 25 stdlib materials.
    // Lower values gain almost nothing (1e-10 separates only ~30 accidental
    // near-coincidences in the whole stdlib, with no CPU or memory effect),
    // and families anyway get the average values of their members.

    std::size_t max_buffered_points = 0;// Max hkl points buffered before they
                                        // are merged into families (0 means
                                        // 16MB worth). Only used when a space
                                        // group number is NOT available, and
                                        // mostly intended for testing.

    //For specialised expert usage only, all Debye Waller factors can be forced
    //to be unity. If not set, the default is false unless overridden by an
    //environment variables:
    Optional<bool> use_unit_debye_waller_factor = NullOpt;
  };

  HKLList calculateHKLPlanes( const StructureInfo&,
                              const AtomInfoList&,
                              FillHKLCfg = {} );

}

#endif

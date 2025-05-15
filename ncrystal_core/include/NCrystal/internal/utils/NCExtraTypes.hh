#ifndef NCrystal_ExtraTypes_hh
#define NCrystal_ExtraTypes_hh

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

#include "NCrystal/core/NCTypes.hh"

///////////////////////////////////////////////////
//                                               //
//  Extra types not intended for the public API  //
//                                               //
///////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  using VectDFM = std::vector<PairDD>;

  struct PreparedPowderInputData : private MoveOnly {

    //Vector of (d-spacing,fsquared*multiplicity) pairs, sorted so larger
    //d-spacing values come first:
    VectDFM d_fm_list;

    //V0 is the unit cell volume in Aa^3 and n_atoms is the number of atoms
    //per unit cell. We only need their multiplied value:
    double v0_times_natoms;

  };

}

#endif

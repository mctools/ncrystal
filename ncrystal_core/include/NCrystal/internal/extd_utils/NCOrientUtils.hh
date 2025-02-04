#ifndef NCrystal_OrientUtils_hh
#define NCrystal_OrientUtils_hh

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

#include "NCrystal/internal/utils/NCRotMatrix.hh"
#include "NCrystal/interfaces/NCSCOrientation.hh"

namespace NCRYSTAL_NAMESPACE {

  class Info;
  struct StructureInfo;

  //Helper function which finally constructs the rotation between the crystal
  //frame and the laboratory frame based on a completed SCOrientation
  //object. Additionally, a reciprocal lattice rotation matrix is needed
  //(typically from getReciprocalLatticeRot(..)). The reciprocal lattice matrix
  //is needed for the case where one or more directions are defined in hkl (RL)
  //coordinates (if is is known that both crystal directions are defined in
  //direct space, it is safe to pass a dummy RotMatrix as the reci_lattice
  //parameter):

  RotMatrix getCrystal2LabRot( const SCOrientation&, const RotMatrix& reci_lattice );

  //Convenience wrapper, providing the reciprocal lattice matrix needed above
  //directly from aStructureInfo object:
  RotMatrix getReciprocalLatticeRot( const StructureInfo& );

}

#endif

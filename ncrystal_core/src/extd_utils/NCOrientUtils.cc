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

#include "NCrystal/internal/extd_utils/NCOrientUtils.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/utils/NCLatticeUtils.hh"
#include <iomanip>

namespace NC = NCrystal;

NC::RotMatrix NC::getCrystal2LabRot( const NC::SCOrientation& sco,
                                     const NC::RotMatrix& reci_lattice )
{
  auto data = sco.getData();
  return verifyLatticeOrientDefAndConstructCrystalRotation( data.dir1,
                                                            data.dir2,
                                                            data.dirtol,
                                                            reci_lattice );
}

NC::RotMatrix NC::getReciprocalLatticeRot( const NC::StructureInfo& si )
{
  return getReciprocalLatticeRot( si.lattice_a, si.lattice_b, si.lattice_c,
                                  si.alpha*kDeg, si.beta*kDeg, si.gamma*kDeg );
}

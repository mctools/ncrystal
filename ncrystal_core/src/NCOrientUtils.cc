////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2018 NCrystal developers                                   //
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

#include "NCOrientUtils.hh"
#include "NCrystal/NCInfo.hh"
#include "NCLatticeUtils.hh"
#include <iomanip>

NCrystal::RotMatrix NCrystal::getCrystal2LabRot( const NCrystal::SCOrientation& sco,
                                                 const NCrystal::RotMatrix& reci_lattice )
{
  if (!sco.isComplete())
    NCRYSTAL_THROW(BadInput,"Incomplete SCOrientation object - must set both primary and secondary directions.");

  Vector dirc[2];
  for (size_t i=0; i < 2; ++i) {
    if (sco.getCrysIsHKL(i)) {
      dirc[i] = reci_lattice * asVect(sco.getCrysDir(i));
    } else {
      dirc[i] = asVect(sco.getCrysDir(i));
    }
  }

  Vector dirl[2];
  dirl[0] = asVect(sco.getLabDir(0));
  dirl[1] = asVect(sco.getLabDir(1));

  if (dirc[0].isParallel(dirc[1],1.0e-6))
    NCRYSTAL_THROW(BadInput,"Chosen SCOrientation directions in the crystal reference frame are too parallel.");

  if (dirl[0].isParallel(dirl[1],1.0e-6))
    NCRYSTAL_THROW(BadInput,"Chosen SCOrientation directions in the laboratory frame are too parallel.");

  const double anglec = dirc[0].angle(dirc[1]);
  const double anglel = dirl[0].angle(dirl[1]);
  if ( ncabs(anglec-anglel)>sco.getTolerance() ) {
    NCRYSTAL_THROW2(BadInput,"Chosen SCOrientation directions in the lab frame are "<<std::setprecision(8)
                    <<anglel*kToDeg<<" deg apart, while the chosen directions in the crystal frame"
                    " are "<<anglec*kToDeg<<" deg apart. This is not within the specified"
                    " tolerance of "<<sco.getTolerance()<<" rad. = "<<sco.getTolerance()*kToDeg<<" deg.");
  }
  //We are within the tolerance, but now ensure exact anglec==anglel by removing
  //components of secondary direction parallel to the primary direction:

  for (size_t i=0; i < 2; ++i) {
    dirc[i].normalise();
    dirl[i].normalise();
  }
  dirc[1] -= dirc[0] * dirc[1].dot(dirc[0]);
  dirl[1] -= dirl[0] * dirl[1].dot(dirl[0]);
  dirc[1].normalise();
  dirl[1].normalise();

  return RotMatrix(dirl[0],dirc[0],dirl[1],dirc[1]);
}

NCrystal::RotMatrix NCrystal::getReciprocalLatticeRot( const NCrystal::Info& cinfo )
{
  if (!cinfo.hasStructureInfo())
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks Structure information.");

  const StructureInfo & si = cinfo.getStructureInfo();
  return getReciprocalLatticeRot( si.lattice_a, si.lattice_b, si.lattice_c,
                                  si.alpha*kDeg, si.beta*kDeg, si.gamma*kDeg );

}

#ifndef NCrystal_PlaneProvider_hh
#define NCrystal_PlaneProvider_hh

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

#include "NCrystal/internal/NCVector.hh"

namespace NCrystal {

  class MatInfo;

  class PlaneProvider {
  public:

    //Interface class providing loops over crystal (demi) planes, each loop
    //giving dspacing, F-squared and a demi normal.

    PlaneProvider();
    virtual ~PlaneProvider();

    //Method used to loop over planes:
    virtual bool getNextPlane(double& dspacing, double& fsq, Vector& demi_normal) = 0;

    //Rewind the looping to prepare for a new loop with getNextPlane (does not
    //need to be called for the first loop):
    virtual void prepareLoop() = 0;

    //Whether or not it is safe to call getNextPlane and prepareLoop. If calling
    //anyway, exceptions will be thrown by conforming implementations. A false
    //return value here usually indicates incomplete information for normals to
    //be provided:
    virtual bool canProvide() const = 0;
  };

  //Creates standard plane provider from Info object, which will attempt various
  //means of producing the HKL normals (preferring actual deminormals if
  //available, then expanded hkl info and finally falling back to attempting
  //their construction based on space group and multiplicity info):

  //Version partaking in lifetime management of MatInfo:
  std::unique_ptr<PlaneProvider> createStdPlaneProvider(shared_obj<const MatInfo>);

  //Version in which caller guarantees MatInfo object will remain alive as long
  //as any plane provider methods are called:
  std::unique_ptr<PlaneProvider> createStdPlaneProvider(const MatInfo*);

}


#endif

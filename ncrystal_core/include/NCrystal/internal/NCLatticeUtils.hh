#ifndef NCrystal_LatticeUtils_hh
#define NCrystal_LatticeUtils_hh

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

#include "NCrystal/internal/NCRotMatrix.hh"

namespace NCrystal {

  //Construct cell rotation matrix from lattice parameters (in angstrom and
  //radians):
  RotMatrix getLatticeRot( double lattice_a, double lattice_b, double lattice_c,
                           double alpha, double beta, double gamma );

  //Reciprocal lattice rotation (inverted lattice transform, multiplied by 2pi):
  RotMatrix getReciprocalLatticeRot( double lattice_a, double lattice_b, double lattice_c,
                                     double alpha, double beta, double gamma );

  //Given a set of lattice parameters, translate dcutoff into required maximum
  //values of |h|, |k| and |l|:
  struct MaxHKL{ int h, k, l; };
  MaxHKL estimateHKLRange( double dcutoff,
                           double lattice_a, double lattice_b, double lattice_c,
                           double alpha, double beta, double gamma );

  //Validate that lattice lengths are compatible with given spacegroup. For
  //space groups where a==b or a==c, it is allowed to provide b=0 or c=0, and
  //the function will then update the values of b and/or c accordingly. Any
  //inconsistencies or errors will result in BadInput exceptions thrown.
  void checkAndCompleteLattice( unsigned spacegroup, double a, double& b, double & c );

  //Calculate d-spacing from Miller index and reciprocal lattice rotation:
  double dspacingFromHKL( int h, int k, int l, const RotMatrix& rec_lat );


}

#endif

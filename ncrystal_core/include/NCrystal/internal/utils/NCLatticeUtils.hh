#ifndef NCrystal_LatticeUtils_hh
#define NCrystal_LatticeUtils_hh

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

namespace NCRYSTAL_NAMESPACE {

  //Construct cell rotation matrix from lattice parameters (in angstrom and
  //radians):
  RotMatrix getLatticeRot( double lattice_a, double lattice_b, double lattice_c,
                           double alpha, double beta, double gamma );

  // Validate orientation parameters (check against null-vectors, parallel
  // vectors, consistency of angles in the two frames given dirtol,
  // etc.). Throws BadInput exception in case of issues.
  void precheckLatticeOrientTol( double dirtol );
  void precheckLatticeOrientDir( const OrientDir& );
  void verifyLatticeOrientDef( const LabAxis& l1, const CrystalAxis& c1,
                               const LabAxis& l2, const CrystalAxis& c2,
                               double dirtol );

  //Same, but on OrientDir's which might include HKLPoint's, in which case the
  //validation won't be complete (due to the absence of a lattice rotation
  //matrix):
  void precheckLatticeOrientDef( const OrientDir& dir1, const OrientDir& dir2, double dirtol );

  //Fully validate orientation parameters and construct crystal rotation matrix
  //(transforms lab frame -> crystal frame). Requires a reciprocal lattice
  //rotation matrix.
  RotMatrix verifyLatticeOrientDefAndConstructCrystalRotation( const OrientDir& dir1,
                                                               const OrientDir& dir2,
                                                               double dirtol,
                                                               const RotMatrix& reci_lattice );


  //Reciprocal lattice rotation (inverted lattice transform, multiplied by 2pi):
  RotMatrix getReciprocalLatticeRot( double lattice_a, double lattice_b, double lattice_c,
                                     double alpha, double beta, double gamma );

  //Given a set of lattice parameters, translate dcutoff into required maximum
  //values of |h|, |k| and |l|:
  struct MaxHKL{ int h, k, l; };
  MaxHKL estimateHKLRange( double dcutoff,
                           double lattice_a, double lattice_b, double lattice_c,
                           double alpha, double beta, double gamma );

  //Calculate d-spacing from Miller index and reciprocal lattice rotation:
  double dspacingFromHKL( int h, int k, int l, const RotMatrix& rec_lat );

  //Validate that lattice lengths are compatible with given spacegroup. For
  //space groups where a==b or a==c, it is allowed to provide b=0 or c=0, and
  //the function will then update the values of b and/or c accordingly. Any
  //inconsistencies or errors will result in BadInput exceptions thrown. It is
  //allowed to provide spacegroup==0, in which case the function won't do much:
  void checkAndCompleteLattice( unsigned spacegroup, double a, double& b, double& c );

  //Same for angles, except that all three angles might be left zero in case the
  //space group defines them:
  void checkAndCompleteLatticeAngles( unsigned spacegroup, double& alpha, double& beta, double& gamma );

  //Spacegroup utils:
  enum CrystalSystem {
    Triclinic,
    Monoclinic,
    Orthorhombic,
    Tetragonal,
    Trigonal,
    Hexagonal,
    Cubic
  };
  CrystalSystem crystalSystem( int spacegroup );

  //Try to estimate the HKL point based on normal + dspacing, essentially by hkl
  //= lattice_rot * normal / dspacing, but then returning the smallest of
  //(hkl,-hkl).  Will throw CalcError in case the normal does not actually
  //correspond to one coming from an integral hkl point.
  std::tuple<int,int,int> normalAndDSpacingToHKLIndex( const RotMatrix& lattice_rot,
                                                       double dspacing,
                                                       const Vector& normal );


}

#endif

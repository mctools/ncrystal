////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2019 NCrystal developers                                   //
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

#include "NCrystal/NCSCOrientation.hh"
#include "NCVector.hh"

namespace NCrystal {

  SCOrientation::SCOrientation()
    : m_tolerance(-1.0)
  {
    for (int i = 0; i < 2; ++i) {
      m_crystal_is_hkl[i] = false;
      for (int xyz = 0; xyz < 3; ++xyz) {
        m_crystal[i][xyz] = m_lab[i][xyz] = 0.0;
      }
    }
  }

  SCOrientation::~SCOrientation()
  {
  }

  void SCOrientation::setPrimaryDirection( const double (&crystal_direction)[3],
                                           const double (&direction)[3])
  {
    if ( ! asVect(direction).mag2() )
      NCRYSTAL_THROW(BadInput,"Specified lab-direction is a null-vector.");

    m_crystal_is_hkl[0] = false;
    for (int xyz = 0; xyz < 3; ++xyz) {
      m_crystal[0][xyz] = crystal_direction[xyz];
      m_lab[0][xyz] = direction[xyz];
    }

    if (isComplete())
      checkInput();
  }

  void SCOrientation::setSecondaryDirection( const double (&crystal_direction)[3],
                                             const double (&direction)[3],
                                             double tolerance )
  {
    if ( ! asVect(direction).mag2() )
      NCRYSTAL_THROW(BadInput,"Specified lab-direction is a null-vector.");

    m_crystal_is_hkl[1] = false;
    for (int xyz = 0; xyz < 3; ++xyz) {
      m_crystal[1][xyz] = crystal_direction[xyz];
      m_lab[1][xyz] = direction[xyz];
    }

    m_tolerance = tolerance;

    if (isComplete())
      checkInput();
  }

  void SCOrientation::setPrimaryDirection( double h, double k, double l,
                                           const double (&direction)[3] )
  {
    if ( ! asVect(direction).mag2() )
      NCRYSTAL_THROW(BadInput,"Specified lab-direction is a null-vector.");

    m_crystal_is_hkl[0] = true;
    m_crystal[0][0] = h;
    m_crystal[0][1] = k;
    m_crystal[0][2] = l;
    for (int xyz = 0; xyz < 3; ++xyz)
      m_lab[0][xyz] = direction[xyz];

    if (isComplete())
      checkInput();
  }

  void SCOrientation::setSecondaryDirection( double h, double k, double l,
                                             const double (&direction)[3],
                                             double tolerance )
  {
    if ( ! asVect(direction).mag2() )
      NCRYSTAL_THROW(BadInput,"Specified lab-direction is a null-vector.");

    m_crystal_is_hkl[1] = true;
    m_crystal[1][0] = h;
    m_crystal[1][1] = k;
    m_crystal[1][2] = l;
    for (int xyz = 0; xyz < 3; ++xyz)
      m_lab[1][xyz] = direction[xyz];

    m_tolerance = tolerance;

    if (isComplete())
      checkInput();
  }

  bool SCOrientation::isComplete() const
  {
    return asVect(m_lab[0]).mag2() && asVect(m_lab[1]).mag2();
  }

  SCOrientation::SCOrientation(const SCOrientation& o)
  {
    *this = o;
  }

  SCOrientation& SCOrientation::operator=(const SCOrientation& o)
  {
    m_tolerance = o.m_tolerance;
    for (int i = 0; i < 2; ++i) {
      m_crystal_is_hkl[i] = o.m_crystal_is_hkl[i];
      for (int xyz = 0; xyz < 3; ++xyz) {
        m_crystal[i][xyz] = o.m_crystal[i][xyz];
        m_lab[i][xyz] = o.m_lab[i][xyz];
      }
    }
    return *this;
  }

  void SCOrientation::checkInput() const
  {
    nc_assert_always(isComplete());

    //Check that user provided sensible input: Tolerance in (0,pi] and
    //non-parallel and non-null vectors. Note that we can not check mixed-mode
    //input (one direction as normal of hkl plane and one directly in crystal
    //coordinates) here, and neither can we check the tolerances in case of
    //hkl-modes, so it is not guaranteed that there are no remaining problems
    //after the checks here.

    for (int i = 0; i < 2; ++i) {
      if ( ! asVect(m_crystal[i]).mag2() ) {
        NCRYSTAL_THROW(BadInput, m_crystal_is_hkl[i]
                       ? "Specified point in hkl space is a null-vector"
                       : "Specified direction in crystal frame is a null-vector");
      }
      nc_assert_always( asVect(m_lab[i]).mag2() );//should have been tested already
    }

    if ( m_tolerance <= 0.0 || m_tolerance > kPi )
      NCRYSTAL_THROW(BadInput, "Tolerance must be in interval (0.0,pi]");

    if ( asVect(m_lab[0]).isParallel( asVect(m_lab[1]), 1.0e-6 ) )
      NCRYSTAL_THROW(BadInput, "Specified primary and secondary lab directions are parallel");

    if ( m_crystal_is_hkl[0] == m_crystal_is_hkl[1] ) {
      if ( asVect(m_crystal[0]).isParallel( asVect(m_crystal[1]), 1.0e-6 ) ) {
        NCRYSTAL_THROW(BadInput, m_crystal_is_hkl[0]
                       ? "Specified primary and secondary hkl points have planes with parallel normals"
                       : "Specified primary and secondary directions in the crystal frame are parallel" );
      }
    }
  }
}

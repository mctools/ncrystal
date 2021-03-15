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

#include "NCrystal/NCSCOrientation.hh"
#include "NCrystal/internal/NCVector.hh"

namespace NC = NCrystal;

struct NC::SCOrientation::Impl {

  static const char * crystalAxisName( const CrystalAxis& ) { return "CrystalAxis"; }
  static const char * crystalAxisName( const HKLPoint& ) { return "HKLPoint"; }

  template <class T>
  static inline void setDir(SCOrientation* THIS, unsigned idx, const T& c, const LabAxis& l, Optional<double> opt_tol )
  {
    static_assert(std::is_same<T,HKLPoint>::value || std::is_same<T,CrystalAxis>::value, "");
#if __cplusplus >= 201703L
    static_assert(!std::is_aggregate<SCOrientation>::value);
#endif
    nc_assert( idx < 2 );
    const char * psname = idx == 0 ? "primary" : "secondary";
    if (!c.template as<Vector>().mag2())
      NCRYSTAL_THROW2(BadInput,"Specified "<<psname<<" "<<crystalAxisName(c)<<" is a null-vector.");
    if (!l.as<Vector>().mag2())
      NCRYSTAL_THROW2(BadInput,"Specified "<<psname<<" LabAxis is a null-vector.");

    if ( opt_tol.has_value() ) {
      nc_assert( idx == 1 );
      double tol = opt_tol.value();
      if ( !( tol > 0.0 ) || tol > kPi )
        NCRYSTAL_THROW(BadInput, "Tolerance must be in interval (0.0,pi]");
      THIS->m_tolerance = opt_tol.value();
    }

    THIS->m_crystal[idx] = c;
    if (!std::is_same<T,HKLPoint>::value)
      THIS->m_crystal[idx].get<CrystalAxis>().as<Vector>().normalise();

    THIS->m_lab[idx] = l;
    THIS->m_lab[idx].value().as<Vector>().normalise();

    if (THIS->isComplete()) {

      //Check that user did not provide parallel vectors as input.

      //Note that we can not check mixed-mode input (one direction as normal of hkl
      //plane and one directly in crystal coordinates) here, and neither can we
      //check the tolerances in case of hkl-modes, so it is not guaranteed that
      //there are no remaining problems after the checks here.

      constexpr double paralleltol = 1e-6;
      if ( THIS->m_lab[0].value().as<Vector>().isParallel( THIS->m_lab[1].value().as<Vector>(), paralleltol ) )
        NCRYSTAL_THROW(BadInput, "Specified primary and secondary lab directions are parallel");

      if ( THIS->m_crystal[0].has_value<CrystalAxis>() && THIS->m_crystal[1].has_value<CrystalAxis>() ) {
        //Both CrystalAxis:
        if ( THIS->m_crystal[0].get<CrystalAxis>().as<Vector>().isParallel( THIS->m_crystal[1].get<CrystalAxis>().as<Vector>(), paralleltol ) )
          NCRYSTAL_THROW(BadInput, "Specified primary and secondary directions in the crystal frame are parallel" );
      } else if ( THIS->m_crystal[0].has_value<HKLPoint>() && THIS->m_crystal[1].has_value<HKLPoint>() ) {
        //Both HKLPoint:
        if ( THIS->m_crystal[0].get<HKLPoint>().as<Vector>().isParallel( THIS->m_crystal[1].get<HKLPoint>().as<Vector>(), paralleltol ) )
          NCRYSTAL_THROW(BadInput, "Specified primary and secondary hkl points have planes with parallel normals" );
      } else {
        //mixed-mode, can't check.
        nc_assert( !THIS->m_crystal[0].empty() && !THIS->m_crystal[1].empty() );
      }
    }
  }
};

void NC::SCOrientation::setPrimaryDirection( const CrystalAxis& c, const LabAxis& l ) { Impl::setDir(this,0,c,l,NullOpt); }
void NC::SCOrientation::setSecondaryDirection( const CrystalAxis& c, const LabAxis& l, double t ) { Impl::setDir(this,1,c,l,t); }
void NC::SCOrientation::setPrimaryDirection( const HKLPoint& c, const LabAxis& l ) { Impl::setDir(this,0,c,l,NullOpt); }
void NC::SCOrientation::setSecondaryDirection( const HKLPoint& c, const LabAxis& l, double t ) { Impl::setDir(this,1,c,l,t); }

bool NC::SCOrientation::isComplete() const
{
  nc_assert(m_lab[0].has_value() == !m_crystal[0].empty());
  nc_assert(m_lab[1].has_value() == !m_crystal[1].empty());
  return m_lab[0].has_value() && m_lab[1].has_value();
}

void NC::SCOrientation::clear()
{
  m_crystal[0].clear();
  m_crystal[1].clear();
  m_lab[0].reset();
  m_lab[1].reset();
  m_tolerance = default_tolerance;
}

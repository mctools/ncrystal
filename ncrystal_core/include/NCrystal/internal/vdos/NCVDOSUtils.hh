#ifndef NCrystal_VDOSUtils_hh
#define NCrystal_VDOSUtils_hh

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

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/internal/utils/NCRect.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace VDOS {

    ///////////////////////////////////////////////////
    // Various utilities needed for VDOS processing. //
    ///////////////////////////////////////////////////

    //Interval where f(x) = x^n*exp(-x) is above eps*fpeak.
    PairDD rangeXNexpMX(unsigned n, double eps, double accuracy = 1e-13 );

    // Returns the intersection between the provided Rectangle in the alpha-beta
    // plane, and the kinematically available phasespace for a neutron of a
    // given E/kT (i.e. the set of points satisfying (alpha-beta)^2 <
    // 4*E_div_kT*alpha). The returned rectangle is empty when the intersection
    // has no area. Boundary-only intersections are therefore excluded, but of
    // course the usual caveats of floating point arithmetic applies.
    Rectangle findABExtentWithinKB( const Rectangle&, double E_div_kT );
  }
}

#endif

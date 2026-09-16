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

#include "NCrystal/internal/utils/NCRect.hh"
#include "NCrystal/internal/utils/NCSpan.hh"

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

    // In-place spacing out of the positive and negative parts of a finite,
    // strictly increasing grid. For adjacent nonzero points a < b with the same
    // sign, the points are considered too close when (0.0<=rtol<1.0):
    //
    //   b - a < rtol * min(abs(a), abs(b)).
    //
    // The first and last point of both positive and negative subgrids remain
    // fixed, while interior points are moved in place as needed. If there is
    // not enough room between fixed points to space out intermediate points, an
    // error is raised. In case of any errors, the grid might be left in an
    // invalid state.
    void spaceOutGrid(Span<double>, double rtol);

  }
}

#endif

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

namespace NCRYSTAL_NAMESPACE {

  namespace VDOS {

    ///////////////////////////////////////////////////
    // Various utilities needed for VDOS processing. //
    ///////////////////////////////////////////////////

    //Interval where f(x) = x^n*exp(-x) is above eps*fpeak.
    PairDD rangeXNexpMX(unsigned n, double eps, double accuracy = 1e-13 );

    //Find the extreme (as in highest alpha, lowest beta) kinematically
    //accessible point in the provided rectangular region in (alpha,beta) space
    //for a neutron with a given energy/kT. Returns NullOpt in case no point is
    //accessible (or the area only touches the phasespace). Note that we on
    //purpose consider only the kinematic edge given by the alpha+(beta) and
    //beta=-E/kT curves, ignoring the alpha-(beta) curve.
    //
    //NB: with e=E/kT with have:
    //  alpha+-(beta) = 2*e+beta +- 2*sqrt(e*(e+beta))
    Optional<PairDD>
    findExtremeSABPointWithinAlphaPlusCurve( double E_div_kT,
                                             PairDD alphaRange,
                                             PairDD betaRange );

    //Same as findExtremeSABPointWithinAlphaPlusCurve, but only testing whether
    //or not a given single point is accessible, in the sense that it has
    //beta>=-E/kT and alpha<alpha+(beta). Note that this deliberately ignores
    //the alpha-(beta) curve.:
    bool sabPointWithinAlphaPlusCurve( double E_div_kT,
                                       double alpha,
                                       double beta );
  }
}

////////////////////////////
// Inline implementations //
////////////////////////////

inline bool NCrystal::VDOS::sabPointWithinAlphaPlusCurve( double E_div_kT,
                                                          double alpha,
                                                          double beta )
{
  nc_assert( alpha >= 0.0 );
  nc_assert( E_div_kT > 0.0 );
  const double e = E_div_kT;
  const double epb = e + beta;
  if ( epb < 0.0 )
    return false;
  //With e=E/kT, "within" means:
  //  alpha+(beta)>=alpha <=> sqrt(e*(e+beta))>=(alpha-beta)/2-e
  const double t = 0.5 * ( alpha - beta ) - e;
  return t <= 0.0 || e*epb >= t*t;
}

#endif

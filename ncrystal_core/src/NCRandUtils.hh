#ifndef NCrystal_RandUtils_hh
#define NCrystal_RandUtils_hh

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

#include "NCrystal/NCDefs.hh"
#include <vector>

namespace NCrystal {

  //Todo for NC2: create Vector interfaces for interfaces below (can be inlined wrappers):
  double randIsotropicScatterAngle( RandomBase * );
  double randIsotropicScatterMu( RandomBase * );
  void randIsotropicDirection( RandomBase *, double (&)[3]);//result will be unit vector
  void randDirectionGivenScatterMu( RandomBase *, double mu/*=cos(scatangle)*/, const double(&in)[3], double(&out)[3]);//outdir will be unit vector
  void randPointOnUnitCircle( RandomBase *,  double & x, double& y );//Sample a random point on the unit circle
  double randNorm( RandomBase * );//sample single value from unit Gaussian
  void randNorm( RandomBase *, double&g1, double&g2);//sample two independent values from unit Gaussian.
  std::size_t pickRandIdxByWeight( RandomBase *, const std::vector<double>& commulvals);//pick index according to weights (values must be commulative)

}


////////////////////////////
// Inline implementations //
////////////////////////////

inline double NCrystal::randIsotropicScatterAngle( NCrystal::RandomBase * rand )
{
  //acos() is not fast, but hard to come up with something faster. We could
  //consider downgrading to single-precision results, and using acosf instead of
  //acos (~twice as fast).
  const double mu = -1.0+rand->generate()*2.0;
  nc_assert(mu>=-1.0&&mu<=1.0);
  return std::acos(mu);
}

inline double NCrystal::randIsotropicScatterMu( NCrystal::RandomBase * rand )
{
  //... or, just use randIsotropicScatterMu instead
  const double mu = -1.0+rand->generate()*2.0;
  nc_assert(mu>=-1.0&&mu<=1.0);
  return mu;
}

#endif

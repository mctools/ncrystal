////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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

#include "NCrystal/NCBkgdExtCurve.hh"
#include "NCrystal/NCInfo.hh"
#include "NCRandUtils.hh"

namespace NC = NCrystal;

NC::BkgdExtCurve::BkgdExtCurve( const Info* ci )
  : ScatterIsotropic("BkgdExtCurve"),
    m_ci(ci)
{
  nc_assert_always(ci);
  if (!ci->providesNonBraggXSects())
    NCRYSTAL_THROW(MissingInfo,"BkgdExtCurve: Passed Info object lacks NonBraggXSects needed for cross sections.");
  validate();
}

NC::BkgdExtCurve::~BkgdExtCurve() = default;

double NC::BkgdExtCurve::crossSectionNonOriented(double ekin) const
{
  return m_ci.obj()->xsectScatNonBragg(ekin2wl(ekin));
}

void NC::BkgdExtCurve::generateScatteringNonOriented( double, double& angle, double& de ) const
{
  angle = randIsotropicScatterAngle(getRNG());
  de = 0.0;
}

void NC::BkgdExtCurve::generateScattering( double, const double (&)[3],
                                           double (&outdir)[3], double& de ) const
{
  randIsotropicDirection(getRNG(),outdir);
  de = 0.0;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

#include "NCrystal/NCScatterXSCurve.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCException.hh"
#include "NCMath.hh"

NCrystal::ScatterXSCurve::ScatterXSCurve(const Info* ci, const char * calcname, bool thermalise )
  : NonOrientedScatter(calcname),
    m_ci(ci),
    m_tempk(0)
{
  nc_assert_always(ci);
  ci->ref();
  if (thermalise) {
    if (ci->hasTemperature()) {
      //mass dependency cancels for thermal energy distribution, so using neutron mass=1.0
      m_tempk = ci->getTemperature();
    } else {
      //missing information will result in exception if asked to generate scattering:
      m_tempk = -1.0;
    }
  }
  validate();
}

NCrystal::ScatterXSCurve::~ScatterXSCurve()
{
  m_ci->unref();
}

double NCrystal::ScatterXSCurve::calcDeltaE(double ekin) const
{
  if (!m_tempk)
    return 0;//always elastic

  //thermalise immediately, so new energy will be given by maxwell distribution.
  if (m_tempk<0)
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks Temperature information needed"
                   " to generate scatterings when thermalise flag is set.");

  //TODO for NC2: From comment on DGSW-192: "According to XX, the thermalise=true option
  //can generate impossible (p,E non-conserving) scenarios. So we should check
  //that the result is OK and if not throw another (up to N times, then default
  //to de=0).". Have to double-check equations to see if this is true.

  return genThermalNeutronEnergy(m_tempk, this->rand()) - ekin;
}

void NCrystal::ScatterXSCurve::generateScatteringNonOriented( double ekin,
                                                          double& angle_radians, double& delta_ekin_eV ) const
{
  angle_radians = randIsotropicScatterAngle();
  delta_ekin_eV = calcDeltaE(ekin);
}

void NCrystal::ScatterXSCurve::generateScattering( double ekin, const double (&)[3],
                                               double (&resulting_neutron_direction)[3], double& delta_ekin_eV ) const
{
  randIsotropicDirection(resulting_neutron_direction);
  delta_ekin_eV = calcDeltaE(ekin);
}


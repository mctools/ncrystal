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

#include "NCrystal/NCSimpleBkgd.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCException.hh"
#include "NCMath.hh"

NCrystal::SimpleBkgd::SimpleBkgd(const Info* ci, bool thermalise)
  : NonOrientedScatter(thermalise?"SimpleBkgdT":"SimpleBkgdE"),
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

NCrystal::SimpleBkgd::~SimpleBkgd()
{
  m_ci->unref();
}

double NCrystal::SimpleBkgd::crossSectionNonOriented(double ekin) const
{
  if (!m_ci->providesNonBraggXSects())
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks NonBraggXSects needed for cross sections.");
  return m_ci->xsectScatNonBragg(ekin2wl(ekin));
}

double NCrystal::SimpleBkgd::calcDeltaE(double ekin) const
{
  if (!m_tempk)
    return 0;//always elastic

  //thermalise immediately, so new energy will be given by maxwell distribution.
  if (m_tempk<0)
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks Temperature information needed"
                   " to generate scatterings when thermalise flag is set.");

  return genThermalNeutronEnergy(m_tempk, this->rand()) - ekin;
}

void NCrystal::SimpleBkgd::generateScatteringNonOriented( double ekin,
                                                          double& angle_radians, double& delta_ekin_eV ) const
{
  angle_radians = randIsotropicScatterAngle();
  delta_ekin_eV = calcDeltaE(ekin);
}

void NCrystal::SimpleBkgd::generateScattering( double ekin, const double (&)[3],
                                               double (&resulting_neutron_direction)[3], double& delta_ekin_eV ) const
{
  randIsotropicDirection(resulting_neutron_direction);
  delta_ekin_eV = calcDeltaE(ekin);
}


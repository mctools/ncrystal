#ifndef NCrystal_PCBragg_hh
#define NCrystal_PCBragg_hh

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

#include "NCrystal/NCScatterIsotropic.hh"
#include <vector>

namespace NCrystal {

  class Info;

  class PCBragg : public ScatterIsotropic {
  public:

    //Calculates Bragg diffraction in a polycrystalline/powdered material.

    //Constructor:
    PCBragg(const Info*);

    //The cross-section (in barns):
    virtual double crossSectionNonOriented(double ekin) const;

    //There is a maximum wavelength at which Bragg diffraction is possible,
    //so ekin_low will be set to reflect this (ekin_high will be set to infinity):
    virtual void domain(double& ekin_low, double& ekin_high) const { ekin_low = m_threshold_ekin; ekin_high = infinity; }

    //Generate scatter angle according to Bragg diffraction (defaulting to
    //isotropic if the provided wavelength is above threshold()). This is
    //elastic scattering and will always result in delta_ekin_eV=0:
    virtual void generateScatteringNonOriented( double ekin,
                                                double& angle_radians, double& delta_ekin_eV ) const;

    virtual void generateScattering( double ekin, const double (&neutron_direction)[3],
                                     double (&resulting_neutron_direction)[3], double& delta_ekin_eV ) const;

  protected:
    virtual ~PCBragg();
    void genSinThetaBragg(double,double&) const;
    double m_threshold_wl;
    double m_threshold_ekin;
    double m_xsectfact;
    std::vector<double> m_2d;
    std::vector<double> m_fdm_commul;
  };
}

#endif

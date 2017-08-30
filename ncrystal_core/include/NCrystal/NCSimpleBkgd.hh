#ifndef NCrystal_SimpleBkgd_hh
#define NCrystal_SimpleBkgd_hh

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

#include "NCrystal/NCNonOrientedScatter.hh"

namespace NCrystal {

  class Info;

  class SimpleBkgd : public NonOrientedScatter {
  public:

    //Calculates non-Bragg scattering in a crystal, based on the XSectProvider
    //info (using xsectScatNonBragg()) info in the passed Info object.
    //
    //Scatter angles will be isotropic.
    //
    //As accurate energy transfer information is unavailable, it will, depending
    //on the flag passed in the constructor, either be approximated as always 0,
    //or as always completely thermalising immediately according to the
    //temperature of the crystal (requires temperature info availability).

    //Constructor (getName() afterwards returns "SimpleBkgdT" if thermalise==true, else "SimpleBkgdE")
    SimpleBkgd(const Info*, bool thermalise = true );

    virtual double crossSectionNonOriented(double ekin) const;

    virtual void generateScatteringNonOriented( double ekin_wavelength_aangstrom,
                                                double& angle_radians, double& delta_ekin_eV ) const;

    virtual void generateScattering( double ekin, const double (&neutron_direction)[3],
                                     double (&resulting_neutron_direction)[3], double& delta_ekin_eV ) const;

  protected:
    virtual ~SimpleBkgd();
    const Info* m_ci;
    double m_tempk;
    double calcDeltaE(double) const;
  };
}

#endif

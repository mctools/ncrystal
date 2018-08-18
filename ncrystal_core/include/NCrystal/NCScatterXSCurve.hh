#ifndef NCrystal_ScatterXSCurve_hh
#define NCrystal_ScatterXSCurve_hh

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

namespace NCrystal {

  class Info;

  class NCRYSTAL_API ScatterXSCurve : public ScatterIsotropic {
  public:

    //Scatter base class for models which provides only cross-sections and no
    //specific ability to generate scatterings.
    //
    //Scatter angles will be modelled as isotropic. As accurate energy transfer
    //information is unavailable, it will, depending on the value of the
    //thermaise flag passed in the constructor, either be approximated as always
    //0 (elastic), or as always completely thermalising immediately according to
    //the temperature of the crystal (requires temperature info availability).
    //
    //Derived classes should always implement the crossSectionNonOriented(..)
    //method, and possibly also the domain(..) and crossSection(..) methods if
    //appropriate.

    ScatterXSCurve(const Info*, const char * calcname, bool thermalise );


    virtual void generateScatteringNonOriented( double ekin_wavelength,
                                                double& angle, double& delta_ekin ) const;

    virtual void generateScattering( double ekin, const double (&neutron_direction)[3],
                                     double (&resulting_neutron_direction)[3], double& delta_ekin ) const;

  protected:
    virtual ~ScatterXSCurve();
    const Info* m_ci;
    double m_tempk;
    double calcDeltaE(double) const;
  };
}

#endif

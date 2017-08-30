#ifndef NCrystal_Scatter_hh
#define NCrystal_Scatter_hh

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

#include "NCrystal/NCProcess.hh"

/////////////////////////////////////////////////////////////////////
// Base class for calculations of scattering in materials, adding  //
// generateScattering methods in addition to the crossSection      //
// methods from the Process class.                                 //
//                                                                 //
// Note that the unit for kinetic energy in the calls below is eV  //
// (electronvolt) and angles are in radians.                       //
/////////////////////////////////////////////////////////////////////

namespace NCrystal {

  class Scatter : public Process {
  public:

    Scatter(const char * calculator_type_name);

    //Assuming a scattering took place, the following method generate an energy
    //transfer and new direction for the neutron.
    virtual void generateScattering( double ekin, const double (&neutron_direction)[3],
                                     double (&resulting_neutron_direction)[3], double& delta_ekin ) const = 0;

    //For non-oriented scatter calculators, callers can either use the methods
    //above to access cross-sections and scatterings (with free choice of
    //reference frame for the direction vectors), or use the following methods
    //for convenience (note that the returned angle is the angle from incoming
    //to outgoing neutron):
    virtual void generateScatteringNonOriented( double ekin,
                                                double& angle, double& delta_ekin ) const;

  protected:
    virtual ~Scatter();
  };

}

#endif

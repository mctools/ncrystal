#ifndef NCrystal_FreeGas_hh
#define NCrystal_FreeGas_hh

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

#include "NCrystal/NCScatterIsotropic.hh"
#include "NCrystal/NCAtomData.hh"

namespace NCrystal {

  class FreeGas final : public ScatterIsotropic {
  public:

    ////////////////////////////////////////////////////////////////////////////////
    //                                                                            //
    // Thermal scattering with a free-gas model. Specifically this carefully      //
    // implements cross-sections and sampling of scatterings based on the neutron //
    // scattering function, S, given in eq. 19 of the NCrystal sampling paper     //
    // (10.1016/j.jcp.2018.11.043) and is thus a fully quantum mechanical model.  //
    //                                                                            //
    // See NCFreeGasUtils.hh for more details.                                    //
    //                                                                            //
    ////////////////////////////////////////////////////////////////////////////////

    //Explicitly provide parameters:
    enum class SigmaType { FREE, BOUND };
    FreeGas( double temp_kelvin,
             double target_mass_amu,
             double sigma_barn,
             SigmaType sigma_type = SigmaType::FREE );

    //Take parameters from AtomData object:
    FreeGas( double temp_kelvin, const AtomData& );

    double crossSectionNonOriented(double ekin) const override;
    double crossSection(double ekin, const double (&neutron_direction)[3] ) const override;

    void generateScatteringNonOriented( double ekin, double& angle, double& delta_ekin ) const override;

    void generateScattering( double ekin, const double (&neutron_direction)[3],
                             double (&resulting_neutron_direction)[3], double& delta_ekin ) const override;

  protected:
    virtual ~FreeGas();
    struct Impl;
    Pimpl<Impl> m_impl;
  };
}

#endif

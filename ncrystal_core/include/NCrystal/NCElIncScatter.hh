#ifndef NCrystal_ElIncScatter_hh
#define NCrystal_ElIncScatter_hh

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

namespace NCrystal {

  class ElIncXS;
  class Info;

  class ElIncScatter : public ScatterIsotropic {
  public:

    //Model elastic-incoherent scatterings based on the Debye model with
    //isotropic atomic displacements. For each element in the crystalline
    //material, the mean-squared-displacement must be provided, along with the
    //corresponding bound incoherent scattering cross-section and the scale. The
    //scale will often be the fraction of the element (by count). For more
    //details, see section 2.3 of https://doi.org/10.1016/j.cpc.2019.07.015 and
    //comments in the NCElIncXS.hh header file.

    //Construct from Info:
    ElIncScatter( const Info* );

    //Constructor similar to the ElIncXS constructor:
    ElIncScatter( const VectD& elements_meanSqDisp,
                  const VectD& elements_boundincohxs,
                  const VectD& elements_scale );

    double crossSectionNonOriented(double ekin) const override;
    void generateScatteringNonOriented( double ekin, double& angle, double& delta_ekin ) const override;
    void generateScattering( double ekin, const double (&neutron_direction)[3],
                             double (&resulting_neutron_direction)[3], double& delta_ekin ) const override;

  protected:
    virtual ~ElIncScatter();
    std::unique_ptr<ElIncXS> m_elincxs;
  };
}

#endif

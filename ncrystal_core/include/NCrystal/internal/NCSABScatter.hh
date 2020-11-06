#ifndef NCrystal_SABScatter_hh
#define NCrystal_SABScatter_hh

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

  class DI_ScatKnl;
  class SABData;
  namespace SAB {
    class SABScatterHelper;
  }

  class SABScatter : public ScatterIsotropic {
  public:

    //Provides cross-sections and samplings based on an S(alpha,beta) scattering
    //kernel.

    //Construct from SABData and (optionally) an energy grid. The first
    //constructor takes a DI_ScatKnl object from the DynamicInfo list on an an
    //Info object, along with (optionally) the unique id of the Info object. If
    //an energy grid is not supplied, a reasonable default value will be
    //used. When possible, caching will be enabled by default - making sure that
    //multiple SABScatter instances based on the same input object will avoid
    //duplicated resource consumption.
    //
    //The vdoslux parameter has no effect if input is not a VDOS.
    SABScatter( const DI_ScatKnl&, unsigned vdoslux = 3, bool useCache = true );
    SABScatter( SABData &&,
                const VectD& energyGrid = VectD() );
    SABScatter( std::shared_ptr<const SABData>,
                std::shared_ptr<const VectD> energyGrid = nullptr );
    SABScatter( std::shared_ptr<const SAB::SABScatterHelper> );
    SABScatter( SAB::SABScatterHelper&& );

    virtual ~SABScatter();

    double crossSectionNonOriented(double ekin) const final;
    void generateScatteringNonOriented( double ekin, double& angle, double& delta_ekin ) const final;
    void generateScattering( double ekin, const double (&neutron_direction)[3],
                             double (&resulting_neutron_direction)[3], double& delta_ekin ) const final;

  protected:
    struct Impl;
    Pimpl<Impl> m_impl;
    const SAB::SABScatterHelper * m_sh;
  };

}

#endif

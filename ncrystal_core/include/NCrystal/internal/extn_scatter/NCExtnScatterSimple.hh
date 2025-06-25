#ifndef NCrystal_ExtnScatterSimple_hh
#define NCrystal_ExtnScatterSimple_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
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

#include "NCrystal/interfaces/NCProcImpl.hh"
#include "NCrystal/internal/utils/NCExtraTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  class Info;

  namespace Extn {

    class ExtnScatterSimple final : public ProcImpl::ScatterIsotropicMat {
    public:

      //A simple Extinction model for primary extinction in isotropic materials
      //(based on Sabine's model with rectangular tilt function). (fixme revisit
      //description).

      ExtnScatterSimple( PowderBraggInput::Data&& data, Length domainSize );

      const char * name() const noexcept override { return "ExtnScatterSimple"; }
      EnergyDomain domain() const noexcept override;
      CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ) const override;
      ScatterOutcomeIsotropic sampleScatterIsotropic(CachePtr&,
                                                     RNG&,
                                                     NeutronEnergy ) const override;
    protected:
      Optional<std::string> specificJSONDescription() const override;
    private:
      NeutronEnergy m_threshold = NeutronEnergy{kInfinity};
      PowderBraggInput::Data m_data;
      double m_domainSizeAa;
    };

  }

}

#endif

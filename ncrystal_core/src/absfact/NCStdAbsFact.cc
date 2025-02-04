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

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/absoov/NCAbsOOV.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;

//////////////////////////////////////////////////////////////////
//
// The standard absorption factory, handling both single-phase and multi-phase
// materials, with the simple 1/v model. Note that multiple phases might be
// defined at the MatCfg level, but can also appear only at the Info level
// (e.g. if a single NCMAT file produced multiple phases). This factory
// trivially deals with both kinds since the resulting Info object is always
// guaranteed to provide the combined composition and hence average absorption
// cross section.
//
//////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  class NCStdAbsFact final : public FactImpl::AbsorptionFactory {
  public:
    const char * name() const noexcept override { return "stdabs"; }

    MultiPhaseCapability multiPhaseCapability() const override
    {
      return MultiPhaseCapability::Both;
    }

    Priority query( const FactImpl::AbsorptionRequest& ) const override
    {
      return Priority{100};
    }

    ProcImpl::ProcPtr produce( const FactImpl::AbsorptionRequest& cfg ) const override
    {
      return makeSO<AbsOOV>( cfg.info().getXSectAbsorption() );
    }

  };
}

//Finally, a function which can be used to enable the above factory. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_stdabs_factory)()
{
  NC::FactImpl::registerFactory( std::make_unique<NC::NCStdAbsFact>() );
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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

#include "NCrystal/NCFactImpl.hh"
#include "NCrystal/internal/NCAbsOOV.hh"

namespace NC = NCrystal;

namespace NCrystal {

  class NCStdAbsFact final : public FactImpl::AbsorptionFactory {
  public:
    const char * name() const noexcept override { return "stdabs"; }

    Priority query( const MatCfg& cfg ) const override
    {
      auto info = FactImpl::createInfo( cfg );
      if ( ! info->hasXSectAbsorption() )
        return Priority::Unable;
      return Priority{100};
    }

    ProcImpl::ProcPtr produce( const MatCfg& cfg ) const override
    {
      auto info = FactImpl::createInfo( cfg );
      return makeSO<AbsOOV>( info );
    }

  };
}

//Finally, a function which can be used to enable the above factory. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void ncrystal_register_stdabs_factory()
{
  NC::FactImpl::registerFactory( std::make_unique<NC::NCStdAbsFact>(),
                                 NC::FactImpl::RegPolicy::IGNORE_IF_EXISTS );
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

#ifndef NCPlugin_Factory_hh
#define NCPlugin_Factory_hh

#include "NCrystal/NCPluginBoilerplate.hh"//Common stuff (includes NCrystal
                                          //public API headers, sets up
                                          //namespaces and aliases)

namespace NCPluginNamespace {

  //Factory which implements logic of how the physics model provided by the
  //plugin should be combined with existing models in NCrystal:

  class PluginFactory final : public NC::FactImpl::ScatterFactory {
  public:
    const char * name() const noexcept override;

    NC::Priority query( const NC::FactImpl::ScatterRequest& ) const override;
    NC::ProcImpl::ProcPtr produce( const NC::FactImpl::ScatterRequest& ) const override;
  };
}

#endif

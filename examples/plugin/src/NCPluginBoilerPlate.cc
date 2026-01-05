
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

//Define exported symbols used by NCrystal to activate the plugin. This is
//mostly generic, except for the registerPlugin function which can potentially
//be modified in some use-cases.

#define NCPLUGIN_BOILERPLATE_CC
#include "NCrystal/NCPluginBoilerplate.hh"

#include "NCPluginFactory.hh"
#include "NCTestPlugin.hh"

void NCP::registerPlugin()
{
  //This function is required for the plugin to work. It should register
  //factories, and potentially other stuff as appropriate for the plugin (like
  //adding in-mem data files, adding test functions, ...).
  NC::FactImpl::registerFactory(std::make_unique<NCP::PluginFactory>());
  NC::Plugins::registerPluginTestFunction( std::string("test_") + pluginName(),
                                           customPluginTest );
};

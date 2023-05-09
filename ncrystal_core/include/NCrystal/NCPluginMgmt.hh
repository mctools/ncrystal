#ifndef NCrystal_PluginMgmt_hh
#define NCrystal_PluginMgmt_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2023 NCrystal developers                                   //
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

#include "NCrystal/NCDefs.hh"

namespace NCrystal {

  namespace Plugins {
    //
    // The loadDynamicPlugin function attempts to load a plugin by specifying a path to
    // a shared library which defines (with C-mangling) two functions:
    //
    //   void ncplugin_register();
    //   const char* ncplugin_getname();
    //
    // The ncplugin_getname function should simply provide the plugin name
    // whenever called. When ncplugin_register is called, the plugin is expected
    // to enable itself by calling global NCrystal functions (usually this means
    // calling registerFactory(..), but in principle it could also be functions
    // registering i/o managers, in-memory file data, or changing the default
    // RNG).
    //
    // The loadBuiltinPlugin does the same, but directly providing the plugin
    // name and registration function to call.
    //
    // Normally the plugin interface is not accessed. Rather, dynamic plugins
    // are enabled by adding them to the NCRYSTAL_PLUGIN_LIST environment
    // variable, while the loading of builtin plugins is handled by NCrystal's
    // cmake configuration.

    //
    enum class PluginType { Dynamic, Builtin, Undefined };
    struct NCRYSTAL_API PluginInfo {
      //Simple information about loaded plugin. For now this is just the path to
      //the library and the pluginname as returned by ncplugin_getname (see
      //below). For builtin plugins the fileName will be empty.
      std::string pluginName;
      std::string fileName;
      PluginType pluginType = PluginType::Undefined;
    };

    //Load plugins (returns the resulting plugin info):
    NCRYSTAL_API PluginInfo loadDynamicPlugin( std::string path_to_shared_lib );
    NCRYSTAL_API PluginInfo loadBuiltinPlugin( std::string pluginName,
                                               std::function<void()> regfct );

    //Query loaded plugins:
    NCRYSTAL_API std::vector<PluginInfo> loadedPlugins();

    //Call this to ensure plugins (both builtin and those in
    //NCRYSTAL_PLUGIN_LIST) are loaded. Multiple calls to this function will
    //have no effect, and it will be called automatically when users query the
    //global createXXX(..) functions from NCFactory.hh:
    NCRYSTAL_API void ensurePluginsLoaded();

  }
}

#endif

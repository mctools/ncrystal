#ifndef NCrystal_PluginBoilerplate_hh
#define NCrystal_PluginBoilerplate_hh

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

//Header file which should be included by code creating plugins (and nowhere
//else). It contains macros for creating required namespaces and symbols. It
//exists purely as a matter of convenience, to avoid duplicating code in a lot
//of plugins, and to ensure all plugins have a C++ namespace which contains the
//name of the plugin (to avoid symbol clashes). Of course, usage of the present
//header file is not strictly required: the actual basic requirements of a
//plugin is just a shared object (module) which defines two functions with
//C-mangling: void ncplugin_register(void) and const char * ncplugin_getname().

#ifndef NCPLUGIN_NAME
#  error NCPLUGIN_NAME undefined
#endif

//Internal macros needed below:
#define ncpluginhh_str2(s) #s
#define ncpluginhh_str(s) ncpluginhh_str2(s)
#define ncpluginhh_nsgen0(m) NCrystalPlugin_##m
#define ncpluginhh_nsgen(m) ncpluginhh_nsgen0(m)

//Macros which can/should be used in plugin code. In fact, ALL plugin code
//should be in the NCPluginNamespace!
#define NCPLUGIN_NAME_CSTR ncpluginhh_str(NCPLUGIN_NAME)
#define NCPluginNamespace ncpluginhh_nsgen(NCPLUGIN_NAME)

//All plugins must implement the following registerPlugin function:
namespace NCPluginNamespace {
  void registerPlugin();
}

//For convenience we include the public NCrystal API and create a few aliases:
#include "NCrystal/NCrystal.hh"

namespace NC = NCrystal;
namespace NCP = NCPluginNamespace;

//In addition to the NCPLUGIN_NAME_CSTR macro above, the plugin name is also
//available via functions return strings (including upper/lower cased):
namespace NCPluginNamespace {
  const std::string& pluginName();
  const std::string& pluginNameUpperCase();
  const std::string& pluginNameLowerCase();
}

////////////////////////////////////////////////////////////////////////////////////
//Exactly one .cc file in a plugin must define NCPLUGIN_BOILERPLATE_CC before
//including the file, which will then generate the C-mangled hooks needed for
//dynamic access to the plugin:
#if defined(NCPLUGIN_BOILERPLATE_CC)
//C-mangled hooks (only for dynamic plugins)
#  if defined (_WIN32) || defined (__CYGWIN__) || defined (WIN32)
#    define ncpluginhh_export __declspec(dllexport)
#  elif defined(__GNUC__) || defined(__clang__)
#    define ncpluginhh_export __attribute__ ((visibility ("default")))
#  else
#    define ncpluginhh_export
#  endif

extern "C" {
  ncpluginhh_export void ncplugin_register()
  {
    NCP::registerPlugin();
  }
  ncpluginhh_export const char * ncplugin_getname()
  {
    return NCPLUGIN_NAME_CSTR;
  }
}

const std::string& NCP::pluginName()
{
  static std::string s_name(NCPLUGIN_NAME_CSTR);
  return s_name;
}

const std::string& NCP::pluginNameUpperCase()
{
  static std::string s_name = [](){
    std::string s( NCPLUGIN_NAME_CSTR );
    for (auto& c : s)
      if ( c >= 'a' && c <= 'z' )
        c -= 32;
    return s;
  }();
  return s_name;
}

const std::string& NCP::pluginNameLowerCase()
{
  static std::string s_name = [](){
    std::string s( NCPLUGIN_NAME_CSTR );
    for (auto& c : s)
      if ( c >= 'A' && c <= 'Z' )
        c += 32;
    return s;
  }();
  return s_name;
}

#endif//NCPLUGIN_BOILERPLATE_CC

#endif

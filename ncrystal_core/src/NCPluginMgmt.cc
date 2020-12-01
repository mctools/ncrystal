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

#include "NCrystal/NCPluginMgmt.hh"
#include "NCrystal/internal/NCDynLoader.hh"
#include "NCrystal/internal/NCString.hh"
#include <iostream>

//MT TODO: Do we need to make these thread-safe?

namespace NC = NCrystal;
namespace NCP = NCrystal::Plugins;

namespace NCrystal {
  namespace Plugins {
    namespace {
      std::mutex& getPluginMgmtMutex()
      {
        static std::mutex mtx;
        return mtx;
      }
      std::vector<PluginInfo>& getPLList()
      {
        static std::vector<PluginInfo> thelist;
        return thelist;
      }

      void actualLoadPlugin( PluginInfo pinfo,
                             std::function<void()> regfct )
      {
        //Mutex is already locked when this is called!
        nc_assert_always(pinfo.pluginType==PluginType::Dynamic||pinfo.pluginType==PluginType::Builtin);
        bool verbose = ncgetenv_bool("DEBUG_PLUGIN");
        std::string ptypestr(pinfo.pluginType==PluginType::Dynamic?"dynamic":"builtin");
        if (verbose)
          std::cout<<"NCrystal: Loading "<<ptypestr<<" plugin \""<<pinfo.pluginName<<"\"."<<std::endl;

        for ( const auto& pl : getPLList() ) {
          if ( pl.pluginName == pinfo.pluginName )
            NCRYSTAL_THROW2(CalcError,"ERROR: attempting to load plugin named \""<<pinfo.pluginName<<"\" more than once!");
        }

        try {
          regfct();
        } catch (std::exception& e) {
          std::cout<<"NCrystal ERROR: Problems while loading plugin!"<<std::endl;
          throw;
        }

        getPLList().push_back(pinfo);

        if (verbose)
          std::cout<<"NCrystal: Done loading plugin \""<<pinfo.pluginName<<"\"."<<std::endl;
      }
    }
  }
}

namespace NCrystal {
  namespace Plugins {
    namespace {
      PluginInfo loadDynamicPluginImpl( std::string path_to_shared_lib,
                                        std::string pluginName,
                                        std::string regfctname )
      {
        PluginInfo pinfo;
        pinfo.pluginType = PluginType::Dynamic;
        pinfo.fileName = path_to_shared_lib;
        pinfo.pluginName = pluginName;

        std::lock_guard<std::mutex> guard(getPluginMgmtMutex());

        if (ncgetenv_bool("DEBUG_PLUGIN"))
          std::cout<<"NCrystal: Attempting to loading dynamic library with plugin: "<<pinfo.fileName<<std::endl;

        DynLoader dl( pinfo.fileName,
                      DynLoader::ScopeFlag::local,//reduce symbol interference between plugins
                      DynLoader::LazyFlag::now );//better get the error upon load

        if (pinfo.pluginName.empty())
          pinfo.pluginName = dl.getFunction<const char*()>("ncplugin_getname")();
        auto regfct = dl.getFunction<void()>(regfctname);

        dl.doNotClose();//avoid dlclose (sadly, this is leaky but otherwise we
        //might get segfaults at programme shutdown...)

        actualLoadPlugin( pinfo, std::move(regfct) );

        return pinfo;
      }
    }
  }
}

NCP::PluginInfo NCP::loadDynamicPlugin( std::string path_to_shared_lib )
{
  return loadDynamicPluginImpl( path_to_shared_lib, "", "ncrystal_register" );
}

NCP::PluginInfo NCP::loadBuiltinPlugin( std::string pluginName,
                                        std::function<void()> regfct )
{
  PluginInfo pinfo;
  pinfo.pluginType = PluginType::Builtin;
  pinfo.pluginName = pluginName;
  std::lock_guard<std::mutex> guard(getPluginMgmtMutex());
  actualLoadPlugin( pinfo, std::move(regfct) );
  return pinfo;
}

std::vector<NCP::PluginInfo> NCP::loadedPlugins()
{
  NCP::ensurePluginsLoaded();
  std::vector<NCP::PluginInfo> result;
  {
    std::lock_guard<std::mutex> guard(getPluginMgmtMutex());
    result = getPLList();
  }
  return result;
}

//Automatic enablement of .nxs/.laz support is controlled via the
//NCRYSTAL_ENABLE_NXSLAZ macro, and the nxs/laz factories must obviously have
//been built for this to work. Support for .ncmat files is on the other hand
//ALWAYS enabled by default, unless disabled with the macro
//NCRYSTAL_DISABLE_NCMAT. Likewise, the standard scatter/absorption factories
//are always enabled, unless respectively NCRYSTAL_DISABLE_STDSCAT and
//NCRYSTAL_DISABLE_STDABS is defined.

#ifndef NCRYSTAL_DISABLE_STDSCAT
extern "C" void ncrystal_register_stdscat_factory();
#endif
#ifndef NCRYSTAL_DISABLE_STDABS
extern "C" void ncrystal_register_stdabs_factory();
#endif
#ifndef NCRYSTAL_DISABLE_NCMAT
extern "C" void ncrystal_register_ncmat_factory();
#endif
#ifdef NCRYSTAL_ENABLE_NXSLAZ
extern "C" void ncrystal_register_nxslaz_factories();
#endif

//If NCRYSTAL_HAS_BUILTIN_PLUGINS is defined, it is assumed that the code is
//being built with a void NCrystal::provideBuiltinPlugins() function defined in
//another module. The envisioned use-case is that plugins can be available in
//two modes: 1) as a separate shared library with standard symbols for loading
//the library, and 2) being built into a given NCrystal installation at compile
//time.
#ifdef NCRYSTAL_HAS_BUILTIN_PLUGINS
namespace NCrystal {
  void provideBuiltinPlugins();
}
#endif

void NCP::ensurePluginsLoaded()
{
  static std::atomic<bool> first(true);
  bool btrue(true);
  if (!first.compare_exchange_strong(btrue,false))
    return;//only the first call will proceed past this point. If we instead
           //used a mutex, we would get a deadlock (e.g. when registerFactory
           //calls ensurePluginsLoaded before checking that the new factory name
           //is unique). The MT-safety of the calls below should be the same as
           //for any other call to loadXXXPlugin.

  static bool done = false;
  if (done)
    return;
  done = true;

  //Standard plugins, always as builtin plugins:
#ifndef NCRYSTAL_DISABLE_STDSCAT
  loadBuiltinPlugin("stdscat",ncrystal_register_stdscat_factory);
#endif
#ifndef NCRYSTAL_DISABLE_STDABS
  loadBuiltinPlugin("stdabs",ncrystal_register_stdabs_factory);
#endif
#ifndef NCRYSTAL_DISABLE_NCMAT
  loadBuiltinPlugin("stdncmat",ncrystal_register_ncmat_factory);
#endif
#ifdef NCRYSTAL_ENABLE_NXSLAZ
  loadBuiltinPlugin("nxslaz",ncrystal_register_ncmat_factory);//TODO: As external plugin?
#endif

  //Static custom (builtin) plugins:
#ifdef NCRYSTAL_HAS_BUILTIN_PLUGINS
  provideBuiltinPlugins();
#endif

  //Dynamic custom plugins, as indicated by environment variable:
  for (auto& pluginlib : split2(ncgetenv("PLUGIN_LIST"),0,':')) {
    trim(pluginlib);
    if (pluginlib.empty())
      continue;
    Plugins::loadDynamicPlugin(pluginlib);
  }
}

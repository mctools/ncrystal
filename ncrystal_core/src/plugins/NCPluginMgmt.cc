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

#include "NCrystal/plugins/NCPluginMgmt.hh"
#include "NCrystal/internal/utils/NCDynLoader.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCFileUtils.hh"

//MT TODO: Do we need to make these thread-safe?

namespace NC = NCrystal;
namespace NCP = NCrystal::Plugins;

namespace NCRYSTAL_NAMESPACE {
  namespace Plugins {

    namespace {

      struct PluginTestFctsDB {
        std::vector<std::pair<std::string,voidfct_t>> fcts;
        std::mutex mtx;
      };

      PluginTestFctsDB& getPluginTestFctsDB()
      {
        static PluginTestFctsDB db;
        return db;
      }

      void runRegisteredPluginTestFunctions( std::size_t skip = 0 )
      {
        auto testfcts = getRegisteredPluginTestFunctions();
        auto it = std::next(testfcts.begin(),skip);
        auto itE = testfcts.end();
        for ( ; it!=itE; ++it ) {
          auto& e = *it;
          NCRYSTAL_MSG("Launching plugin test function \""<<e.first<<"\":");
          e.second();
          NCRYSTAL_MSG("End of plugin test function \""<<e.first<<"\".");
        }
      }

      Optional<std::string> executeCommandAndCaptureOutput( const char* cmd )
      {
        constexpr auto nbuf = 4096;
        std::array<char, nbuf> buffer;
        Optional<std::string> result;

#if ( defined (_WIN32) || defined (WIN32) )
        FILE* pipe = _popen(cmd, "r");
#else
        FILE* pipe = popen(cmd, "r");
#endif
        if (pipe) {
          result.emplace();
          while (fgets(buffer.data(), nbuf, pipe) != NULL)
            result.value() += buffer.data();
#if ( defined (_WIN32) || defined (WIN32) )
          auto returnCode = _pclose(pipe);
#else
          auto returnCode = pclose(pipe);
#endif
          if ( returnCode != 0 )
            result.reset();
        }
        return result;
      }

#if ( defined (_WIN32) || defined (WIN32) )
      bool win_pluginmgr_cmd_exists()
      {
        //Check PATH directly, since "2>/dev/null" does not work on windows.
        const char * raw_path = std::getenv("PATH");
        if (!raw_path)
          return false;
        std::string bn1( "ncrystal-pluginmanager.exe" );
        std::string bn2( "ncrystal-pluginmanager.bat" );
        for ( auto& e : StrView(raw_path).splitTrimmedNoEmpty(';') ) {
          std::string dirname = e.to_string();
          if ( file_exists( path_join(dirname,bn1) ) )
            return true;
          if ( file_exists( path_join(dirname,bn2) ) )
            return true;
        }
        return false;
      }
#endif

      struct PlugMgrResults {
        SmallVector<PairSS,4> plugin_datadirs;
        VectS dynlibs;
      };
      PlugMgrResults queryPluginManagerCmd() {
        //Invokes ncrystal-pluginmanager to determine any plugins to use, thus
        //supporting the ability for "pip install ./my/plugin" to work
        //immediately with no further issues.
        PlugMgrResults result;
#ifdef NCRYSTAL_DISABLE_DYNLOAD
        return result;//all plugin loading forbidden
#endif
#ifdef NCRYSTAL_DISABLE_CMDLINEPLUGINMGR
        return result;//just the ncrystal-pluginmanager disabled
#endif
#if ( defined (_WIN32) || defined (WIN32) )
        if (!win_pluginmgr_cmd_exists())
          return result;
        auto out = executeCommandAndCaptureOutput( "ncrystal-pluginmanager" );
#else
        auto out = executeCommandAndCaptureOutput( "ncrystal-pluginmanager 2>/dev/null" );
#endif
        if (out.has_value()) {
          auto parts = StrView(out.value()).splitTrimmedNoEmpty(';');
          for ( auto& e : parts ) {
            if ( e.startswith(":DATA:") ) {
              //NB: Paths can contain ':' on windows (e.g. "C:...")
              auto ee = e.substr(6);
              auto idx = ee.find(':');
              StrView p1, p2;
              if ( idx != StrView::npos ) {
                p1 = ee.substr(0,idx).trimmed();
                p2 = ee.substr(idx+1).trimmed();
              }
              if ( !p1.has_value() || p1.empty() || p2.empty() )
                NCRYSTAL_THROW2( BadInput,
                                 "Invalid ncrystal-pluginmanager"
                                 " output in entry \""<<e<<"\"" );
              result.plugin_datadirs.emplace_back( p1.to_string(),
                                                   p2.to_string() );
            } else {
              result.dynlibs.emplace_back( e.to_string() );
            }
          }
        }
        return result;
      }

      std::mutex& getPluginMgmtMutex()
      {
        static std::mutex mtx;
        return mtx;
      }
      std::vector<PluginInfo>& getSharedLibPLList()
      {
        static std::vector<PluginInfo> thelist;
        return thelist;
      }

      void actualLoadPlugin( PluginInfo pinfo,
                             voidfct_t regfct )
      {
        //Mutex is already locked when this is called!
        nc_assert_always( pinfo.pluginType==PluginType::Dynamic
                          || pinfo.pluginType==PluginType::Builtin );
        bool verbose = ncgetenv_bool("DEBUG_PLUGIN");
        std::string ptypestr( pinfo.pluginType==PluginType::Dynamic
                              ? "dynamic" : "builtin" );
        if (verbose)
          NCRYSTAL_MSG("Loading "<<ptypestr<<" plugin \""<<pinfo.pluginName<<"\".");

        for ( const auto& pl : getSharedLibPLList() ) {
          if ( pl.pluginName == pinfo.pluginName )
            NCRYSTAL_THROW2(CalcError,"ERROR: attempting to load plugin named \""<<pinfo.pluginName<<"\" more than once!");
        }

        try {
          regfct();
        } catch (std::exception& e) {
          NCRYSTAL_RAWOUT("NCrystal ERROR: Problems while loading plugin!\n");
          throw;
        }

        getSharedLibPLList().push_back(pinfo);

        if (verbose)
          NCRYSTAL_MSG("Done loading plugin \""<<pinfo.pluginName<<"\".");
      }
      struct PluginsDataDirDBState {
        std::mutex mtx;
        std::vector<PairSS> data;
      };
      PluginsDataDirDBState& getPluginDataDirDBState()
      {
        static PluginsDataDirDBState db;
        return db;
      }
      void appendPluginDataDirDB( SmallVector<PairSS,4>&& new_entries )
      {
        if (new_entries.empty())
          return;
        auto& db = getPluginDataDirDBState();
        NCRYSTAL_LOCK_GUARD(db.mtx);
        db.data.reserve( db.data.size() + new_entries.size() );
        for (auto& e : new_entries )
          db.data.push_back( std::move(e) );
      }

    }//end anon namespace
    namespace detail {
      std::vector<PairSS> getPluginDataDirDB()
      {
        auto& db = getPluginDataDirDBState();
        NCRYSTAL_LOCK_GUARD(db.mtx);
        return db.data;
      }
    }
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace Plugins {
    namespace {
      std::string resolvePathToShlib( std::string path )
      {
        std::string found;
        if ( endswith(path,".SODYLIB") ) {
          //Special ending, allowing to try both .so/.dylib
          std::string path_base = path.substr(0,path.size()-8);
          for ( auto& e : { ".so", ".dylib", ".dll", ".DLL" } ) {
            if ( file_exists( path_base + e ) ) {
              found = path_base + e;
              break;
            }
          }
        } else {
          if ( file_exists( path ) )
            found = std::move(path);
        }
        if ( found.empty() )
          NCRYSTAL_THROW2( BadInput, "Could not find dynamic plugin"
                           " shared library: " << path );
        return found;
      }

      PluginInfo loadDynamicPluginImpl( std::string path_to_shared_lib,
                                        std::string pluginName = "",
                                        std::string regfctname = "ncplugin_register" )
      {
        path_to_shared_lib = resolvePathToShlib(path_to_shared_lib);

        PluginInfo pinfo;
        pinfo.pluginType = PluginType::Dynamic;
        pinfo.fileName = path_to_shared_lib;
        pinfo.pluginName = pluginName;

        if (ncgetenv_bool("DEBUG_PLUGIN"))
          NCRYSTAL_MSG("Attempting to loading dynamic library with plugin: "<<pinfo.fileName);

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

      PluginInfo loadBuiltinPluginImpl( std::string pluginName,
                                        voidfct_t regfct )
      {
        PluginInfo pinfo;
        pinfo.pluginType = PluginType::Builtin;
        pinfo.pluginName = pluginName;
        actualLoadPlugin( pinfo, std::move(regfct) );
        return pinfo;
      }
    }
  }
}

NCP::PluginInfo NCP::loadBuiltinPlugin( std::string pluginName,
                                        voidfct_t regfct )
{
  NCRYSTAL_LOCK_GUARD(getPluginMgmtMutex());
  return loadBuiltinPluginImpl( std::move(pluginName),
                                std::move(regfct) );

}

NCP::PluginInfo NCP::loadDynamicPlugin( std::string path_to_shared_lib )
{
  PluginInfo res;
  const bool run_test_functions = ncgetenv_bool("PLUGIN_RUNTESTS");
  std::size_t ntestfcts = ( run_test_functions
                            ? 0
                            : getRegisteredPluginTestFunctions().size() );
  {
    NCRYSTAL_LOCK_GUARD(getPluginMgmtMutex());
    res = loadDynamicPluginImpl( path_to_shared_lib );
  }
  if ( run_test_functions )
    runRegisteredPluginTestFunctions( ntestfcts );
  return res;
}

std::vector<NCP::PluginInfo> NCP::loadedPlugins()
{
  NCP::ensurePluginsLoaded();
  std::vector<NCP::PluginInfo> result;
  {
    NCRYSTAL_LOCK_GUARD(getPluginMgmtMutex());
    result = getSharedLibPLList();
  }
  std::vector<PairSS> datadirdb;
  {
    datadirdb = detail::getPluginDataDirDB();
  }

  for ( auto& e : datadirdb ) {
    //Also add any purely static data plugins (i.e. pure python plugins):
    bool found(false);
    for ( auto& r : result ) {
      if ( e.first == r.pluginName ) {
        found = true;
        break;
      }
    }
    if ( !found ) {
      result.emplace_back();
      result.back().pluginName = e.first;
      result.back().fileName = e.second;
      result.back().pluginType = PluginType::Dynamic;
    }
  }

  return result;
}

#ifndef NCRYSTAL_SIMPLEBUILD_DEVEL_MODE

#ifndef NCRYSTAL_DISABLE_STDDATASOURCES
extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_stddatasrc_factory)();
#endif
#ifndef NCRYSTAL_DISABLE_STDSCAT
extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_stdscat_factory)();
#endif
#ifndef NCRYSTAL_DISABLE_STDMPSCAT
extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_stdmpscat_factory)();
#endif
#ifndef NCRYSTAL_DISABLE_STDABS
extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_stdabs_factory)();
#endif
#ifndef NCRYSTAL_DISABLE_STDNCMAT
extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_stdncmat_factory)();
#endif
#ifndef NCRYSTAL_DISABLE_STDLAZ
extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_stdlaz_factory)();
#endif
#ifndef NCRYSTAL_DISABLE_QUICKFACT
extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_quick_factory)();
#endif
#ifndef NCRYSTAL_DISABLE_EXPERIMENTALSCATFACT
extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_experimentalscatfact)();
#endif

#endif

#ifdef NCRYSTAL_SIMPLEBUILD_DEVEL_MODE
namespace NCRYSTAL_NAMESPACE {
  namespace {
    std::string getSBLDPkgLib(const char* pkgname) {
      auto sbld_libdir_raw = std::getenv("SBLD_LIB_DIR");
      if ( !sbld_libdir_raw)
        NCRYSTAL_THROW(CalcError,"SBLD_LIB_DIR not set");
      auto soprefix = std::string(sbld_libdir_raw) + "/libPKG__";
      std::string lib(soprefix + pkgname);
      if ( file_exists(lib+".so") ) {
        lib+=".so";
      } else if (file_exists(lib+".dylib")) {
        lib+=".dylib";
      } else {
        NCRYSTAL_THROW2(FileNotFound,"Could not find"
                        " sbld pkg library (\""
                        <<pkgname<<"\"). Is the package enabled?");
      }
      return lib;
    }
  }
}
#endif

void NCP::ensurePluginsLoaded()
{
  bool run_test_functions = false;
  std::size_t ntestfcts = 0;

  {
    NCRYSTAL_LOCK_GUARD(getPluginMgmtMutex());
    static bool first = true;
    if (!first)
      return;
    first = false;
    run_test_functions = ncgetenv_bool("PLUGIN_RUNTESTS");
    if ( run_test_functions )
      ntestfcts = getRegisteredPluginTestFunctions().size();

#ifdef NCRYSTAL_SIMPLEBUILD_DEVEL_MODE
    //Standard plugins, dynamic in simplebuild mode, except for datasources:
    auto loadstdplugin = []( const char* pkgname,
                             const char* pluginname,
                             std::string regfctname )
    {
      regfctname = std::string(ncrystal_xstr(NCRYSTAL_C_NAMESPACE)) + regfctname;
      loadDynamicPluginImpl( getSBLDPkgLib(pkgname), pluginname, regfctname);
    };
    loadstdplugin("NCFactories","stddatasrc","register_stddatasrc_factory");
    loadstdplugin("NCScatFact","stdscat","register_stdscat_factory");
    loadstdplugin("NCScatFact","stdmpscat","register_stdmpscat_factory");
    loadstdplugin("NCExperimental","stdexpscat","register_experimentalscatfact");
    loadstdplugin("NCFactory_Laz","stdlaz","register_stdlaz_factory");
    loadstdplugin("NCAbsFact","stdabs","register_stdabs_factory");
    loadstdplugin("NCFactory_NCMAT","stdncmat","register_stdncmat_factory");
    loadstdplugin("NCQuickFact","stdquick","register_quick_factory");
    //loadstdplugin("NCNXSFactories","legacynxslaz","register_nxslaz_factories");
#endif

#ifndef NCRYSTAL_SIMPLEBUILD_DEVEL_MODE
    //Standard plugins, always as builtin plugins:
#  ifndef NCRYSTAL_DISABLE_STDDATASOURCES
    loadBuiltinPluginImpl("stddatasrc",NCRYSTAL_APPLY_C_NAMESPACE(register_stddatasrc_factory));
#  endif
#  ifndef NCRYSTAL_DISABLE_STDSCAT
    loadBuiltinPluginImpl("stdscat",NCRYSTAL_APPLY_C_NAMESPACE(register_stdscat_factory));
#  endif
#  ifndef NCRYSTAL_DISABLE_STDMPSCAT
    loadBuiltinPluginImpl("stdmpscat",NCRYSTAL_APPLY_C_NAMESPACE(register_stdmpscat_factory));
#  endif
#  ifndef NCRYSTAL_DISABLE_EXPERIMENTALSCATFACT
    loadBuiltinPluginImpl("stdexpscat",NCRYSTAL_APPLY_C_NAMESPACE(register_experimentalscatfact));
#  endif
#  ifndef NCRYSTAL_DISABLE_STDLAZ
    loadBuiltinPluginImpl("stdlaz",NCRYSTAL_APPLY_C_NAMESPACE(register_stdlaz_factory));
#  endif
#  ifndef NCRYSTAL_DISABLE_STDABS
    loadBuiltinPluginImpl("stdabs",NCRYSTAL_APPLY_C_NAMESPACE(register_stdabs_factory));
#  endif
#  ifndef NCRYSTAL_DISABLE_NCMAT
    loadBuiltinPluginImpl("stdncmat",NCRYSTAL_APPLY_C_NAMESPACE(register_stdncmat_factory));
#  endif
#  ifndef NCRYSTAL_DISABLE_QUICKFACT
    loadBuiltinPluginImpl("stdquick",NCRYSTAL_APPLY_C_NAMESPACE(register_quick_factory));
#  endif
#endif

    //Dynamic custom plugins, as indicated by environment variable:
#ifdef NCRYSTAL_SIMPLEBUILD_DEVEL_MODE
    bool first_plugin = false;
#endif

    //Get list of dynamic plugins to load via the ncrystal-pluginmanager command
    //and NCRYSTAL_PLUGIN_LIST environment variables:
    auto plugmgrcmd_results = queryPluginManagerCmd();
    VectS dynplugin_list = std::move(plugmgrcmd_results.dynlibs);
    for ( auto& e : split2(ncgetenv("PLUGIN_LIST"),0,':') )
      dynplugin_list.push_back( trim2(std::move(e)) );

    //However, NCRYSTAL_DISABLE_DYNLOAD env var and NCRYSTAL_DISABLE_DYNLOAD
    //define can be used to prevent any third party dynamic plugins:
#if defined(NCRYSTAL_DISABLE_DYNLOAD)
    const bool cfg_disable_dynload = true;
#else
    const bool cfg_disable_dynload = ncgetenv_bool("DISABLE_DYNLOAD");

#endif
    if ( cfg_disable_dynload )
      dynplugin_list.clear();

    std::set<std::string> dynplugins_already_loaded;
    for ( auto& pluginlib : dynplugin_list ) {
      if ( pluginlib.empty() || dynplugins_already_loaded.count(pluginlib) )
        continue;
      dynplugins_already_loaded.insert(pluginlib);

#ifdef NCRYSTAL_SIMPLEBUILD_DEVEL_MODE
      //I am not 100% sure the following is still needed, but keeping it just in
      //case for now:
      if (first_plugin) {
        first_plugin=false;
        DynLoader( getSBLDPkgLib("NCFactories"),
                   DynLoader::ScopeFlag::global ).doNotClose();
      }
#endif
      Plugins::loadDynamicPluginImpl(pluginlib);
    }

    //Data directories from ncrystal-pluginmanager:
    SmallVector<PairSS,4> ddirs = std::move( plugmgrcmd_results.plugin_datadirs );
    //Data directories from env var: TODO: Consider alternative to ':' in all
    //path vars, since they do not work on windows (since there are colons in
    //places like: "C:\...")
    for ( auto& e : split2(ncgetenv("PLUGIN_DATADIRS"),0,':') ) {
      auto p = StrView(e).splitTrimmedNoEmpty('@');
      if ( p.empty() )
        continue;
      if ( p.size()!=2 || p.at(0).empty() || p.at(1).empty() )
        NCRYSTAL_THROW2( BadInput,
                         "Invalid NCRYSTAL_PLUGIN_DATADIRS"
                         " content in entry \""<<e<<"\"" );
      ddirs.emplace_back( p.at(0).to_string(),
                          p.at(1).to_string() );
    }
    if ( !cfg_disable_dynload )
      appendPluginDataDirDB( std::move( ddirs ) );

  }//release mutex

  auto required_plugins = ncgetenv("REQUIRED_PLUGINS");
  if (!required_plugins.empty()) {
    auto avail_plugins = loadedPlugins();
    for ( auto& required_plugin : split2(required_plugins,0,':') ) {
      std::string found_ptypestr;
      for ( auto& pinfo : avail_plugins ) {
        if ( pinfo.pluginName == required_plugin ) {
          found_ptypestr = ( pinfo.pluginType == PluginType::Dynamic
                             ? "dynamic" : "builtin" );
          break;
        }
      }
      if (found_ptypestr.empty())
        NCRYSTAL_THROW2( LogicError, "Required plugin was not loaded: \""
                         << required_plugin << '"' );
      NCRYSTAL_MSG("Required plugin \""<<required_plugin
                   <<"\" was indeed available ("<<found_ptypestr<<").");
    }
  }

  if ( run_test_functions )
    runRegisteredPluginTestFunctions( ntestfcts );
}

void NCP::registerPluginTestFunction( std::string test_name,
                                      voidfct_t test_fct )
{
  auto& db = getPluginTestFctsDB();
  NCRYSTAL_LOCK_GUARD( db.mtx );
  db.fcts.emplace_back( std::move(test_name), std::move(test_fct) );
}

std::vector<std::pair<std::string,NC::voidfct_t>>
NCP::getRegisteredPluginTestFunctions()
{
  std::vector<std::pair<std::string,NC::voidfct_t>> res;
  {
    auto& db = getPluginTestFctsDB();
    NCRYSTAL_LOCK_GUARD( db.mtx );
    res = db.fcts;
  }
  return res;
}

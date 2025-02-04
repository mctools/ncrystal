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

#include "NCrystal/factories/NCDataSources.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCFileUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/plugins/NCPluginMgmt.hh"

namespace NC = NCrystal;
namespace NCD = NCrystal::DataSources;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    struct AtomicBoolSetValResult { bool wasChanged; };
    AtomicBoolSetValResult setValue( std::atomic<bool>& var, bool newvalue ) {
      bool b = newvalue;
      bool oldval = var.exchange(b);
      return { newvalue != oldval };
    }
  }
  namespace DataSources {
    using BrowseEntry = FactImpl::TextDataFactory::BrowseEntry;
    std::vector<BrowseEntry> browseDir( const std::string& dirname, Priority priority ) {
      //Browse for known extensions in directory and any direct sub-dirs without
      //"." in the name.
      std::vector<BrowseEntry> out;
      out.reserve(16);
      VectS subdirs;
      VectS extensionsV = recognisedFileExtensions();
      const std::set<std::string> extensions(extensionsV.begin(),extensionsV.end());
      for ( const auto& f : ncglob(dirname+"/*") ) {
        auto bn = basename(f);
        if ( !contains(bn,'.') ) {
          //might be subdir (or file without '.' in the name, but the glob below
          //will get no hits in that case).
          subdirs.push_back(bn);
        } else if ( extensions.count(getfileext(bn))!=0 ) {
          out.push_back( { bn, dirname, priority } );
        }
      }
      //Deal with subdirs:
      for ( const auto& subdir : subdirs ) {
        auto subdir_fullpath = path_join( dirname, subdir );
        for ( const auto& f : ncglob(subdir_fullpath+"/*.*") ) {
          if ( extensions.count(getfileext(f))!=0 )
            out.push_back( { path_join(subdir,basename(f)), dirname, priority } );
        }
      }
      out.shrink_to_fit();
      return out;
    }

    namespace {
      static std::atomic<bool> s_was_called_enableAbsolutePaths = {false};
      static std::atomic<bool> s_was_called_enableRelativePaths = {false};
      static std::atomic<bool> s_was_called_enableStandardDataLibrary = {false};
      static std::atomic<bool> s_was_called_enableStandardSearchPath = {false};
      static std::atomic<bool> s_was_called_enablePluginSearchPaths = {false};
    }

    static const auto factNameStdLib = "stdlib";
    static const auto factNameStdSearchPath = "stdpath";
    static const auto factNamePluginSearchPaths = "plugins";
    static const auto factNameRelPath = "relpath";
    static const auto factNameAbsPath = "abspath";
    static const auto factNameVirtualFiles = "virtual";
    static const auto factNameCustomSearchDirs = "customdirs";

    class TDFact_AbsPath final : public FactImpl::TextDataFactory {
    public:
      const char * name() const noexcept override
      {
        return factNameAbsPath;
      }
      Priority query( const TextDataPath& p ) const override
      {
        return ( path_is_absolute(p.path()) && file_exists(p.path()) )
          ? Priority{default_priority_abspath}
          : Priority::Unable;
      }
      TextDataSource produce( const TextDataPath& p ) const override
      {
        nc_assert( path_is_absolute(p.path()) && file_exists(p.path()) );
        return TextDataSource::createFromOnDiskPath( p.path() );
      }
      //Can't browse absolute paths:
      std::vector<BrowseEntry> browse() const override { return {}; }
    };

    class TDFact_PluginDirs final : public FactImpl::TextDataFactory {
      struct ParsedReq {
        StrView plugin_name;
        StrView file_name;
      };
      ParsedReq parsePath( const TextDataPath& p ) const
      {
        ParsedReq result;
        StrView pathstr( p.path() );
        if ( pathstr.contains_any(":#~\\") )
          return result;
        auto idx = pathstr.find("/");
        if ( idx == StrView::npos )
          return result;
        auto pluginname = pathstr.substr(0,idx).trimmed();
        auto filename = pathstr.substr(idx+1).trimmed();
        if ( filename.empty() || pluginname.empty() || filename.contains('/') )
          return result;
        result.plugin_name = pluginname;
        result.file_name = filename;
        return result;
      }
      std::string lookupFile( const TextDataPath& p ) const
      {
        auto request = parsePath(p);
        if (!request.plugin_name.has_value())
          return {};
        auto db = Plugins::detail::getPluginDataDirDB();
        if ( db.empty() )
          return {};
        std::string filename = request.file_name.to_string();
        for ( auto& e : db ) {
          if ( request.plugin_name == e.first ) {
            auto full_path = path_join( e.second, filename );
            if ( file_exists( full_path ) )
              return full_path;
          }
        }
        return {};
      }
    public:
      const char * name() const noexcept override
      {
        return factNamePluginSearchPaths;
      }
      Priority query( const TextDataPath& p ) const override
      {
        auto full_path = lookupFile(p);
        return ( full_path.empty()
                 ? Priority::Unable
                 : Priority{Priority::OnlyOnExplicitRequest} );
      }
      TextDataSource produce( const TextDataPath& p ) const override
      {
        auto full_path = lookupFile(p);
        if (full_path.empty())
          NCRYSTAL_THROW2( DataLoadError,"File disappeared suddenly"
                           " during request: " << p );
        return TextDataSource::createFromOnDiskPath( full_path );
      }
      std::vector<BrowseEntry> browse() const override {
        //In case a plugin has more than one directory with the same filename
        //inside, we would only serve the first entry. We faithfully reproduce
        //this in the browsed list, but also emit a warning:
        std::vector<BrowseEntry> out;
        std::set<PairSS> seen;
        auto db = Plugins::detail::getPluginDataDirDB();
        const auto priority = Priority{Priority::OnlyOnExplicitRequest};
        for ( auto& db_entry : db ) {
          for ( auto& be : browseDir( db_entry.second, priority ) ) {
            PairSS key{ be.name, be.source };
            std::string fullname = db_entry.first + "/" + be.name;
            if ( seen.count( key ) ) {
              NCRYSTAL_WARN("Ignoring multiple results for plugin data file: "
                            <<fullname)
              continue;
            }
            seen.insert(std::move(key));
            be.name = std::move(fullname);
            out.push_back( std::move( be ) );
          }
        }
        return out;
      }
    };
  }
}


void NCD::enableAbsolutePaths( bool doEnable )
{
  s_was_called_enableAbsolutePaths.store(true);
  static std::atomic<bool> s_enabled = {false};
  if ( !setValue(s_enabled,doEnable).wasChanged )
    return;
  if ( !doEnable ) {
    FactImpl::removeTextDataFactoryIfExists( factNameAbsPath );
  } else {
    FactImpl::registerFactory( std::make_unique<TDFact_AbsPath>() );
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace DataSources {

    class TDFact_RelPath final : public FactImpl::TextDataFactory {
      std::string resolve( const TextDataPath& p ) const {
        if ( path_is_absolute( p.path() ) )
          return {};
        if ( file_exists(p.path()) )
          return p.path();
        return {};
      }
    public:
      const char * name() const noexcept override
      {
        return factNameRelPath;
      }
      Priority query( const TextDataPath& p ) const override
      {
        return resolve(p).empty() ? Priority::Unable : Priority{default_priority_relpath};
      }
      TextDataSource produce( const TextDataPath& p ) const override
      {
        std::string s = resolve(p);
        if (s.empty())
          NCRYSTAL_THROW2(DataLoadError,"File disappeared suddenly during request: "<<p.path());
        return TextDataSource::createFromOnDiskPath( s );
      }
      std::vector<BrowseEntry> browse() const override
      {
        return browseDir( ncgetcwd(), Priority{default_priority_relpath} );
      }
    };

    class TDFact_DirList : public FactImpl::TextDataFactory {
      VectS m_dirList;
      Priority m_priority;
      std::string m_name;
      std::string resolve( const TextDataPath& p ) const {
        if ( path_is_absolute( p.path() ) )
          return {};
        //NB: Do NOT allow "../" usage here, it would be too weird if you could
        //"step up" of all searched directories. To make the check bullet proof
        //we just disallow ".." anywhere in the path (thus throwing out some
        //potentially legal but highly dubious use-cases).
        if ( contains( p.path(), ".." ) )
          return {};
        for ( auto& dir : m_dirList ) {
          std::string s = path_join( dir, p.path() );
          if ( file_exists(s) )
            return s;
        }
        return {};
      }
    public:
      TDFact_DirList( VectS&& dirs, std::string name, Priority priority )
        : m_dirList(std::move(dirs)),
          m_priority(priority),
          m_name(std::move(name))
      {
      }

      const char * name() const noexcept override
      {
        return m_name.c_str();
      }
      Priority query( const TextDataPath& p ) const override
      {
        return resolve(p).empty() ? Priority::Unable : m_priority;
      }
      TextDataSource produce( const TextDataPath& p ) const override
      {
        std::string s = resolve(p);
        if (s.empty())
          NCRYSTAL_THROW2(DataLoadError,"File disappeared suddenly during request: "<<p.path());
        return TextDataSource::createFromOnDiskPath( s );
      }
      std::vector<BrowseEntry> browse() const override
      {
        std::vector<BrowseEntry> out;
        for ( auto& dir : m_dirList ) {
          auto b = browseDir( dir, m_priority );
          out.insert(out.end(), b.begin(), b.end());
        }
        return out;
      }
    };

    //custom directories (keeping state in shared structure):
    struct CustomDirList {
      std::mutex mtx;
      std::vector<std::pair<Priority,std::string> > dirList;//sorted by priority (higher first).
    };
    CustomDirList& getCustomDirList() { static CustomDirList cdl; return cdl; }

    class TDFact_CustomDirList : public FactImpl::TextDataFactory {
      std::pair<Priority,std::string> resolve( const TextDataPath& p ) const {
        //NB: Do NOT allow absolute paths or "../" usage here, as in TDFact_DirList
        if ( path_is_absolute( p.path() ) || contains( p.path(), ".." ) )
          return { Priority::Unable, {} };
        auto & cdl = getCustomDirList();
        NCRYSTAL_LOCK_GUARD(cdl.mtx);
        for ( auto& e : cdl.dirList ) {
          std::string s = path_join( e.second, p.path() );
          if ( file_exists(s) )
            return {e.first,s};
        }
        return { Priority::Unable, {} };
      }
    public:
      const char * name() const noexcept override { return factNameCustomSearchDirs; }
      Priority query( const TextDataPath& p ) const override
      {
        return resolve(p).first;
      }
      TextDataSource produce( const TextDataPath& p ) const override
      {
        auto s = resolve(p);
        if ( s.second.empty() )
          NCRYSTAL_THROW2(DataLoadError,"File disappeared suddenly during request: "<<p.path());
        return TextDataSource::createFromOnDiskPath( s.second );
      }
      std::vector<BrowseEntry> browse() const override
      {
        std::vector<BrowseEntry> out;
        auto & cdl = getCustomDirList();
        NCRYSTAL_LOCK_GUARD(cdl.mtx);
        for ( auto& e : cdl.dirList ) {
          auto b = browseDir( e.second, e.first );
          out.insert(out.end(), b.begin(), b.end());
        }
        return out;
      }
    };
  }
}

void NCD::enableRelativePaths( bool doEnable )
{
  s_was_called_enableRelativePaths.store(true);
  static std::atomic<bool> s_enabled = {false};
  if ( !setValue(s_enabled,doEnable).wasChanged )
    return;
  if ( !doEnable ) {
    FactImpl::removeTextDataFactoryIfExists( factNameRelPath );
  } else {
    FactImpl::registerFactory( std::make_unique<TDFact_RelPath>() );
  }
}

namespace NCRYSTAL_NAMESPACE {
#ifdef NCRYSTAL_STDCMAKECFG_EMBED_DATA_ON
  //If NCrystal is installed using the standard CMake with -DEMBED_DATA=ON, an
  //autogenerated function NCrystal::AutoGenNCMAT::registerStdNCMAT will have
  //been compiled in, and when invoking it we will trigger a series of calls to
  //registerEmbeddedNCMAT with the virtual filenames and addresses of the
  //embedded data.
  namespace AutoGenNCMAT { void registerStdNCMAT(); }//fwd declared - linked in elsewhere

  namespace DataSources {
    struct StdDataLibInMemDB {
      std::map<std::string,TextDataSource> virtFileMap;
      std::mutex mtx;
    };
    StdDataLibInMemDB& getStdDataLibInMemDB()
    {
      static StdDataLibInMemDB s_db;
      return s_db;
    }
  }
  namespace internal {
    void registerEmbeddedNCMAT( const char* name, const char* static_data )
    {
      //This function is called repeatedly from autogenerated code's
      //registerStdNCMAT fct when we invoke it. For simplicity we lock/unlock the
      //mutex on each call here (there will anyway be no contention):
      auto& db = NCD::getStdDataLibInMemDB();
      NCRYSTAL_LOCK_GUARD(db.mtx);
      nc_map_force_emplace( db.virtFileMap, name,
                            TextDataSource::createFromInMemData( RawStrData( RawStrData::static_data_ptr_t(), static_data) ) );
    }
  }
#else
  //On-disk standard data library:
  Optional<std::string> getStdDataLibDir()
  {
    //If NCRYSTAL_DATADIR environment variable is set and non-empty we always
    //use that:
    const std::string envvar = ncgetenv("DATADIR");
    if ( !envvar.empty() )
      return envvar;
    //Otherwise use the hardwired (non-relocatable) path:
#ifdef NCRYSTAL_DATADIR
#  define NCRYSTAL_str(s) #s
#  define NCRYSTAL_xstr(s) NCRYSTAL_str(s)
    const std::string hardwired{NCRYSTAL_xstr(NCRYSTAL_DATADIR)};
    if (!hardwired.empty())
      return hardwired;
#endif
    //Hmm... we don't know where it is!
    return NullOpt;
  }

#endif

}

void NCD::enableStandardDataLibrary( bool doEnable, Optional<std::string> requested )
{
  s_was_called_enableStandardDataLibrary.store(true);
  if ( requested.has_value() ) {
    auto rp = tryRealPath( requested.value() );
    if ( !rp.empty() )
      requested = std::move(rp);
  }

  if ( requested.has_value() && !doEnable )
    NCRYSTAL_THROW(BadInput, "Do not provide path to enableStandardDataLibrary when the first argument is false");

  const Priority thePriority = Priority{default_priority_stdlib};

  static std::mutex mtx_esdl;
  NCRYSTAL_LOCK_GUARD(mtx_esdl);
  static bool s_enabled = false;
  static Optional<std::string> s_requested = NullOpt;

  if ( s_enabled == doEnable && s_requested == requested )
    return;//no changes

  //Always remove (might need re-init):
  FactImpl::removeTextDataFactoryIfExists( factNameStdLib );
  s_requested.set(requested);
  s_enabled = doEnable;
  if ( !s_enabled )
    return;

#ifdef NCRYSTAL_STDCMAKECFG_EMBED_DATA_ON
  if ( !s_requested.has_value() ) {
    //Register in-mem factory:
    auto& db = getStdDataLibInMemDB();
    {
      static std::mutex mtx_embed;
      NCRYSTAL_LOCK_GUARD(mtx_embed);
      static bool first = true;
      if ( first ) {
        first = false;
        AutoGenNCMAT::registerStdNCMAT();//trigger db.virtFileMap filling
      }
    }
    NCRYSTAL_LOCK_GUARD(db.mtx);
    //Copy entries map:
    decltype(db.virtFileMap) virtFileMapCopy;
    for ( const auto& e :  db.virtFileMap )
      nc_map_force_emplace( virtFileMapCopy, e.first, e.second );
    registerNamedVirtualDataSource( factNameStdLib, std::move(virtFileMapCopy), thePriority );
    return;
  }
  const std::string phys_dir = ( s_requested.has_value()
                                 ? s_requested.value()
                                 : "" );
#else
  static const Optional<std::string> s_std_path_dir = getStdDataLibDir();
  const std::string phys_dir = ( s_requested.has_value()
                                 ? s_requested.value()
                                 : ( s_std_path_dir.has_value() ? s_std_path_dir.value() : std::string() )
                                 );
#endif
  if (!phys_dir.empty())
    FactImpl::registerFactory( std::make_unique<TDFact_DirList>( VectS{phys_dir}, factNameStdLib, thePriority) );
  //silently register nothing if phys_dir is empty.
}

void NCD::addCustomSearchDirectory( std::string dirname, Priority pr )
{
  if ( !pr.canServiceRequest() || pr.needsExplicitRequest() )
    NCRYSTAL_THROW(BadInput,"addCustomSearchDirectory needs ordinary priority value");

  {
    auto rp = tryRealPath( dirname );
    if ( !rp.empty() )
      dirname = std::move(rp);
  }

  auto & cdl = getCustomDirList();
  NCRYSTAL_LOCK_GUARD(cdl.mtx);
  //check if already added:
  bool added = false;
  for ( auto& e : cdl.dirList ) {
    if ( e.second == dirname ) {
      e.first = pr;
      added = true;
    }
  }
  if (!added)
    cdl.dirList.emplace_back(pr,std::move(dirname));

  //Sort, higher priority first:
  if ( cdl.dirList.size()>1 ) {
    using E_t = decltype(cdl.dirList)::value_type;
    std::stable_sort(cdl.dirList.begin(),cdl.dirList.end(),
                     []( const E_t& a, const E_t& b )
                     {
                       //descending sort
                       nc_assert(a.second!=b.second);
                       nc_assert(a.first.priority()>0);
                       nc_assert(b.first.priority()>0);
                       return a.first.priority() > b.first.priority();
                     });
    nc_assert(cdl.dirList.at(0).first.priority() > cdl.dirList.at(1).first.priority() );
  }

  //Make sure relevant factory is present:
  if (!FactImpl::currentlyHasFactory( FactImpl::FactoryType::TextData,
                                      factNameCustomSearchDirs) )
    FactImpl::registerFactory( std::make_unique<TDFact_CustomDirList>() );
}

void NCD::removeCustomSearchDirectories()
{
  auto & cdl = getCustomDirList();
  NCRYSTAL_LOCK_GUARD(cdl.mtx);
  cdl.dirList.clear();
  FactImpl::removeTextDataFactoryIfExists(factNameCustomSearchDirs);
}

void NCD::enablePluginSearchPaths( bool doEnable )
{
  s_was_called_enablePluginSearchPaths.store(true);
  static std::atomic<bool> s_enabled = {false};
  if ( !setValue(s_enabled,doEnable).wasChanged )
    return;
  if ( !doEnable ) {
    FactImpl::removeTextDataFactoryIfExists( factNamePluginSearchPaths );
  } else {
    FactImpl::registerFactory( std::make_unique<TDFact_PluginDirs>() );
  }
}

void NCD::enableStandardSearchPath( bool doEnable )
{
  s_was_called_enableStandardSearchPath.store(true);
  static std::atomic<bool> s_enabled = {false};
  if ( !setValue(s_enabled,doEnable).wasChanged )
    return;
  if ( !doEnable ) {
    FactImpl::removeTextDataFactoryIfExists( factNameStdSearchPath );
  } else {
    //Determine standard search path:
    VectS dirs;
    auto addEntry = [&dirs](std::string e)
    {
      trim(e);
      if ( e.empty() )
        return;
      auto it = std::find(dirs.begin(),dirs.end(),e);
      if ( it == dirs.end() )
        dirs.emplace_back(std::move(e));
    };
    for ( auto& e : split2(ncgetenv("DATA_PATH"),0,':') )
      addEntry(e);
#ifdef NCRYSTAL_DATA_PATH
    const std::string hardwired{NCRYSTAL_xstr(NCRYSTAL_DATA_PATH)};
    for ( auto& e : split2(hardwired,0,':') )
      addEntry( e );
#endif
    FactImpl::registerFactory( std::make_unique<TDFact_DirList>( std::move(dirs), factNameStdSearchPath, Priority{default_priority_stdpath} ) );
  }
}

namespace NCRYSTAL_NAMESPACE {

  namespace DataSources {
    //Data shared by all virtual files (we only actually register a single factory):
    struct VirtFilesSharedData {
      std::mutex mtx;
      std::map<std::string,std::pair<TextDataSource,Priority>> virtualFiles;
    };

    VirtFilesSharedData& virtualFilesSharedData()
    {
      static VirtFilesSharedData sd; return sd;
    }

    class TDFact_VirtualFiles final : public FactImpl::TextDataFactory {
    public:

      const char * name() const noexcept override
      {
        return factNameVirtualFiles;
      }
      Priority query( const TextDataPath& p ) const override
      {
        auto& vfs = virtualFilesSharedData();
        NCRYSTAL_LOCK_GUARD(vfs.mtx);
        auto it = vfs.virtualFiles.find(p.path());
        return it == vfs.virtualFiles.end() ? Priority::Unable : it->second.second;
      }
      TextDataSource produce( const TextDataPath& p ) const override
      {
        auto& vfs = virtualFilesSharedData();
        NCRYSTAL_LOCK_GUARD(vfs.mtx);
        auto it = vfs.virtualFiles.find(p.path());
        if ( it == vfs.virtualFiles.end() )
          NCRYSTAL_THROW2(DataLoadError,"Virtual file disappeared suddenly during request: "<<p.path());
        return it->second.first;
      }
      std::vector<BrowseEntry> browse() const override {
        auto& vfs = virtualFilesSharedData();
        NCRYSTAL_LOCK_GUARD(vfs.mtx);
        std::vector<BrowseEntry> v;
        v.reserve( vfs.virtualFiles.size() );
        const std::string source = name();
        for ( const auto& e : vfs.virtualFiles )
          v.push_back( BrowseEntry{ e.first, source, e.second.second } );
        return v;
      }

    };

    void validateVirtFilename( const std::string& fn )
    {
      if ( fn.empty() )
        NCRYSTAL_THROW2(BadInput,"Empty file names are not allowed");
      std::string fntrimmed(fn);
      trim(fntrimmed);
      if ( fn != fntrimmed || contains( fn, ' ' ) || contains( fn, '\t' )|| contains( fn, '\r' )|| contains( fn, '\n' ) )
        NCRYSTAL_THROW2(BadInput,"White space is not allowed in file names: \""<<fn<<"\"");
      if ( contains(fn, "::" ) )
        NCRYSTAL_THROW2(BadInput,"Double-semicolons, ::, are not allowed in file names: "<<fn);
    }

    void registerVirtualDataSource( const std::string& virtualFilename,
                                    TextDataSource tsd,
                                    Priority priority )
    {
      validateVirtFilename(virtualFilename);
      auto& vfs = virtualFilesSharedData();
      NCRYSTAL_LOCK_GUARD(vfs.mtx);
      const bool was_empty = vfs.virtualFiles.empty();
      nc_map_force_emplace( vfs.virtualFiles, virtualFilename, std::move(tsd), priority );
      if ( was_empty )
        FactImpl::registerFactory(std::make_unique<TDFact_VirtualFiles>());
    }
  }

  class TDFact_VirtualDataSource final : public FactImpl::TextDataFactory {
  public:
    using VirtFileMap = std::map<std::string,TextDataSource>;
    TDFact_VirtualDataSource( std::string name, VirtFileMap&& vf, Priority priority )
      : m_virtFileMap(std::move(vf)), m_name(name), m_priority(priority)
    {
    }
    const char * name() const noexcept override
    {
      return m_name.c_str();
    }
    Priority query( const TextDataPath& p ) const override
    {
      auto it = m_virtFileMap.find(p.path());
      return it == m_virtFileMap.end() ? Priority::Unable : m_priority;
    }
    TextDataSource produce( const TextDataPath& p ) const override
    {
      auto it = m_virtFileMap.find(p.path());
      nc_assert( it != m_virtFileMap.end() );
      return it->second;
    }
    std::vector<BrowseEntry> browse() const override {
      std::vector<BrowseEntry> v;
      v.reserve( m_virtFileMap.size() );
      const std::string source = m_name;
      for ( const auto& e : m_virtFileMap )
        v.push_back( { e.first, source, m_priority } );
      return v;
    }

  private:
    VirtFileMap m_virtFileMap;
    std::string m_name;
    Priority m_priority;
  };


}

void NCD::registerNamedVirtualDataSource( const std::string& factoryName,
                                          std::map<std::string,TextDataSource>&& virtualFiles,
                                          Priority priority )
{
  if ( !priority.canServiceRequest() )
    NCRYSTAL_THROW(BadInput,"Virtual data sources can not be added with Priority::Unable");
  for ( const auto& e : virtualFiles )
    validateVirtFilename(e.first);

  static std::mutex mtx;
  NCRYSTAL_LOCK_GUARD(mtx);

  auto new_fact = std::make_unique<TDFact_VirtualDataSource>( factoryName,
                                                              std::move(virtualFiles),
                                                              priority );
  if ( FactImpl::currentlyHasFactory(FactImpl::FactoryType::TextData,
                                     new_fact->name()) )
    FactImpl::removeTextDataFactoryIfExists( new_fact->name() );
  FactImpl::registerFactory( std::move(new_fact) );
}

void NCD::registerInMemoryFileData( std::string virtualFileName,
                                    std::string&& data,
                                    Priority priority )
{
  registerVirtualDataSource( std::move(virtualFileName),
                             TextDataSource::createFromInMemData( RawStrData(std::move(data)) ),
                             priority );
}

void NCD::registerInMemoryStaticFileData( std::string virtualFileName,
                                          const char* static_data,
                                          Priority priority )
{
  registerVirtualDataSource( std::move(virtualFileName),
                             TextDataSource::createFromInMemData( RawStrData(RawStrData::static_data_ptr_t(),
                                                                             static_data) ),
                             priority );
}

void NCD::registerVirtualFileAlias( std::string virtualFileName,
                                    std::string realAbsPathFileName,
                                    Priority priority )
{
  if ( !priority.canServiceRequest() )
    NCRYSTAL_THROW(BadInput,"Virtual data sources can not be added with Priority::Unable");
  auto rp = tryRealPath( realAbsPathFileName );
  if ( !rp.empty() )
    realAbsPathFileName = std::move(rp);
  registerVirtualDataSource( std::move(virtualFileName),
                             TextDataSource::createFromOnDiskPath( std::move(realAbsPathFileName) ),
                             priority );
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    struct ExtensionsDB
    {
      std::mutex mtx;
      VectS list;
    };
    ExtensionsDB& extensionsDB() { static ExtensionsDB db; return db; }
  }
}

void NCD::addRecognisedFileExtensions( std::string s )
{
  if ( s.empty() )
    return;
  if ( s.at(0) == '.' )
    s = s.substr(1);
  auto& db = extensionsDB();
  NCRYSTAL_LOCK_GUARD(db.mtx);
  if ( std::find( db.list.begin(), db.list.end(), s )  == db.list.end() )
    db.list.push_back(s);
}

NC::VectS NCD::recognisedFileExtensions()
{
  Plugins::ensurePluginsLoaded();
  auto& db = extensionsDB();
  NCRYSTAL_LOCK_GUARD(db.mtx);
  return db.list;
}

std::vector<NCD::FileListEntry> NCD::listAvailableFiles()
{
  auto facts = FactImpl::getTextDataFactoryList();
  std::vector<NCD::FileListEntry> v;
  v.reserve(128);
  for ( const auto& f : FactImpl::getTextDataFactoryList() ) {
    const std::string factName = f->name();
    auto l = f->browse();
    std::stable_sort(l.begin(),l.end(),[](const decltype(l)::value_type& a,const decltype(l)::value_type& b)
    {
      if ( a.priority != b.priority ) {
        auto pa = a.priority;
        auto pb = b.priority;
        if (!pa.canServiceRequest()||!pb.canServiceRequest())
          NCRYSTAL_THROW2(LogicError,"Factory "<<(!pa.canServiceRequest()?a.name:b.name)
                          <<" browse() method returns entries with Priority::Unable");
        uint_fast32_t pav = pa.needsExplicitRequest() ? 0 : pa.priority();
        uint_fast32_t pbv = pb.needsExplicitRequest() ? 0 : pb.priority();
        return pav > pbv;//descending sort
      }
      if ( a.name != b.name )
        return a.name < b.name;
      return a.source < b.source;
    });
    for ( auto& b : l )
      v.push_back(FileListEntry{std::move(b.name),
                                std::move(b.source),
                                factName,
                                b.priority});
  }
  return v;
}

#ifndef NCRYSTAL_DISABLE_STDDATASOURCES
extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_stddatasrc_factory)()
{
  //Only load if not explicitly modified already:
  if (!NC::DataSources::s_was_called_enableAbsolutePaths.load())
    NC::DataSources::enableAbsolutePaths(true);

  if (!NC::DataSources::s_was_called_enableRelativePaths.load())
    NC::DataSources::enableRelativePaths(true);

  if (!NC::DataSources::s_was_called_enableStandardDataLibrary.load())
    NC::DataSources::enableStandardDataLibrary(true);

  if (!NC::DataSources::s_was_called_enableStandardSearchPath.load())
    NC::DataSources::enableStandardSearchPath(true);

  if (!NC::DataSources::s_was_called_enablePluginSearchPaths.load())
    NC::DataSources::enablePluginSearchPaths(true);
}
#endif

void NCD::removeAllDataSources()
{
  enableAbsolutePaths(false);
  enableRelativePaths(false);
  enableStandardDataLibrary(false);
  enableStandardSearchPath(false);
  removeCustomSearchDirectories();
  {
    auto& vfs = virtualFilesSharedData();
    NCRYSTAL_LOCK_GUARD(vfs.mtx);
    vfs.virtualFiles.clear();
  }
  ::NCrystal::clearCaches();
  //NB: We do not clear out the extensionsDB(), as the Info factories are not
  //unloaded and users might add new data sources. If stdlib is embedded we also
  //leave it alone, so it can be enabled again.
}

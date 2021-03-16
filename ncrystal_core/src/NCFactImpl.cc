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
#include "NCrystal/NCPluginMgmt.hh"
#include "NCrystal/internal/NCFactoryUtils.hh"
#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCFileUtils.hh"

namespace NC = NCrystal;
namespace NCF = NCrystal::FactImpl;

namespace NCrystal {

  namespace FactImpl {

    namespace {

      static std::atomic<bool> s_cache_enabled( ! ncgetenv_bool("NOCACHE") );

      static_assert(Priority{Priority::Unable}.canServiceRequest()==false,"");
      static_assert(Priority{Priority::Unable}.needsExplicitRequest()==false,"");
      static_assert(Priority{Priority::Unable}.priority()==0,"");
      static_assert(Priority{Priority::OnlyOnExplicitRequest}.canServiceRequest()==true,"");
      static_assert(Priority{Priority::OnlyOnExplicitRequest}.needsExplicitRequest()==true,"");
      static_assert(Priority{Priority::OnlyOnExplicitRequest}.priority()==0,"");
#if __cplusplus >= 201703L
      static_assert(Priority{17}.canServiceRequest()==true,"");
      static_assert(Priority{17}.needsExplicitRequest()==false,"");
      static_assert(Priority{17}.priority()==17,"");
#endif

      template<class FactDef>
      class FactDB final : public CachedFactoryBase<typename FactDef::key_type,
                                                    typename FactDef::produced_type,
                                                    FactDef::nstrongrefs_kept,
                                                    typename FactDef::TKeyThinner> {
      public:
        //Since base class depends on the template parameter FactDef, we don't
        //automatically inherit all names, so we need a bunch of using statements:
        using base = CachedFactoryBase<typename FactDef::key_type,typename FactDef::produced_type,
                                       FactDef::nstrongrefs_kept,typename FactDef::TKeyThinner>;
        using typename base::key_type;
        using typename base::value_type;
        using typename base::ShPtr;
        using base::cleanup;
        using FactoryClass = typename FactDef::pubfactory_type;
        using FactoryClassUPtr = std::unique_ptr<const FactoryClass>;
        using FactoryClassShPtr = shared_obj<const FactoryClass>;

        std::string keyToString( const key_type& key ) const final
        {
          return key.toString();
        };
        const char* factoryName() const final
        {
          static std::string name = std::string(FactDef::name())+"FactoryDB";
          return name.c_str();
        }

        ShPtr createWithOrWithoutCache(const key_type& key)
        {
          Plugins::ensurePluginsLoaded();
          return s_cache_enabled ? this->create(key) : this->createWithoutCache(key);
        }
      protected:
        ShPtr actualCreate(const key_type& key) const final
        {
          return FactDef::transformTProdRVToShPtr( searchAndCreateTProdRV(key) );
        }

      public:
        typename FactDef::TProdRV searchAndCreateTProdRV(const key_type& key) const
        {
          const bool verbose = getFactoryVerbosity();

          class VerboseOutput {
            bool m_verbose;
            std::chrono::time_point<std::chrono::steady_clock> m_t0;
            const typename FactDef::key_type* m_keyptr;
          public:
            VerboseOutput(const typename FactDef::key_type*keyptr) : m_verbose(getFactoryVerbosity()), m_keyptr(keyptr)
            {
              if (m_verbose) {
                nc_assert_always(m_keyptr!=nullptr);
                std::cout<<"NCrystal::FactImpl requested to create "<<FactDef::name()<<" based on key "<<m_keyptr->toString()<<std::endl;
                m_t0 = std::chrono::steady_clock::now();
              }
            }
            ~VerboseOutput() {
              if (m_verbose) {
                auto t1 = std::chrono::steady_clock::now();
                double dtsec = std::chrono::duration<double,std::milli>(t1-m_t0).count()*0.001;
                std::cout<<"NCrystal::FactImpl creation of "<<FactDef::name()<<" object based on key "
                         <<m_keyptr->toString()<<" took "<<dtsec<<"s"<<std::endl;
              }
            }
          };

          VerboseOutput produceVerboseOutput(&key);

          key.validate();

          //First consider any specific requests:
          auto requested = FactDef::extractRequestedFactoryName(key);
          if ( requested.excluded.count(requested.specific) )
            NCRYSTAL_THROW2(LogicError,"Factory is both requested and excluded: "<<requested.specific);

          //Get our own list of factories. This locks m_dbmutex briefly, but
          //unlocks it again. This is important since calls to query(..) or
          //produce(..) below might internally invoke other creation calls
          //(e.g. when a scatter factory calls createScatter from its produce
          //method). Of course, this could in principle mean that the factory
          //list can get modified between the top-level creation call and the
          //secondary createScatter. There are ways around this issue, but the
          //current implementation is at least thread-safe - and it is anyway
          //asking for trouble to modify the factory list in one thread while
          //concurrently using the factories in another thread!
          auto factories = getFactoryList();

          //For convenience below, create list of non-excluded factories:
          std::vector<const FactoryClass*> db;
          db.reserve(factories.size());
          for ( auto& e : factories ) {
            if ( !requested.excluded.count(e->name()) )
              db.push_back(&*e);
          }

          //Deal with request for specific factory first:
          if (!requested.specific.empty()) {
            for (auto f : db) {
              if (f->name() == requested.specific) {
                auto priority = f->query( key.getUserFactoryKey() );
                if ( !priority.canServiceRequest() ) {
                  FactDef::produceCustomNoFactFoundError( key, requested.specific );
                  NCRYSTAL_THROW2(BadInput,"Requested "<<FactDef::name()<<" factory \""<<requested.specific
                                  <<"\" does not actually have capability to service request: \""<<key.toString()<<"\"");
                }
                if ( verbose )
                  std::cout<<"NCrystal::FactImpl selected factory [specific request] \""<<f->name()
                           <<"\" to create "<<FactDef::name()<<" based on key "<<key.toString()<<std::endl;
                return f->produce(key.getUserFactoryKey());
              }
            }
            FactDef::produceCustomNoSpecificFactAvail( key, requested.specific );
            NCRYSTAL_THROW2(BadInput,"Specific "<<FactDef::name()<<" factory requested which is unavailable: \""
                            <<requested.specific<<"\"");
          }

          //Nothing specific requested, query all non-excluded factories for
          //their priorities:
          const FactoryClass* best = nullptr;
          Priority best_priority{Priority::Unable};
          for (auto f : db) {
            auto priority = f->query(key.getUserFactoryKey());
            const bool unable = ( !priority.canServiceRequest() || priority.needsExplicitRequest() );
            if ( verbose ) {
              std::cout<<"NCrystal::FactImpl "<<FactDef::name()<<" factory \""<<f->name()
                       <<"\" responded to request for \""<< key.toString()<<"\" with priority: ";
              if ( unable ) {
                std::cout << "UNABLE";
                if ( priority.needsExplicitRequest() )
                  std::cout<<" (NeedsExplicitRequest)";
              } else {
                std::cout << priority.priority();
              }
              std::cout<<std::endl;
            }
            if ( unable )
              continue;
            if ( best == nullptr || best_priority.priority() < priority.priority() ) {
              best = f;
              best_priority = priority;
            }
          }
          if ( best == nullptr ) {
            FactDef::produceCustomNoFactFoundError( key );//give possibility to throw custom msg
            NCRYSTAL_THROW2(BadInput,"Could not find factory to service "<<FactDef::name()
                            <<" creation request for \""<<key.toString()<<"\" ("<<db.size()<<" factories considered)");
          }
          if ( verbose )
            std::cout<<"NCrystal::FactImpl selected factory [highest priority] \""<<best->name()
                     <<"\" to create "<<FactDef::name()<<" based on key "<<key.toString()<<std::endl;
          return best->produce(key.getUserFactoryKey());
        }

      private:
        std::vector<FactoryClassShPtr> m_db;
        mutable std::mutex m_dbmutex;//Must lock whenever accessing m_db
      public:
        void removeFactoryIfExists(const std::string& name)
        {
          Plugins::ensurePluginsLoaded();
          NCRYSTAL_LOCK_GUARD(m_dbmutex);//lock while accessing m_db
          auto it(m_db.begin()), itE(m_db.end());
          for ( ; it!=itE; ++it ) {
            if ( (*it)->name() == name ) {
              //Remove!
              for ( ; std::next(it)!=itE; ++it )
                *it = std::move(*std::next(it));
              m_db.pop_back();
              cleanup();//invalidate all caches
            }
          }
        }

        bool hasFactory(const std::string& name)
        {
          Plugins::ensurePluginsLoaded();
          NCRYSTAL_LOCK_GUARD(m_dbmutex);//lock while accessing m_db
          for ( auto& e : m_db )
            if ( name == e->name() )
              return true;
          return false;
        }

        std::vector<FactoryClassShPtr> getFactoryList() const {
          Plugins::ensurePluginsLoaded();
          NCRYSTAL_LOCK_GUARD(m_dbmutex);//lock while accessing m_db
          std::vector<FactoryClassShPtr> v(m_db);
          return v;
        }

        void addFactory(FactoryClassUPtr f, RegPolicy rp)
        {
          nc_assert_always(!!f);
          Plugins::ensurePluginsLoaded();
          std::string newname(f->name());
          NCRYSTAL_LOCK_GUARD(m_dbmutex);//lock while accessing m_db
          bool inserted(false);
          for ( auto& f_existing : m_db ) {
            if ( newname == f_existing->name() ) {
              if ( rp == RegPolicy::ERROR_IF_EXISTS )
                NCRYSTAL_THROW2(CalcError,"Trying to add "<<FactDef::name()<<" factory \""<<newname
                                <<"\"but existing factory with that name is already registered"
                                " and RegPolicy was set to ERROR_IF_EXISTS");
              if ( rp == RegPolicy::IGNORE_IF_EXISTS )
                return;//don't add!
              nc_assert( rp == RegPolicy::OVERRIDE_IF_EXISTS );
              f_existing = std::move( f );
              inserted = true;
              break;
            }
          }
          if ( !inserted )
            m_db.push_back(std::move(f));
          cleanup();//invalidate all caches
        }
      };

      class DBKey_TextDataPath {
        TextDataPath m_path;
      public:
        DBKey_TextDataPath( const TextDataPath& p ) : m_path(p) {}
        const TextDataPath& getUserFactoryKey() const { return m_path; }
        void validate() const {}//Nothing here, the TextDataPath already validated.
        std::string toString() const { return m_path.toString(); }
        bool operator<(const DBKey_TextDataPath&o) const { return m_path < o.m_path; }
      };

      class DBKey_MatCfg {
        MatCfg m_cfg;
      public:
        using userfact_keytype = MatCfg;
        DBKey_MatCfg( const MatCfg& cfg ) : m_cfg(cfg) {}
        const MatCfg& getUserFactoryKey() const { return m_cfg; }
        void validate() const { m_cfg.checkConsistency(); }
        std::string toString() const { return m_cfg.toStrCfg(); }
        bool operator<(const DBKey_MatCfg&o) const { return m_cfg < o.m_cfg; }
        DBKey_MatCfg cloneThinned() const { return m_cfg.cloneThinned(); }
      };

      class DBKey_MatInfoCfg {
        MatInfoCfg m_cfg;
      public:
        using userfact_keytype = MatInfoCfg;
        DBKey_MatInfoCfg( const MatInfoCfg& cfg ) : m_cfg(cfg) {}
        const MatInfoCfg& getUserFactoryKey() const { return m_cfg; }
        void validate() const { m_cfg.checkConsistency(); }
        std::string toString() const { return m_cfg.toStrCfg(); }
        bool operator<(const DBKey_MatInfoCfg&o) const { return m_cfg < o.m_cfg; }
        DBKey_MatInfoCfg cloneThinned() const { return m_cfg.cloneThinned(); }
      };

      template<class TKey = DBKey_MatCfg>
      struct DBKeyThinner {
        //Thinning (DBKey_)MatCfg objects so we don't have strong TextData refs in the cache map keys.
        using key_type = TKey;
        using thinned_key_type = TKey;
        template <class TMap>
        static typename TMap::mapped_type& cacheMapLookup( TMap& map, const key_type& key, Optional<thinned_key_type>& tkey )
        {
          if ( !tkey.has_value() )
            tkey = key.cloneThinned();
          return map[tkey.value()];
        }
        template <class TMap>
        static typename TMap::iterator cacheMapFind( TMap& map, const key_type& key, Optional<thinned_key_type>& tkey )
        {
          if ( !tkey.has_value() )
            tkey = key.cloneThinned();
          return map.find(tkey.value());
        }
      };

      struct FactDefTextData {
        static constexpr const char* name() { return "TextData"; }
        constexpr static unsigned nstrongrefs_kept = 0;//not used anyway
        using key_type = DBKey_TextDataPath;
        using produced_type = TextDataSource;
        using pubfactory_type = TextDataFactory;
        static MatCfg::FactRequested extractRequestedFactoryName( const key_type& key )
        {
          MatCfg::FactRequested fr;
          fr.specific = key.getUserFactoryKey().fact();
          //Special: paths starting with "./" should correspond to a "relpath' request.
          if ( fr.specific.empty() && startswith(key.getUserFactoryKey().path(),"./") )
            fr.specific = "relpath";
          if ( fr.specific.empty() && path_is_absolute(key.getUserFactoryKey().path()) )
            fr.specific = "abspath";
          //fr.exluded = ...;//TODO: not possible to exclude textdata factories by name!
          return fr;
        }
        static void produceCustomNoSpecificFactAvail( const key_type& key, const std::string& fact_requested )
        {
          if ( fact_requested == "abspath" )
            NCRYSTAL_THROW2(FileNotFound,"Can not load absolute file path since absolute path input has been disabled: \""<<key.getUserFactoryKey().path()<<"\"");
          if ( fact_requested == "relpath" )
            NCRYSTAL_THROW2(FileNotFound,"Can not load relative file path since relative path input has been disabled: \""<<key.getUserFactoryKey().path()<<"\"");
          if ( fact_requested == "stdlib" )
            NCRYSTAL_THROW2(FileNotFound,"Requested \"stdlib\" factory but the standard NCrystal Data library is unavailable or disabled (requested: \""<<key.getUserFactoryKey().path()<<"\")");
          if ( fact_requested == "stdpath" )
            NCRYSTAL_THROW2(FileNotFound,"Requested \"stdpath\" factory but the standard NCrystal search path is unavailable or disabled (requested: \""<<key.getUserFactoryKey().path()<<"\")");
          NCRYSTAL_THROW2(FileNotFound,"Requested factory \""<<fact_requested<<"\" is not available (requested: \""<<key.getUserFactoryKey().path()<<"\").");
        }

        static void produceCustomNoFactFoundError( const key_type& key, const std::string& fact_requested = {} )
        {
          if ( fact_requested.empty() )
            NCRYSTAL_THROW2(FileNotFound,"Could not find data: \""<<key.toString()<<"\"");
          if ( fact_requested == "abspath" || fact_requested == "relpath" )
            NCRYSTAL_THROW2(FileNotFound,"No such file: \""<<key.getUserFactoryKey().path()<<"\"");
          NCRYSTAL_THROW2(FileNotFound,"Requested factory \""<<fact_requested<<"\" can not provide data: \""<<key.getUserFactoryKey().path()<<"\".");
        }
        //Does not produce shared objects:
        using TProdRV = produced_type;
        //Not used, just needs to be here:
        static shared_obj<const TProdRV> transformTProdRVToShPtr( TProdRV ) { return optional_shared_obj<const TProdRV>{nullptr}; }
        using TKeyThinner = CFB_Unthinned_t<key_type>;
      };

      struct FactDefInfo {
        static constexpr const char* name() { return "Info"; }
        constexpr static unsigned nstrongrefs_kept = 20;
        using key_type = DBKey_MatInfoCfg;
        using produced_type = Info;
        using pubfactory_type = FactImpl::InfoFactory;
        static MatCfg::FactRequested extractRequestedFactoryName( const key_type& key )
        {
          MatCfg::FactRequested fr;
          fr.specific = key.getUserFactoryKey().get_infofact_name();
          //fr.exluded = ...;//TODO: not currently possible to exclude info factories by name!
          return fr;
        }
        using TKeyThinner = DBKeyThinner<key_type>;//no strong refs to TextData objects
        static void produceCustomNoSpecificFactAvail( const key_type&, const std::string& ) {}
        static void produceCustomNoFactFoundError( const key_type&, const std::string& = {} ) {}
        //produces shared objects directly:
        using TProdRV = shared_obj<const produced_type>;
        static TProdRV transformTProdRVToShPtr( TProdRV o ) { return o; }
      };

      struct FactDefScatter {
        static constexpr const char* name() { return "Scatter"; }
        constexpr static unsigned nstrongrefs_kept = 20;
        using key_type = DBKey_MatCfg;
        using produced_type = ProcImpl::Process;
        using pubfactory_type = FactImpl::ScatterFactory;
        static MatCfg::FactRequested extractRequestedFactoryName( const key_type& key )
        {
          return key.getUserFactoryKey().get_scatfactory_parsed();
        }
        using TKeyThinner = DBKeyThinner<key_type>;//no strong refs to TextData objects
        static void produceCustomNoSpecificFactAvail( const key_type&, const std::string& ) {}
        static void produceCustomNoFactFoundError( const key_type&, const std::string& = {} ) {}
        //produces shared objects directly:
        using TProdRV = shared_obj<const produced_type>;
        static TProdRV transformTProdRVToShPtr( TProdRV o ) { return o; }
      };

      struct FactDefAbsorption {
        static constexpr const char* name() { return "Absorption"; }
        constexpr static unsigned nstrongrefs_kept = 5;
        using key_type = DBKey_MatCfg;
        using produced_type = ProcImpl::Process;
        using pubfactory_type = FactImpl::AbsorptionFactory;
        static MatCfg::FactRequested extractRequestedFactoryName( const key_type& key )
        {
          return key.getUserFactoryKey().get_absnfactory_parsed();
        }
        using TKeyThinner = DBKeyThinner<key_type>;//no strong refs to TextData objects
        static void produceCustomNoSpecificFactAvail( const key_type&, const std::string& ) {}
        static void produceCustomNoFactFoundError( const key_type&, const std::string& ={} ) {}
        //produces shared objects directly:
        using TProdRV = shared_obj<const produced_type>;
        static TProdRV transformTProdRVToShPtr( TProdRV o ) { return o; }
      };

      //The actual global DB instances (avoid global statics and encapsulate in
      //functions which is MT-safe and gives deterministic init. order):
      FactDB<FactDefTextData>& textDataDB() { static FactDB<FactDefTextData> db; return db; }
      FactDB<FactDefInfo>& infoDB() { static FactDB<FactDefInfo> db; return db; }
      FactDB<FactDefScatter>& scatterDB() { static FactDB<FactDefScatter> db; return db; }
      FactDB<FactDefAbsorption>& absorptionDB() { static FactDB<FactDefAbsorption> db; return db; }
    }

    void setCachingEnabled(bool b) {
      if (getFactoryVerbosity())
        std::cout<<"NCrystal::Factory - called setCachingEnabled("<<(b?"true":"false")<<")."<<std::endl;
      s_cache_enabled = b;
      if (!s_cache_enabled) {
        textDataDB().cleanup();
        infoDB().cleanup();
        scatterDB().cleanup();
        absorptionDB().cleanup();
      }
    }

    bool getCachingEnabled() {
      return s_cache_enabled;
    }
  }
}

void NCF::removeTextDataFactoryIfExists( const std::string& name )
{
  return textDataDB().removeFactoryIfExists( name );
}

void NCF::registerFactory( std::unique_ptr<const TextDataFactory> f, RegPolicy rp ) { textDataDB().addFactory(std::move(f),rp); }
void NCF::registerFactory( std::unique_ptr<const NCF::InfoFactory> f, NCF::RegPolicy rp ) { infoDB().addFactory(std::move(f),rp); }
void NCF::registerFactory( std::unique_ptr<const NCF::ScatterFactory> f, NCF::RegPolicy rp ) { scatterDB().addFactory(std::move(f),rp); }
void NCF::registerFactory( std::unique_ptr<const NCF::AbsorptionFactory> f, NCF::RegPolicy rp ) { absorptionDB().addFactory(std::move(f),rp); }

bool NCF::hasFactory( FactoryType ft, const std::string& name )
{
  switch(ft) {
  case FactoryType::TextData:
    return textDataDB().hasFactory(name);
  case FactoryType::Info:
    return infoDB().hasFactory(name);
  case FactoryType::Scatter:
    return scatterDB().hasFactory(name);
  case FactoryType::Absorption:
    return absorptionDB().hasFactory(name);
  };
  nc_assert_always(false);
  return false;//should not happen
}

std::vector<NC::shared_obj<const NCF::TextDataFactory>> NCF::getTextDataFactoryList() { return textDataDB().getFactoryList(); }
std::vector<NC::shared_obj<const NCF::InfoFactory>> NCF::getInfoFactoryList() { return infoDB().getFactoryList(); }
std::vector<NC::shared_obj<const NCF::ScatterFactory>> NCF::getScatterFactoryList() { return scatterDB().getFactoryList(); }
std::vector<NC::shared_obj<const NCF::AbsorptionFactory>> NCF::getAbsorptionFactoryList() { return absorptionDB().getFactoryList(); }

namespace NCrystal {
  namespace FactImpl {
    //Fwd declare fct implemented in NCTDProd.cc:
   TextDataSP produceTextDataSP_PreferPreviousObject( const TextDataPath&, TextDataSource&& );
  }
}

NC::shared_obj<const NC::TextData> NCF::createTextData( const TextDataPath& path )
{
  //Always recheck the source without cache (file might have changed on-disk,
  //process might have changed working directory, ...):
  auto textDataSource = textDataDB().searchAndCreateTProdRV( path );

  //But whenever identical (and possible given memory constraints of caching),
  //the object is the same as the one returned on previous calls:
  return produceTextDataSP_PreferPreviousObject( path, std::move(textDataSource) );
}

NC::shared_obj<const NC::Info> NCF::createInfo( const MatCfg& cfg )
{
  return infoDB().createWithOrWithoutCache( { cfg.createInfoCfg() } );
}

NC::shared_obj<const NC::ProcImpl::Process> NCF::createScatter( const MatCfg& cfg )
{
  auto p = scatterDB().createWithOrWithoutCache( { cfg } );
  auto pt = p->processType();
  if ( pt != ProcessType::Scatter )
    NCRYSTAL_THROW2(CalcError,"Scatter factory created "<<pt<<" process!");
  return p;
}

NC::shared_obj<const NC::ProcImpl::Process> NCF::createAbsorption( const MatCfg& cfg )
{
  auto p = absorptionDB().createWithOrWithoutCache( { cfg } );
  auto pt = p->processType();
  if ( pt != ProcessType::Absorption )
    NCRYSTAL_THROW2(CalcError,"Absorption factory created "<<pt<<" process!");
  return p;
}

NC::ProcImpl::ProcPtr NCF::ScatterFactory::globalCreateScatter( const MatCfg& cfg, bool allowself ) const
{
  auto cfg2 = cfg.clone();
  if ( ! allowself ) {
    std::string ourname=name();
    auto factrequest = cfg2.get_scatfactory_parsed();
    factrequest.excluded.insert(ourname);
    if (factrequest.specific==ourname)
      factrequest.specific.clear();
    cfg2.set_scatfactory(factrequest);
  }
  return ::NCrystal::FactImpl::createScatter(cfg2);
}

NC::shared_obj<const NC::Info> NCF::ScatterFactory::globalCreateInfo( const MatCfg& cfg )
{
  return ::NCrystal::FactImpl::createInfo(cfg);
}

NC::shared_obj<const NC::Info> NCF::ScatterFactory::createInfo( const MatCfg& cfg )
{
  return ::NCrystal::FactImpl::createInfo(cfg);
}

std::string NCF::guessDataType( const RawStrData& data,
                                const std::string& filename )
{
  //Figure out data type. We are able to recognise NCMAT content from the
  //data itself, and otherwise we look at the file extension.
  if ( 0 == std::strncmp( data.begin(), "NCMAT", 5 ) )
    return "ncmat"_s;
  auto ext = getfileext(filename);
  if ( !ext.empty() && isAlphaNumeric(ext) )
    return lowerCase(ext);
  return {};
}

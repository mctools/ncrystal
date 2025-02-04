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
#include "NCrystal/plugins/NCPluginMgmt.hh"
#include "NCrystal/internal/fact_utils/NCFactoryUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCFileUtils.hh"
#include "NCrystal/internal/infobld/NCInfoBuilder.hh"
#include "NCrystal/internal/cfgutils/NCCfgManip.hh"
#include "NCrystal/internal/extd_utils/NCProcCompBldr.hh"
#include "NCrystal/threads/NCFactThreads.hh"
#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include <list>

namespace NC = NCrystal;
namespace NCF = NCrystal::FactImpl;

namespace NCRYSTAL_NAMESPACE {

  namespace FactImpl {

    namespace {

      static_assert(Priority{Priority::Unable}.canServiceRequest()==false,"");
      static_assert(Priority{Priority::Unable}.needsExplicitRequest()==false,"");
      static_assert(Priority{Priority::Unable}.priority()==0,"");
      static_assert(Priority{Priority::OnlyOnExplicitRequest}.canServiceRequest()==true,"");
      static_assert(Priority{Priority::OnlyOnExplicitRequest}.needsExplicitRequest()==true,"");
      static_assert(Priority{Priority::OnlyOnExplicitRequest}.priority()==0,"");
#if nc_cplusplus >= 201703L
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

        ShPtr loadPluginsAndCreate(const key_type& key)
        {
          Plugins::ensurePluginsLoaded();
          return this->create(key);
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
                NCRYSTAL_MSG("FactImpl requested to create "<<FactDef::name()<<" based on key "<<m_keyptr->toString());
                m_t0 = std::chrono::steady_clock::now();
              }
            }
            ~VerboseOutput() {
              if (m_verbose) {
                auto t1 = std::chrono::steady_clock::now();
                double dtsec = std::chrono::duration<double,std::milli>(t1-m_t0).count()*0.001;
                NCRYSTAL_MSG("FactImpl creation of "<<FactDef::name()<<" object based on key "
                             <<m_keyptr->toString()<<" took "<<dtsec<<"s");
              }
            }
          };

          VerboseOutput produceVerboseOutput(&key);

          //First consider any specific requests:
          auto requested = FactDef::extractRequestedFactory(key);
          nc_assert_always( !requested.hasSpecificRequest()
                            || !requested.excludes(requested.specificRequest()) );//already checked by FactNameRequest::Parser

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
            if ( !requested.excludes(e->name()) )
              db.push_back(&*e);
          }

          //Deal with request for specific factory first:
          if ( requested.hasSpecificRequest() ) {
            for (auto f : db) {
              if (f->name() == requested.specificRequest()) {
                auto priority = FactDef::isSuitableForKey(*f,key) ? f->query( key.getUserFactoryKey() ) : Priority::Unable;
                if ( !priority.canServiceRequest() ) {
                  FactDef::produceCustomNoFactFoundError( key, requested.specificRequest() );
                  NCRYSTAL_THROW2(BadInput,"Requested "<<FactDef::name()<<" factory \""<<requested.specificRequest()
                                  <<"\" does not actually have capability to service request: \""<<key.toString()<<"\"");
                }
                if ( verbose )
                  NCRYSTAL_MSG("FactImpl selected factory [specific request] \""<<f->name()
                               <<"\" to create "<<FactDef::name()<<" based on key "<<key.toString());
                return f->produce(key.getUserFactoryKey());
              }
            }
            FactDef::produceCustomNoSpecificFactAvail( key, requested.specificRequest() );
            NCRYSTAL_THROW2(BadInput,"Specific "<<FactDef::name()<<" factory requested which is unavailable: \""
                            <<requested.specificRequest()<<"\"");
          }

          //Nothing specific requested, query all non-excluded factories for
          //their priorities:
          const FactoryClass* best = nullptr;
          Priority best_priority{Priority::Unable};
          for (auto f : db) {
            auto priority = FactDef::isSuitableForKey(*f,key) ? f->query( key.getUserFactoryKey() ) : Priority::Unable;
            const bool unable = ( !priority.canServiceRequest() || priority.needsExplicitRequest() );
            if ( verbose ) {
              std::ostringstream ss;
              ss<<"FactImpl "<<FactDef::name()<<" factory \""<<f->name()
                       <<"\" responded to request for \""<< key.toString()<<"\" with priority: ";
              if ( unable ) {
                ss << "UNABLE";
                if ( priority.needsExplicitRequest() )
                  ss<<" (NeedsExplicitRequest)";
              } else {
                ss << priority.priority();
              }
              Msg::outputMsg(ss.str());
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
            NCRYSTAL_MSG("FactImpl selected factory [highest priority] \""<<best->name()
                         <<"\" to create "<<FactDef::name()<<" based on key "<<key.toString());
          return best->produce(key.getUserFactoryKey());
        }

      private:
        std::vector<FactoryClassShPtr> m_db;
        mutable std::mutex m_dbmutex;//Must lock whenever accessing m_db
      public:
        void removeFactoryIfExists(const std::string& name)
        {
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

        void addFactory(FactoryClassUPtr f)
        {
          nc_assert_always(!!f);
          std::string newname(f->name());
          NCRYSTAL_LOCK_GUARD(m_dbmutex);//lock while accessing m_db
          bool inserted(false);
          for ( auto& f_existing : m_db ) {
            if ( newname == f_existing->name() ) {
              NCRYSTAL_THROW2(CalcError,"Trying to add "<<FactDef::name()
                              <<" factory \""<<newname <<"\"but existing"
                              " factory with that name is already registered");
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
        std::string toString() const { return m_path.toString(); }
        bool operator<(const DBKey_TextDataPath&o) const { return m_path < o.m_path; }
      };

      template<class TXXXRequest>
      class DBKey_XXXRequest {
        TXXXRequest m_cfg;
      public:
        using userfact_keytype = TXXXRequest;
        DBKey_XXXRequest( const TXXXRequest& cfg ) : m_cfg(cfg) {}
        const TXXXRequest& getUserFactoryKey() const { return m_cfg; }
        std::string toString() const { std::ostringstream ss; m_cfg.stream(ss); return ss.str(); }
        bool operator<(const DBKey_XXXRequest&o) const { return m_cfg < o.m_cfg; }
        DBKey_XXXRequest cloneThinned() const { return m_cfg.cloneThinned(); }
      };

      using DBKey_InfoRequest = DBKey_XXXRequest<InfoRequest>;
      using DBKey_ScatterRequest = DBKey_XXXRequest<ScatterRequest>;
      using DBKey_AbsorptionRequest = DBKey_XXXRequest<AbsorptionRequest>;

      template<class TKey>
      struct DBKeyThinner {
        //Thinning key objects so we don't have strong TextData/Info refs in the cache map keys.
        using key_type = TKey;
        using thinned_key_type = TKey;
        template <class TMap>
        static typename TMap::mapped_type& cacheMapLookup( TMap& map, const key_type& key, Optional<thinned_key_type>& tkey )
        {
          if ( !tkey.has_value() )
            tkey = key.cloneThinned();
          return map[tkey.value()];
        }
      };

      struct FactDefTextData {
        static constexpr const char* name() { return "TextData"; }
        constexpr static unsigned nstrongrefs_kept = 0;//not used anyway
        using key_type = DBKey_TextDataPath;
        using produced_type = TextDataSource;
        using pubfactory_type = TextDataFactory;
        static bool isSuitableForKey( const pubfactory_type&, const key_type& ) { return true; }
        static Cfg::FactNameRequest extractRequestedFactory( const key_type& key )
        {
          std::string specific = key.getUserFactoryKey().fact();
          //Special: paths starting with "./" should correspond to a "relpath' request.
          if ( specific.empty() && startswith(key.getUserFactoryKey().path(),"./") )
            specific = "relpath";
          if ( specific.empty() && path_is_absolute(key.getUserFactoryKey().path()) )
            specific = "abspath";
          //fr.exluded = ...;//TODO: not possible to exclude textdata factories by name!
          return Cfg::FactNameRequest::Parser::doParse(specific);
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
        using key_type = DBKey_InfoRequest;
        using produced_type = Info;
        using pubfactory_type = FactImpl::InfoFactory;
        static bool isSuitableForKey( const pubfactory_type&, const key_type& ) { return true; }
        static Cfg::FactNameRequest extractRequestedFactory( const key_type& key )
        {
          auto fnrstr = key.getUserFactoryKey().get_infofactory();
          return Cfg::FactNameRequest::Parser::doParse(fnrstr);
        }
        using TKeyThinner = DBKeyThinner<key_type>;//no strong refs to TextData objects
        static void produceCustomNoSpecificFactAvail( const key_type&, const std::string& ) {}
        static void produceCustomNoFactFoundError( const key_type&, const std::string& = {} ) {}
        //produces shared objects directly:
        using TProdRV = shared_obj<const produced_type>;
        static TProdRV transformTProdRVToShPtr( TProdRV o ) { return o; }
      };

      template<class TXXXRequest>
      bool singleMultiPhaseSuitability(MultiPhaseCapability cap, const TXXXRequest& cfg )
      {
        //Info-lvl multiphase
        return ( cfg.isMultiPhase()
                 ? ( cap==MultiPhaseCapability::MPOnly || cap == MultiPhaseCapability::Both )
                 : ( cap==MultiPhaseCapability::SPOnly || cap == MultiPhaseCapability::Both ) );
      }

      struct FactDefScatter {
        static constexpr const char* name() { return "Scatter"; }
        constexpr static unsigned nstrongrefs_kept = 20;
        using key_type = DBKey_ScatterRequest;
        using produced_type = ProcImpl::Process;
        using pubfactory_type = FactImpl::ScatterFactory;
        static bool isSuitableForKey( const pubfactory_type& pf, const key_type& key )
        {
          return singleMultiPhaseSuitability(pf.multiPhaseCapability(),key.getUserFactoryKey());
        }
        static Cfg::FactNameRequest extractRequestedFactory( const key_type& key )
        {
          auto fnrstr = key.getUserFactoryKey().get_scatfactory();
          return Cfg::FactNameRequest::Parser::doParse(fnrstr);
        }

        using TKeyThinner = DBKeyThinner<key_type>;
        static void produceCustomNoSpecificFactAvail( const key_type&, const std::string& ) {}
        static void produceCustomNoFactFoundError( const key_type&, const std::string& = {} ) {}
        //produces shared objects directly:
        using TProdRV = shared_obj<const produced_type>;
        static TProdRV transformTProdRVToShPtr( TProdRV o ) { return o; }
      };

      struct FactDefAbsorption {
        static constexpr const char* name() { return "Absorption"; }
        constexpr static unsigned nstrongrefs_kept = 5;
        using key_type = DBKey_AbsorptionRequest;
        using produced_type = ProcImpl::Process;
        using pubfactory_type = FactImpl::AbsorptionFactory;
        static bool isSuitableForKey( const pubfactory_type& pf, const key_type& key )
        {
          return singleMultiPhaseSuitability(pf.multiPhaseCapability(),key.getUserFactoryKey());
        }
        static Cfg::FactNameRequest extractRequestedFactory( const key_type& key )
        {
          auto fnrstr = key.getUserFactoryKey().get_absnfactory();
          return Cfg::FactNameRequest::Parser::doParse(fnrstr);
        }
        using TKeyThinner = DBKeyThinner<key_type>;
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
  }
}

void NCF::removeTextDataFactoryIfExists( const std::string& name )
{
  return textDataDB().removeFactoryIfExists( name );
}

void NCF::registerFactory( std::unique_ptr<const TextDataFactory> f ) { textDataDB().addFactory(std::move(f)); }
void NCF::registerFactory( std::unique_ptr<const NCF::InfoFactory> f ) { infoDB().addFactory(std::move(f)); }
void NCF::registerFactory( std::unique_ptr<const NCF::ScatterFactory> f ) { scatterDB().addFactory(std::move(f)); }
void NCF::registerFactory( std::unique_ptr<const NCF::AbsorptionFactory> f ) { absorptionDB().addFactory(std::move(f)); }

bool NCF::currentlyHasFactory( FactoryType ft, const std::string& name )
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

bool NCF::hasFactory( FactoryType ft, const std::string& name )
{
  Plugins::ensurePluginsLoaded();
  return currentlyHasFactory(ft, name);
}

std::vector<NC::shared_obj<const NCF::TextDataFactory>> NCF::getTextDataFactoryList() { return textDataDB().getFactoryList(); }
std::vector<NC::shared_obj<const NCF::InfoFactory>> NCF::getInfoFactoryList() { return infoDB().getFactoryList(); }
std::vector<NC::shared_obj<const NCF::ScatterFactory>> NCF::getScatterFactoryList() { return scatterDB().getFactoryList(); }
std::vector<NC::shared_obj<const NCF::AbsorptionFactory>> NCF::getAbsorptionFactoryList() { return absorptionDB().getFactoryList(); }

namespace NCRYSTAL_NAMESPACE {
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

namespace NCRYSTAL_NAMESPACE {

  namespace FactImpl {

    namespace {

      template<typename TRequest>
      class CfgLvlMPProc_Key {
      public:
        using value_type = std::pair<double,TRequest>;

        CfgLvlMPProc_Key( const MatCfg::PhaseList& cfg_phases )
        {
          //Must prepare TRequest objects for each phase, but must also at this
          //point calculate the relative contribution of each phase to the number
          //density, as this contribution determines the weight factor of the
          //associated physics processes. This information can not be derived from
          //the info objects in the TRequest objects, as that object will be
          //independent of the MatCfg level density state.
          //
          //Also we are careful to NOT pass non-trivial (e.g. with density
          //state) cfg objects to TRequest constructors.
          nc_assert(cfg_phases.size()>=2);
          m_data.reserve(cfg_phases.size());
          StableSum combinedNumberDensity;
          for ( auto& e : cfg_phases ) {
            const double volfrac = e.first;
            const auto& cfg = e.second;
            nc_assert_always(!cfg.isMultiPhase());
            nc_assert( cfg.isSinglePhase() );
            //Always get number density of phase like this, to handle density
            //overrides, phase-choices, etc.:
            auto info = FactImpl::createInfo(cfg);
            const double numdens_contrib = volfrac * info->getNumberDensity().dbl();
            m_data.emplace_back( numdens_contrib,//<--will be normalised below
                                 TRequest{info} );
            nc_assert(!m_data.back().second.isThinned());
            combinedNumberDensity.add( numdens_contrib );
          }
          const double totnd = combinedNumberDensity.sum();
          if ( !(totnd>0.0) ) {
            //shouldn't happen, but clear m_data to signal production of null process
            m_data.clear();
          } else {
            for ( auto& e : m_data )
              e.first /= totnd;
          }
          //For consistency (and to make phase specification order irrelevant
          //for process objects), sort. Note that no sorting should happen
          //BEFORE phasechoices are evaluated, since that would mess with user
          //expectations:
          std::stable_sort( m_data.begin(), m_data.end() );
        }

        const value_type* begin() const { return &m_data[0]; }
        const value_type* end() const { return &m_data[0] + m_data.size(); }
        bool empty() const { return m_data.empty(); }

        bool isThinned() const
        {
          return m_data.at(0).isThinned();
        }

        CfgLvlMPProc_Key cloneThinned() const
        {
          CfgLvlMPProc_Key res;
          res.m_data.reserve(m_data.size());
          for ( auto& e : m_data )
            res.m_data.emplace_back( e.first, e.second.cloneThinned() );
          return res;
        }

        bool operator<(const CfgLvlMPProc_Key& o ) const
        {
          auto n = m_data.size();
          if ( n != o.m_data.size() )
            return n < o.m_data.size();
          for ( auto i : ncrange(n) ) {
            nc_assert_always( !std::isnan( m_data[i].first ) );
            nc_assert_always( !std::isnan( o.m_data[i].first ) );
            if ( m_data[i].first != o.m_data[i].first )
              return m_data[i].first < o.m_data[i].first;
          }
          for ( auto i : ncrange(n) ) {
            nc_assert( (m_data[i].second == o.m_data[i].second ) == ( !(m_data[i].second < o.m_data[i].second) && !(o.m_data[i].second < m_data[i].second)) );
            if ( ! (m_data[i].second == o.m_data[i].second ) )
              return m_data[i].second < o.m_data[i].second;
          }
          //equal:
          return false;
        }
      private:
        CfgLvlMPProc_Key() = default;
        std::vector<value_type> m_data;
      };

      template<class TRequest, unsigned nstrongrefs_kept = 20>
      class MPProcCacheDB : public CachedFactoryBase<CfgLvlMPProc_Key<TRequest>,
                                                     ProcImpl::Process,
                                                     nstrongrefs_kept,
                                                     FactImpl::DBKeyThinner<CfgLvlMPProc_Key<TRequest>>> {
      public:
        using key_type = CfgLvlMPProc_Key<TRequest>;
        std::string keyToString( const key_type& key ) const final
        {
          //NB: No attempt to print common parameters at the end - it would be overkill for this usage.
          std::ostringstream ss;
          ss << "MPProcRequest<";
          bool first = true;
          for ( auto& e : key ) {
            if (!first)
              ss<<'&';
            first = false;
            ss<<e.first<<'*'<<e.second;
          }
          ss<<'>';
          return ss.str();
        }
        const char* factoryName() const final
        {
          static_assert( std::is_same<TRequest,ScatterRequest>::value
                         || std::is_same<TRequest,AbsorptionRequest>::value,"" );
          return ( std::is_same<TRequest,ScatterRequest>::value
                   ? "CfgLvlMultiPhaseScatterBuilder"
                   : "CfgLvlMultiPhaseAbsorptionBuilder" );
        }
      protected:
        std::shared_ptr<const ProcImpl::Process> actualCreate(const key_type& key) const final
        {
          static_assert( std::is_same<TRequest,ScatterRequest>::value
                         || std::is_same<TRequest,AbsorptionRequest>::value,"" );
          constexpr ProcessType processType = ( std::is_same<TRequest,ScatterRequest>::value
                                                ? ProcessType::Scatter
                                                : ProcessType::Absorption );
          if ( key.empty() )
            return ProcImpl::getGlobalNullProcess(processType);
#if 0
          ProcImpl::ProcComposition::ComponentList proclist;
          for ( auto& e : key )
            proclist.push_back( ProcImpl::ProcComposition::Component{ e.first, FactImpl::create( e.second ) } );
          return ProcImpl::ProcComposition::consumeAndCombine( std::move(proclist), processType );
#else
          Utils::ProcCompBldr proclist;
          //FTUtils::ComponentAggregator proclist;
          for ( auto& scale_and_request : key ) {
            proclist.addfct_cl( [scale_and_request]() {
              ProcImpl::ProcComposition::ComponentList cl;
              cl.emplace_back( scale_and_request.first,
                               FactImpl::create( scale_and_request.second ) );
              return cl;
            });
          }
          return ProcImpl::ProcComposition::consumeAndCombine( proclist.finalise(), processType );
#endif
        }
      };

      template<ProcessType processType, class TRequest>
      ProcImpl::ProcPtr createProcFromMPCfg(const MatCfg::PhaseList& cfg_phases)
      {
        nc_assert_always(cfg_phases.size()>=2);
        //Create via cache:
        CfgLvlMPProc_Key<TRequest> key(cfg_phases);
        static MPProcCacheDB<TRequest> s_db;
        return s_db.create(key);
      }

      //Mini-factory for the creation of Info objects from MatCfg objects with
      //isMultiPhase()=true. These can not be created from the other Info
      //factory infrastructure since that one works on InfoRequests which are
      //always based on single-phase objects (at the cfg-level).

      struct MultiPhaseMatCfgCache {
        std::mutex mtx;
        std::map<MatCfg,std::weak_ptr<const Info>> db;
        //Strong refs:
        std::list<InfoPtr> strong_refs;
        void cropStrongRefs() {
          //Expire the least accessed strong refs if
          //we are holding too many:
          while ( strong_refs.size() > 20 )
            strong_refs.pop_front();
        };
      };
      MultiPhaseMatCfgCache& getMultiPhaseMatCfgCache()
      {
        static MultiPhaseMatCfgCache db;
        return db;
      }
      void clearMPCfgInfoCache() {
        auto& cache = getMultiPhaseMatCfgCache();
        NCRYSTAL_LOCK_GUARD(cache.mtx);
        cache.db.clear();
        cache.strong_refs.clear();
      }
      InfoPtr createInfoFromMultiPhaseCfg(const MatCfg& cfg)
      {
        const bool verbose = getFactoryVerbosity();

        nc_assert( cfg.isMultiPhase() && cfg.phases().size()>1 );
        auto& cache = getMultiPhaseMatCfgCache();
        {
          //Acquire lock and check cache:
          NCRYSTAL_LOCK_GUARD(cache.mtx);
          auto it = cache.db.find(cfg);
          if ( it != cache.db.end() ) {
            auto res = it->second.lock();
            if ( res != nullptr ) {
              //Great, it was already in the cache! Only thing left is to make
              //sure it is placed at the back of the strong_refs list.
              auto& SR = cache.strong_refs;
              if ( SR.empty() || SR.back() != res ) {
                auto itsr = std::find(SR.begin(), SR.end(), res);
                if ( itsr == SR.end() ) {
                  SR.push_back(res);
                  cache.cropStrongRefs();
                } else {
                  //Move existing entry to the end:
                  SR.splice( SR.end(), SR, itsr );
                }
              }
              if ( verbose )
                NCRYSTAL_MSG("FactImpl (thread_"<<thread_details::currentThreadIDForPrint()
                             <<") Returning existing MatCfg-level multiphase Info object from key "<<cfg);
              return res;
            }
            //weakptr apparently expired:
            cache.db.erase(it);
          }
          //Good place to register cache cleanup fct:
          {
            static bool first = true;
            if ( first ) {
              first = false;
              //We register to be cleaned up whenever the Info factory is cleaned
              //up:
              FactImpl::infoDB().registerCleanupCallback(clearMPCfgInfoCache);
              registerCacheCleanupFunction(clearMPCfgInfoCache);
            }
          }
        }

        //Lock released again.

        //Nothing in the cache, build the object from scratch:
        if ( verbose )
          NCRYSTAL_MSG("FactImpl (thread_"<<thread_details::currentThreadIDForPrint()
                       <<") creating (from scratch) MatCfg-level multiphase Info object from key "<<cfg);

        InfoBuilder::MultiPhaseBuilder mp_builder;
        const auto& cfg_phases = cfg.phases();
        mp_builder.phases.reserve(cfg_phases.size());
#if 0
        //old:
        for ( auto& cfg_ph : cfg_phases )
          mp_builder.phases.emplace_back(cfg_ph.first,FactImpl::createInfo(cfg_ph.second));
#else
        //MT:
        {
          FactoryJobs jobs;
          SmallVector<OptionalInfoPtr,6> tmp_phases;
          tmp_phases.resize(cfg_phases.size());
          unsigned iphase = 0;
          for ( auto& cfg_ph : cfg_phases ) {
            OptionalInfoPtr * tgt = &tmp_phases[iphase++];
            jobs.queue([tgt,&cfg_ph]()
            {
              *tgt = FactImpl::createInfo(cfg_ph.second);
            });
          }
          jobs.waitAll();
          iphase = 0;
          for ( auto& cfg_ph : cfg_phases )
            mp_builder.phases.emplace_back(cfg_ph.first,std::move(tmp_phases.at(iphase++)));
        }
#endif
        auto res = InfoBuilder::buildInfoPtr( std::move(mp_builder) );

        {
          //Acquire lock again and check cache again:
          NCRYSTAL_LOCK_GUARD(cache.mtx);
          auto& entry_wp = cache.db[cfg];
          auto entry_sp = entry_wp.lock();
          if ( entry_sp != nullptr ) {
            //Looks like another thread beat us to the creation. We simply return
            //that one and discard the one we just created (we don't bother
            //updating strong_ref here since the other thread will have done so
            //recently).
            if ( verbose )
              NCRYSTAL_MSG("FactImpl (thread_"<<thread_details::currentThreadIDForPrint()
                           <<") Discarding MatCfg-level multiphase Info object from key "
                           <<cfg<<" (competing thread beat us to it)");
            return entry_sp;
          }
          //Ok, we have to return the one we created. First update cache:
          entry_wp = res;
          //Update strong refs:
          cache.strong_refs.push_back(res);
          cache.cropStrongRefs();
          //Return:
          return res;
        }
      }

    }
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace detail {
#ifndef NCRYSTAL_DISABLE_THREADS
    bool factThreadsEnableCalledExplicitly();//fwd declare fct from NCFactThreads.cc

    namespace {
      void factThreads_checkEnvVar()
      {
        static std::atomic<bool> first(true);
        bool btrue(true);
        if ( !first.compare_exchange_strong(btrue,false) )
          return;
        std::int64_t nthreads_raw = ncgetenv_int64("FACTORY_THREADS",-1);
        if ( nthreads_raw >= 0 && !factThreadsEnableCalledExplicitly()) {
          auto nthreads = ThreadCount{ nthreads_raw > 9999
                                       ? 9999
                                       : static_cast<unsigned>(nthreads_raw) };
          FactoryThreadPool::enable( nthreads );
        }
      }
    }
#else
    namespace { void factThreads_checkEnvVar() {} }
#endif
  }
}

NC::shared_obj<const NC::Info> NCF::createInfo( const InfoRequest& cfg )
{
  ::NC::detail::factThreads_checkEnvVar();
  return infoDB().loadPluginsAndCreate( { cfg } );
}

NC::ProcImpl::ProcPtr NCF::createScatter( const ScatterRequest& cfg )
{
  ::NC::detail::factThreads_checkEnvVar();
  auto p = scatterDB().loadPluginsAndCreate( { cfg } );
  nc_assert( p!=nullptr );
  if ( p->processType() != ProcessType::Scatter )
    NCRYSTAL_THROW2(CalcError,"Scatter factory created "<<p->processType()<<" process!");
  //All null processes use same global instances:
  return p->isNull() ? ProcImpl::getGlobalNullProcess(p->processType()) : ProcImpl::ProcPtr{p};
}

NC::ProcImpl::ProcPtr NCF::createAbsorption( const AbsorptionRequest& cfg )
{
  ::NC::detail::factThreads_checkEnvVar();
  auto p = absorptionDB().loadPluginsAndCreate( { cfg } );
  nc_assert( p!=nullptr );
  if ( p->processType() != ProcessType::Absorption )
    NCRYSTAL_THROW2(CalcError,"Absorption factory created "<<p->processType()<<" process!");
  //All null processes use same global instances:
  return p->isNull() ? ProcImpl::getGlobalNullProcess(p->processType()) : ProcImpl::ProcPtr{p};
}

NC::shared_obj<const NC::Info> NCF::createInfo( const MatCfg& cfg )
{
  ::NC::detail::factThreads_checkEnvVar();
  ///////////////////////////////////////////////////////////////////////////////////////////////
  // First deal with phase-choices (before all other things, which is in line
  // with the documentation's promise that the effect of phase-choice parameter
  // is applied after Info objects are loaded, based on all other parameters),
  // making sure that only cfg objects without top-level phase-choices gets
  // looked up below (to avoid populating the cache with essentially identical
  // objects):
  MatCfg::PhaseChoices phaseChoices = cfg.getPhaseChoices();
  if ( !phaseChoices.empty() ) {
    auto info = createInfo( cfg.cloneWithoutPhaseChoices() );
    //Recursively zoom in to specific chosen phase:
    for ( auto iphasechoice : phaseChoices ) {
      if ( ! info->isMultiPhase() || !( iphasechoice < info->getPhases().size() ) )
        NCRYSTAL_THROW(BadInput,"Invalid phase choice.");
      info = info->getPhases().at(iphasechoice).second;
    }
    return info;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // Next deal with special density requests:

  if ( cfg.hasDensityOverride() ) {
    using DType = DensityState::Type;
    auto info_underlying = createInfo( cfg.cloneWithoutDensityState() );
    auto ds = cfg.get_density();
    if ( info_underlying->isSinglePhase() ) {
      //Single phase => Simply override or scale the density:
      if ( ds.type == DType::SCALEFACTOR )
        return InfoBuilder::buildInfoPtrWithScaledDensity( info_underlying,ds.value );
      if ( ds.type == DType::NUMBERDENSITY )
        return InfoBuilder::buildInfoPtr(info_underlying,NumberDensity{ ds.value });
      nc_assert( ds.type == DType::DENSITY );
      return InfoBuilder::buildInfoPtr(info_underlying,Density{ ds.value });
    } else {
      //Multi-phase => figure out density scale, apply it to all phases, and put together in a new object.
      nc_assert( info_underlying->getNumberDensity().dbl() > 0.0 );
      nc_assert( info_underlying->getDensity().dbl() > 0.0 );
      double scale_factor = ds.value;
      if ( ds.type == DType::NUMBERDENSITY )
        scale_factor = ds.value / info_underlying->getNumberDensity().dbl();
      if ( ds.type == DType::DENSITY )
        scale_factor = ds.value / info_underlying->getDensity().dbl();
      InfoBuilder::MultiPhaseBuilder mpbldr;
      mpbldr.phases.reserve(info_underlying->getPhases().size());
      for ( auto& ph : info_underlying->getPhases() ) {
        mpbldr.phases.emplace_back(ph.first,
                                   InfoBuilder::buildInfoPtrWithScaledDensity(ph.second,scale_factor));

      }
      return InfoBuilder::buildInfoPtr( std::move(mpbldr) );
    }
  }
  nc_assert( !cfg.hasDensityOverride() );

  ///////////////////////////////////////////////////////////////////////////////////////////////
  // Deal with multi-phase MatCfg objects:

  nc_assert(cfg.getPhaseChoices().empty());
  if ( cfg.isMultiPhase() ) {
    nc_assert_always( cfg.phases().size() >= 2 );
    return createInfoFromMultiPhaseCfg(cfg);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////
  //Finally cfg is here is a single-phase object without phase choices or
  //density overrides. Request this from the InfoRequest-based infrastructure:

  //NB: It is possible that the created Info object will be a multi-phase
  //object, if e.g. a single .ncmat file provides multiple phases.

  nc_assert(cfg.isTrivial());
  auto info = createInfo( InfoRequest(cfg) );

  return InfoBuilder::recordCfgDataOnInfoObject( std::move(info), cfg.rawCfgData() );
}

NC::ProcImpl::ProcPtr NCF::createScatter( const MatCfg& cfg )
{
  ::NC::detail::factThreads_checkEnvVar();
  if ( cfg.hasDensityOverride() )
    return createScatter( cfg.cloneWithoutDensityState() );//never matters for a process
  MatCfg::PhaseChoices phaseChoices = cfg.getPhaseChoices();
  if ( !phaseChoices.empty() ) {
    //Let createInfo deal with phasechoices, and let us just transform it into a scatter request:
    return FactImpl::create( ScatterRequest(createInfo(cfg)) );
  }
  nc_assert_always( cfg.getPhaseChoices().empty() );
  nc_assert_always( !cfg.hasDensityOverride() );
  auto cfg_phases = cfg.phases();
  nc_assert_always( cfg_phases.size() != 1 );
  return ( cfg_phases.empty()
           ? create( ScatterRequest(cfg) )
           : createProcFromMPCfg<ProcessType::Scatter,ScatterRequest>(cfg_phases) );
}

NC::ProcImpl::ProcPtr NCF::createAbsorption( const MatCfg& cfg )
{
  ::NC::detail::factThreads_checkEnvVar();
  if ( cfg.hasDensityOverride() )
    return createAbsorption( cfg.cloneWithoutDensityState() );//never matters for a process

  MatCfg::PhaseChoices phaseChoices = cfg.getPhaseChoices();
  if ( !phaseChoices.empty() ) {
    //Let createInfo deal with phasechoices, and let us just transform it into a scatter request:
    return FactImpl::create( AbsorptionRequest(createInfo(cfg) ));
  }
  auto cfg_phases = cfg.phases();
  nc_assert_always( cfg_phases.size() != 1 );
  return ( cfg_phases.empty()
           ? create( AbsorptionRequest(cfg) )
           : createProcFromMPCfg<ProcessType::Absorption,AbsorptionRequest>(cfg_phases) );
}

NC::ProcImpl::ProcPtr NCF::ScatterFactory::globalCreateScatter( const ScatterRequest& req ) const
{
  auto scatfact_strview = Cfg::CfgManip::get_scatfactory(req.rawCfgData());
  auto fnreq = Cfg::FactNameRequest::Parser::doParse( scatfact_strview );
  StrView ourname = name();
  if (!fnreq.excludes(ourname))
    fnreq = fnreq.withAdditionalExclude(ourname);
  if ( ourname == fnreq.specificRequest() )
    fnreq = fnreq.withNoSpecificRequest();
  std::string modstr;
  modstr.reserve(128);
  constexpr StrView varname = Cfg::varInfo( Cfg::VarId::scatfactory ).nameSV();
  varname.appendToString(modstr);
  modstr += '=';
  modstr.append(fnreq.to_string());
  return ::NCrystal::FactImpl::createScatter(req.modified(modstr));
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

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

#include "NCrystal/NCFactory.hh"
#include "NCrystal/NCFactoryRegistry.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCScatterComp.hh"
#include "NCrystal/NCAbsorption.hh"
#include "NCrystal/NCFile.hh"
#include "NCrystal/internal/NCString.hh"
#include <iostream>
#include <atomic>
namespace NC = NCrystal;

namespace NCrystal {

  static std::atomic<bool> s_info_cache_enabled( ! ncgetenv_bool("NOCACHE") );
  static std::atomic<bool> s_debug_factory( ncgetenv_bool("DEBUGFACTORY") );

  struct FactoryCfgSpy : public MatCfg::AccessSpy {
    std::set<std::string> parnames;
    FactoryCfgSpy() {}
    virtual ~FactoryCfgSpy(){};
    virtual void parAccessed(const std::string& parname) { parnames.insert(parname); }
  };

  struct InfoCache {
    std::set<std::string> parnames;
    std::string signature;
    RCHolder<const Info> infoholder;
    bool operator<(const InfoCache&o) const
    {
      //sort by parnames first, which reduces number of calls to
      //getCacheSignature during access.
      if (parnames!=o.parnames)
        return parnames < o.parnames;
      //then signature:
      if (signature==o.signature)
        NCRYSTAL_THROW(LogicError,"Cache inconsistency detected!");
      return signature < o.signature;
    }
  };

  static std::mutex s_infocache_mutex;//For now, should move to CachedFactoryBase implementation!
  static std::map<std::string, std::set<InfoCache> > s_infocache;

  const Info * searchInfoCache(const std::string& key, const MatCfg& cfg) {
    std::map<std::string, std::set<InfoCache> >::const_iterator itKey = s_infocache.find(key);
    if (itKey==s_infocache.end())
      return 0;
    const std::set<InfoCache>& caches = itKey->second;
    //found caches for the factory in question, now search them.
    std::string signature;
    std::set<std::string> signature_parnames;
    std::set<InfoCache>::iterator it(caches.begin()), itE(caches.end());
    for (;it!=itE;++it) {
      if (it->parnames.empty()) {
        nc_assert_always(caches.size()==1);
        //apparently the factory does not read any parameters at all, so we are always valid!
        return it->infoholder.obj();//hit!
      }
      if (signature_parnames!=it->parnames) {
        cfg.getCacheSignature(signature,it->parnames);
        signature_parnames=it->parnames;
      }
      if (signature == it->signature)
        return it->infoholder.obj();//hit!
    }
    //no hit.
    return 0;
  }

}

void NC::clearInfoCaches()
{
  //MT TODO: mutex lock (check this entire file... perhaps migrate to new infrastructure?)
  s_infocache.clear();
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory - clearInfoCaches called."<<std::endl;
}


void NC::disableCaching()
{
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory - disableCaching called."<<std::endl;
  if (!s_info_cache_enabled)
    return;
  s_info_cache_enabled = false;
  clearInfoCaches();
}

void NC::enableCaching()
{
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory - enableCaching called."<<std::endl;
  s_info_cache_enabled = true;
}

const NC::Info* NC::createInfo( const NC::MatCfg& cfg )
{
  std::lock_guard<std::mutex> guard(s_infocache_mutex);

  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createInfo - createInfo( "<<cfg<<" ) called"<<std::endl;


  static bool first = true;
  if (first) {
    registerCacheCleanupFunction( clearInfoCaches );
    first = false;
  }

  cfg.checkConsistency();
  const FactoryList& facts = getFactories();//Access factories
  std::map<int,std::shared_ptr<const FactoryBase>> avail;
  std::shared_ptr<const FactoryBase> chosen;
  std::string specific = cfg.get_infofact_name();
  if (s_debug_factory && !specific.empty())
    std::cout<<"NCrystal::Factory::createInfo - cfg.infofactory=\""<<specific<<"\" so will search for that."<<std::endl;

  for (std::size_t i = 0; i < facts.size(); ++i) {
    auto& f = facts.at(i);
    if (!specific.empty()) {
      if (specific == f->getName()) {
        chosen = std::move(f);
        if (s_debug_factory)
          std::cout<<"NCrystal::Factory::createInfo - about to invoke canCreateInfo on factory \""<<chosen->getName()<<"\""<<std::endl;
        if (!chosen->canCreateInfo(cfg))
          NCRYSTAL_THROW2(BadInput,"Requested infofactory does not actually have capability to service request: \""<<specific<<"\"");
        break;
      } else {
        continue;
      }
    }
    if (s_debug_factory)
      std::cout<<"NCrystal::Factory::createInfo - about to invoke canCreateInfo on factory \""<<f->getName()<<"\""<<std::endl;
    int priority = f->canCreateInfo(cfg);
    if (s_debug_factory) {
      std::cout<<"NCrystal::Factory::createInfo - factory \""<<f->getName()<<"\" canCreateInfo(cfg) returns ";
      if (priority) std::cout<<"YES (priority="<<priority<<")"<<std::endl;
      else std::cout<<"NO"<<std::endl;
    }
    if (priority && avail.find(priority)==avail.end())
      avail[priority] = std::move(f);
  }
  if (!specific.empty() && !chosen)
    NCRYSTAL_THROW2(BadInput,"Specific infofactory requested which is unavailable: \""<<specific<<"\"");
  if (!chosen)
    chosen = avail.empty() ? 0 : avail.rbegin()->second;
  if (!chosen)
    NCRYSTAL_THROW2(BadInput,"Could not find factory to service createInfo request ("<<facts.size()<<" factories registered)");

  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createInfo - factory \""<<chosen->getName()<<"\" chosen to service createInfo request"<<std::endl;

  std::string cachekey;
  if (s_info_cache_enabled) {
    std::stringstream cachekey_stream;
    cachekey_stream<<cfg.getDataFileAsSpecified()<<';'<<chosen->getName();
    cachekey = cachekey_stream.str();
    const Info * cached_info = searchInfoCache(cachekey, cfg);
    if (s_debug_factory)
      std::cout<<"NCrystal::Factory::createInfo - checking cache with key \""<<cachekey<<"\": "<<(cached_info?"found!":"notfound")<<std::endl;
    if (cached_info)
      return cached_info;
  }

  FactoryCfgSpy spy;
  cfg.addAccessSpy(&spy);
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createInfo - invoking createInfo on factory \""<<chosen->getName()<<"\""<<std::endl;
  auto info = chosen->createInfo(cfg);
  cfg.removeAccessSpy(&spy);
  if (!info.obj())
    NCRYSTAL_THROW(BadInput,"Chosen factory could not service createInfo request");


  //to ensure good caching + separation, we enforce dynamically that factories
  //only access a limited subset of the MatCfg parameters during calls to
  //createInfo:
  static std::set<std::string> allowed_info_pars = { "temp", "dcutoff", "dcutoffup", "atomdb", "overridefileext", "infofactory" };
  std::set<std::string>::const_iterator it = spy.parnames.begin();
  for (;it!=spy.parnames.end();++it) {
    if (!allowed_info_pars.count(*it))
      NCRYSTAL_THROW2(LogicError,"Factory \""<<chosen->getName()
                      <<"\" accessed MatCfg parameter \""<<*it<<"\" during createInfo(..) which"
                      " violates caching policies.");
  }

  if ( ! info->isLocked() )
    NCRYSTAL_THROW2(LogicError,"Factory \""<<chosen->getName()<<"\" did not lock created Info object");

  if ( cfg.get_temp()!=-1.0 ) {
    if ( !info->hasTemperature() || info->getTemperature() != cfg.get_temp() )
      NCRYSTAL_THROW2(LogicError,"Factory \""<<chosen->getName()<<"\" did not set temp as required");
  }

  if (info->hasHKLInfo()) {
    if (cfg.get_dcutoff()==-1)
      NCRYSTAL_THROW2(LogicError,"Factory \""<<chosen->getName()
                      <<"\" created HKL info even though dcutoff=-1");
    if ( info->hklDLower() < cfg.get_dcutoff() ||
         info->hklDUpper() > cfg.get_dcutoffup() )
      NCRYSTAL_THROW2(LogicError,"Factory \""<<chosen->getName()
                      <<"\" did not respect dcutoff setting.");
  }

  if (s_info_cache_enabled) {
    //Update cache:
    nc_assert(!cachekey.empty());
    std::string cache_signature;
    cfg.getCacheSignature(cache_signature,spy.parnames);
    if (s_debug_factory)
      std::cout<<"NCrystal::Factory::createInfo - update cache with key \""<<cachekey<<"\" and signature \""<<cache_signature<<"\""<<std::endl;
    InfoCache cachevalue;
    cachevalue.parnames = spy.parnames;
    cachevalue.signature = cache_signature;
    cachevalue.infoholder = info;
    std::map<std::string, std::set<InfoCache> >::iterator itCache = s_infocache.find(cachekey);
    if (itCache==s_infocache.end()) {
      std::set<InfoCache> tmp;
      tmp.insert(cachevalue);
      s_infocache[cachekey]=tmp;
    } else {
      itCache->second.insert(cachevalue);
    }
  }

  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createInfo - createInfo was successful"<<std::endl;
  return info.releaseNoDelete();
}

const NC::Scatter* NC::createScatter( const NC::MatCfg& cfg )
{
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createScatter - createScatter( "<<cfg<<" ) called"<<std::endl;

  cfg.checkConsistency();
  const FactoryList& facts = getFactories();//Access factories
  std::map<int,std::shared_ptr<const FactoryBase>> avail;
  std::shared_ptr<const FactoryBase> chosen;
  std::string specific = cfg.get_scatfactory();
  if (s_debug_factory && !specific.empty())
    std::cout<<"NCrystal::Factory::createScatter - cfg.scatfactory=\""<<specific<<"\" so will search for that."<<std::endl;

  for (std::size_t i = 0; i < facts.size(); ++i) {
    auto& f = facts.at(i);
    if (!specific.empty()) {
      if (specific == f->getName()) {
        chosen = std::move(f);
        if (s_debug_factory)
          std::cout<<"NCrystal::Factory::createScatter - about to invoke canCreateScatter on factory \""<<chosen->getName()<<"\""<<std::endl;
        if (!chosen->canCreateScatter(cfg))
          NCRYSTAL_THROW2(BadInput,"Requested scatfactory does not actually have capability to service request: \""<<specific<<"\"");
        break;
      } else {
        continue;
      }
    }
    if (s_debug_factory)
      std::cout<<"NCrystal::Factory::createScatter - about to invoke canCreateScatter on factory \""<<f->getName()<<"\""<<std::endl;
    int priority = f->canCreateScatter(cfg);
    if (s_debug_factory) {
      std::cout<<"NCrystal::Factory::createScatter - factory \""<<f->getName()<<"\" canCreateScatter(cfg) returns ";
      if (priority) std::cout<<"YES (priority="<<priority<<")"<<std::endl;
      else std::cout<<"NO"<<std::endl;
    }
    if (priority && avail.find(priority)==avail.end())
      avail[priority] = std::move(f);
  }
  if (!specific.empty() && !chosen)
    NCRYSTAL_THROW2(BadInput,"Specific scatfactory requested which is unavailable: \""<<specific<<"\"");
  if (!chosen)
    chosen = avail.empty() ? 0 : avail.rbegin()->second;
  if (!chosen)
    NCRYSTAL_THROW2(BadInput,"Could not find factory to service createScatter request ("<<facts.size()<<" factories registered)");

  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createScatter - factory \""<<chosen->getName()<<"\" chosen to service createScatter request"<<std::endl;

  auto scatter = chosen->createScatter(cfg);
  if (!scatter)
    NCRYSTAL_THROW(BadInput,"Chosen factory could not service createScatter request");
  if (s_debug_factory) {
    const char * prefix = "NCrystal::Factory::createScatter::success ";
    const ScatterComp * scatcomp = dynamic_cast<const ScatterComp*>(scatter.obj());
    std::cout<<prefix<<std::endl;
    std::cout<<prefix<<" createScatter was successful and resulted in "<<scatter->getCalcName()<<" object";
    if ( scatcomp ){
      std::cout<<" with "<<scatcomp->nComponents()<<" components:"<<std::endl;
      for (std::size_t i = 0; i < scatcomp->nComponents(); ++i )
        std::cout<<prefix<<"     => "<<scatcomp->component(i)->getCalcName()<<" (with scale "<<scatcomp->scale(i)<<")"<<std::endl;
    } else {
      std::cout<<std::endl;
    }
    std::cout<<prefix<<std::endl;
  }

  return scatter.releaseNoDelete();
}


//TODO: this function is a cut'n'paste of the above, with scatter->absorption replacement.
const NC::Absorption* NC::createAbsorption( const NC::MatCfg& cfg )
{
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createAbsorption - createAbsorption( "<<cfg<<" ) called"<<std::endl;

  cfg.checkConsistency();
  const FactoryList& facts = getFactories();//Access factories
  std::map<int,std::shared_ptr<const FactoryBase>> avail;
  std::shared_ptr<const FactoryBase> chosen;
  std::string specific = cfg.get_absnfactory();
  if (s_debug_factory && !specific.empty())
    std::cout<<"NCrystal::Factory::createAbsorption - cfg.absnfactory=\""<<specific<<"\" so will search for that."<<std::endl;

  for (std::size_t i = 0; i < facts.size(); ++i) {
    auto& f = facts.at(i);
    if (!specific.empty()) {
      if (specific == f->getName()) {
        chosen = std::move(f);
        if (s_debug_factory)
          std::cout<<"NCrystal::Factory::createAbsorption - about to invoke canCreateAbsorption on factory \""<<chosen->getName()<<"\""<<std::endl;
        if (!chosen->canCreateAbsorption(cfg))
          NCRYSTAL_THROW2(BadInput,"Requested absnfactory does not actually have capability to service request: \""<<specific<<"\"");
        break;
      } else {
        continue;
      }
    }
    if (s_debug_factory)
      std::cout<<"NCrystal::Factory::createAbsorption - about to invoke canCreateAbsorption on factory \""<<f->getName()<<"\""<<std::endl;
    int priority = f->canCreateAbsorption(cfg);
    if (s_debug_factory) {
      std::cout<<"NCrystal::Factory::createAbsorption - factory \""<<f->getName()<<"\" canCreateAbsorption(cfg) returns ";
      if (priority) std::cout<<"YES (priority="<<priority<<")"<<std::endl;
      else std::cout<<"NO"<<std::endl;
    }
    if (priority && avail.find(priority)==avail.end())
      avail[priority] = std::move(f);
  }
  if (!specific.empty() && !chosen)
    NCRYSTAL_THROW2(BadInput,"Specific absnfactory requested which is unavailable: \""<<specific<<"\"");
  if (!chosen)
    chosen = avail.empty() ? 0 : avail.rbegin()->second;
  if (!chosen)
    NCRYSTAL_THROW2(BadInput,"Could not find factory to service createAbsorption request ("<<facts.size()<<" factories registered)");

  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createAbsorption - factory \""<<chosen->getName()<<"\" chosen to service createAbsorption request"<<std::endl;

  auto absorption = chosen->createAbsorption(cfg);
  if (!absorption)
    NCRYSTAL_THROW(BadInput,"Chosen factory could not service createAbsorption request");
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createAbsorption - createAbsorption was successful"<<std::endl;
  return absorption.releaseNoDelete();
}

namespace NCrystal {

#ifdef NCRYSTAL_STDCMAKECFG_EMBED_DATA_ON
  //If NCrystal is installed using the standard CMake with -DEMBED_DATA=ON, an
  //autogenerated function NCrystal::AutoGenNCMAT::registerStdNCMAT will have
  //been compiled in, and we have to call it at initialisation time (after all
  //static objects everywhere have been initialised!)
  namespace AutoGenNCMAT { void registerStdNCMAT(); }//fwd declared - linked in elsewhere
#endif

  namespace {
    class InMemoryFileDB;
    //Keep inmemdb mutex in shared pointer, so it wont be deleted before
    //InMemoryFileDB objects gets deleted (destruction order between static
    //objects in different compilation units is undefined):
    static std::shared_ptr<std::mutex> s_inmemdb_mutex = std::make_shared<std::mutex>();
    static InMemoryFileDB * s_inmemdb = nullptr;

    class InMemoryFileDB : public TextInputManager {
      struct Entry {
        const char * staticData = nullptr;
        std::string data;
      };
      std::map<std::string,Entry> m_db;
      std::shared_ptr<std::mutex> m_mutex;

    public:
      InMemoryFileDB( std::shared_ptr<std::mutex> mm )
        : m_mutex(std::move(mm))
      {
      }
      virtual ~InMemoryFileDB()
      {
        //If deleted in NCFile.cc (which owns a TextInputManager once
        //registered), make sure the non-owning static pointer here is updated
        //as well.
        assert(m_mutex!=nullptr);//NB: Only standard asserts in destructors (can't throw).
        std::lock_guard<std::mutex> guard(*m_mutex);//NB: The mutex lock here is why we need to keep it in a shared pointer!
        s_inmemdb = nullptr;
        //TODO: How can we systematically check other files for this static-destructors-ordering-in-different-units issues??
      }

      void clearCaches(const std::string& name)
      {
        //Clear any existing info caches related to this name
        std::lock_guard<std::mutex> guard(s_infocache_mutex);
        std::string searchpattern(name+";");
        auto itE = s_infocache.end();
        for (auto it = s_infocache.begin(); it!=itE;) {
          if (startswith(it->first,searchpattern)) {
            auto itdel = it;
            ++it;
            s_infocache.erase(itdel);
          } else {
            ++it;
          }
        }
      }

      void addEntry(const std::string& name,
                    std::string&& data)
      {
        Entry entry;
        entry.data = std::move(data);
        m_db[name] = std::move(entry);
        clearCaches(name);
      }

      void addStaticEntry(const std::string& name,
                          const char * staticData)
      {
        Entry entry;
        entry.staticData = staticData;
        m_db[name] = std::move(entry);
        clearCaches(name);
      }

      //Reimplement this custom file searching in this function (can throw
      //FileNotFound in case of problems, but one will in any case be thrown if it
      //doesn't supply a result and if the fallback to the usual search patterns
      //is disallowed or fails):
      std::unique_ptr<TextInputStream> createTextInputStream( const std::string& name ) final
      {
        nc_assert(!!m_mutex);
        std::lock_guard<std::mutex> guard(*m_mutex);
        auto it = m_db.find(name);
        if ( it == m_db.end() )
          return nullptr;//Do not throw FileNotFound here (will prevent on-disk file usage).
        return createTextInputStreamFromBuffer( name,
                                                ( it->second.staticData
                                                  ? it->second.staticData
                                                  : it->second.data ) );
      }

      std::vector<PairSS> getList() const final {
        nc_assert( m_mutex != nullptr );
        std::lock_guard<std::mutex> guard(*m_mutex);//NB: The mutex lock here is why we need to keep it in a shared pointer!
        std::vector<PairSS> list;
        list.reserve(m_db.size());
        for (const auto& e : m_db)
          list.emplace_back(e.first,"<embedded>");
        return list;
      }

    };
    void ensureDBReady() {
      //Assumes mutex is already locked by calling code.
      if (!s_inmemdb) {
        auto tmp = std::make_unique<InMemoryFileDB>(s_inmemdb_mutex);
        s_inmemdb = tmp.get();
        registerTextInputManager(std::move(tmp));
#ifdef NCRYSTAL_STDCMAKECFG_EMBED_DATA_ON
        AutoGenNCMAT::registerStdNCMAT();
#endif
      }
    }
  }

#ifdef NCRYSTAL_STDCMAKECFG_EMBED_DATA_ON
  namespace internal {
    //Other functions needed for the embedding:
    void registerEmbeddedNCMAT( const char* name, const char* static_data )
    {
      //Unsafe version which is intended to be called only when s_inmemdb_mutex
      //is already locked and s_inmemdb is already setup.
      nc_assert_always(s_inmemdb);
      s_inmemdb->addStaticEntry(name,static_data);
    }
    void ensureInMemDBReadyMTSafe() {
      nc_assert(!!s_inmemdb_mutex);
      std::lock_guard<std::mutex> guard(*s_inmemdb_mutex);
      ensureDBReady();
    }
  }
#endif
}

void NC::registerInMemoryFileData( const std::string& name,
                                   const std::string& data )
{
  nc_assert(!!s_inmemdb_mutex);
  std::lock_guard<std::mutex> guard(*s_inmemdb_mutex);
  ensureDBReady();
  s_inmemdb->addEntry(name,std::string(data));
}


void NC::registerInMemoryFileData( const std::string& name,
                                   std::string&& data )
{
  nc_assert(!!s_inmemdb_mutex);
  std::lock_guard<std::mutex> guard(*s_inmemdb_mutex);
  ensureDBReady();
  s_inmemdb->addEntry(name,std::move(data));
}

void NC::registerInMemoryStaticFileData( const std::string& name,
                                         const char* static_data )
{
  nc_assert(!!s_inmemdb_mutex);
  std::lock_guard<std::mutex> guard(*s_inmemdb_mutex);
  ensureDBReady();
  s_inmemdb->addStaticEntry(name,static_data);
}

void NC::ensureEmbeddedDataIsRegistered()
{
#ifdef NCRYSTAL_STDCMAKECFG_EMBED_DATA_ON
  internal::ensureInMemDBReadyMTSafe();
#endif
}

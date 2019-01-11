////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2019 NCrystal developers                                   //
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
#include <map>
#include <set>
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCScatter.hh"
#include "NCrystal/NCAbsorption.hh"
#include <iostream>
#include <cstdlib>

namespace NCrystal {

  static bool s_info_cache_enabled = (std::getenv("NCRYSTAL_NOCACHE") ? false : true);

  static bool s_debug_factory = (std::getenv("NCRYSTAL_DEBUGFACTORY") ? true : false);

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

  static std::map<std::string, std::set<InfoCache> > s_infocache;

  const Info * searchInfoCache(const std::string& key, const NCrystal::MatCfg& cfg) {
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

void NCrystal::clearInfoCaches()
{
  s_infocache.clear();
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory - clearInfoCaches called."<<std::endl;
}


void NCrystal::disableCaching()
{
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory - disableCaching called."<<std::endl;
  if (!s_info_cache_enabled)
    return;
  s_info_cache_enabled = false;
  clearInfoCaches();
}

void NCrystal::enableCaching()
{
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory - enableCaching called."<<std::endl;
  s_info_cache_enabled = true;
}

const NCrystal::Info * NCrystal::createInfo( const NCrystal::MatCfg& cfg )
{
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createInfo - createInfo( "<<cfg<<" ) called"<<std::endl;

  cfg.checkConsistency();
  FactoryList& facts = getFactories();//Access factories
  std::map<int,const FactoryBase*> avail;
  const FactoryBase* chosen = 0;
  std::string specific = cfg.get_infofact_name();
  if (s_debug_factory && !specific.empty())
    std::cout<<"NCrystal::Factory::createInfo - cfg.infofactory=\""<<specific<<"\" so will search for that."<<std::endl;

  for (std::size_t i = 0; i < facts.size(); ++i) {
    const FactoryBase* f = facts.at(i);
    if (!specific.empty()) {
      if (specific == f->getName()) {
        chosen = f;
        if (s_debug_factory)
          std::cout<<"NCrystal::Factory::createInfo - about to invoke canCreateInfo on factory \""<<f->getName()<<"\""<<std::endl;
        if (!f->canCreateInfo(cfg))
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
      avail[priority] = f;
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
    cachekey_stream<<cfg.getDataFile()<<';'<<chosen->getName();
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
  RCHolder<const Info> info(chosen->createInfo(cfg));
  cfg.removeAccessSpy(&spy);
  if (!info.obj())
    NCRYSTAL_THROW(BadInput,"Chosen factory could not service createInfo request");
  if (info.obj()->refCount()!=1)//1 here since RCHolder already incremented
    NCRYSTAL_THROW(BadInput,"Chosen factory returned object with non-zero reference count!");


  //to ensure good caching + separation, we enforce dynamically that factories
  //only access a limited subset of the MatCfg parameters during calls to
  //createInfo:
  static std::set<std::string> allowed_info_pars
#if __cplusplus >= 201103L
    = { "temp", "dcutoff", "dcutoffup", "overridefileext", "infofactory" };
#else
  ;
  if (allowed_info_pars.empty()) {
    allowed_info_pars.insert("temp");
    allowed_info_pars.insert("dcutoff");
    allowed_info_pars.insert("dcutoffup");
    allowed_info_pars.insert("overridefileext");
    allowed_info_pars.insert("infofactory");
  }
#endif
  std::set<std::string>::const_iterator it = spy.parnames.begin();
  for (;it!=spy.parnames.end();++it) {
    if (!allowed_info_pars.count(*it))
      NCRYSTAL_THROW2(LogicError,"Factory \""<<chosen->getName()
                      <<"\" accessed MatCfg parameter \""<<*it<<"\" during createInfo(..) which"
                      " violates caching policies.");
  }

  if ( ! info.obj()->isLocked() )
    NCRYSTAL_THROW2(LogicError,"Factory \""<<chosen->getName()<<"\" did not lock created Info object");

  if ( info.obj()->hasTemperature() && info.obj()->getTemperature() != cfg.get_temp() )
    NCRYSTAL_THROW2(LogicError,"Factory \""<<chosen->getName()<<"\" did not set temp as required");

  if (info.obj()->hasHKLInfo()) {
    if (cfg.get_dcutoff()==-1)
      NCRYSTAL_THROW2(LogicError,"Factory \""<<chosen->getName()
                      <<"\" created HKL info even though dcutoff=-1");
    if ( info.obj()->hklDLower() < cfg.get_dcutoff() ||
         info.obj()->hklDUpper() > cfg.get_dcutoffup() )
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

  //careful when getting the object out that its refcount doesn't drop to zero
  //and trigger cleanup, since the caching above might be disabled.
  const Info * o = info.obj();
  o->ref();
  info.clear();
  o->unrefNoDelete();
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createInfo - createInfo was succesful"<<std::endl;
  return o;
}

const NCrystal::Scatter * NCrystal::createScatter( const NCrystal::MatCfg& cfg )
{
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createScatter - createScatter( "<<cfg<<" ) called"<<std::endl;

  cfg.checkConsistency();
  FactoryList& facts = getFactories();//Access factories
  std::map<int,const FactoryBase*> avail;
  const FactoryBase* chosen = 0;
  std::string specific = cfg.get_scatfactory();
  if (s_debug_factory && !specific.empty())
    std::cout<<"NCrystal::Factory::createScatter - cfg.scatfactory=\""<<specific<<"\" so will search for that."<<std::endl;

  for (std::size_t i = 0; i < facts.size(); ++i) {
    const FactoryBase* f = facts.at(i);
    if (!specific.empty()) {
      if (specific == f->getName()) {
        chosen = f;
        if (s_debug_factory)
          std::cout<<"NCrystal::Factory::createScatter - about to invoke canCreateScatter on factory \""<<f->getName()<<"\""<<std::endl;
        if (!f->canCreateScatter(cfg))
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
      avail[priority] = f;
  }
  if (!specific.empty() && !chosen)
    NCRYSTAL_THROW2(BadInput,"Specific scatfactory requested which is unavailable: \""<<specific<<"\"");
  if (!chosen)
    chosen = avail.empty() ? 0 : avail.rbegin()->second;
  if (!chosen)
    NCRYSTAL_THROW2(BadInput,"Could not find factory to service createScatter request ("<<facts.size()<<" factories registered)");

  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createScatter - factory \""<<chosen->getName()<<"\" chosen to service createScatter request"<<std::endl;

  const Scatter * scatter = chosen->createScatter(cfg);
  if (!scatter)
    NCRYSTAL_THROW(BadInput,"Chosen factory could not service createScatter request");
  if (scatter->refCount()!=0)
    NCRYSTAL_THROW(BadInput,"Chosen factory returned object with non-zero reference count!");
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createScatter - createScatter was succesful"<<std::endl;
  return scatter;
}


//TODO for NC2: this function is a cut'n'paste of the above, with scatter->absorption replacement.
const NCrystal::Absorption * NCrystal::createAbsorption( const NCrystal::MatCfg& cfg )
{
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createAbsorption - createAbsorption( "<<cfg<<" ) called"<<std::endl;

  cfg.checkConsistency();
  FactoryList& facts = getFactories();//Access factories
  std::map<int,const FactoryBase*> avail;
  const FactoryBase* chosen = 0;
  std::string specific = cfg.get_absnfactory();
  if (s_debug_factory && !specific.empty())
    std::cout<<"NCrystal::Factory::createAbsorption - cfg.absnfactory=\""<<specific<<"\" so will search for that."<<std::endl;

  for (std::size_t i = 0; i < facts.size(); ++i) {
    const FactoryBase* f = facts.at(i);
    if (!specific.empty()) {
      if (specific == f->getName()) {
        chosen = f;
        if (s_debug_factory)
          std::cout<<"NCrystal::Factory::createAbsorption - about to invoke canCreateAbsorption on factory \""<<f->getName()<<"\""<<std::endl;
        if (!f->canCreateAbsorption(cfg))
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
      avail[priority] = f;
  }
  if (!specific.empty() && !chosen)
    NCRYSTAL_THROW2(BadInput,"Specific absnfactory requested which is unavailable: \""<<specific<<"\"");
  if (!chosen)
    chosen = avail.empty() ? 0 : avail.rbegin()->second;
  if (!chosen)
    NCRYSTAL_THROW2(BadInput,"Could not find factory to service createAbsorption request ("<<facts.size()<<" factories registered)");

  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createAbsorption - factory \""<<chosen->getName()<<"\" chosen to service createAbsorption request"<<std::endl;

  const Absorption * absorption = chosen->createAbsorption(cfg);
  if (!absorption)
    NCRYSTAL_THROW(BadInput,"Chosen factory could not service createAbsorption request");
  if (absorption->refCount()!=0)
    NCRYSTAL_THROW(BadInput,"Chosen factory returned object with non-zero reference count!");
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory::createAbsorption - createAbsorption was succesful"<<std::endl;
  return absorption;
}

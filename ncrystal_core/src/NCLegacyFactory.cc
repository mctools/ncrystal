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

#include "NCrystal/NCLegacyFactory.hh"
#include "NCrystal/NCFactImpl.hh"
#include "NCrystal/internal/NCString.hh"
#include <iostream>
#include <atomic>
namespace NC = NCrystal;
namespace NCL = NCrystal::Legacy;


namespace NCrystal {
  namespace Legacy {
    namespace {

      static std::atomic<bool> s_info_cache_enabled( ! ncgetenv_bool("NOCACHE") );
      static std::atomic<bool> s_debug_factory( ncgetenv_bool("DEBUGFACTORY") || ncgetenv_bool("DEBUG_FACT") );

      struct LegacyInfoCache {
        std::map<UniqueIDValue,RCHolder<const Legacy::Info>> cache;
        std::mutex mtx;
      };
      LegacyInfoCache& legacyInfoCache() { static LegacyInfoCache s_cache; return s_cache; }
    }
  }
}

const NCL::Scatter* NCL::createScatter( const NC::MatCfg& cfg )
{
  //The legacy scatter factories always returned new objects:
  return wrapModernProcPtrInLegacyScatterClass(FactImpl::createScatter(cfg)).releaseNoDelete();
}

const NCL::Absorption* NCL::createAbsorption( const NC::MatCfg& cfg )
{
  //The legacy absorption factories always returned new objects:
  return wrapModernProcPtrInLegacyAbsorptionClass(FactImpl::createAbsorption(cfg)).releaseNoDelete();
}

const NCL::Info* NCL::createInfo( const NC::MatCfg& cfg )
{
  auto matinfo = FactImpl::createInfo(cfg);
  if ( !s_info_cache_enabled )
    return makeRC<const Info>(std::move(matinfo)).releaseNoDelete();
  auto uid = matinfo->getUniqueID();
  auto& c = legacyInfoCache();
  NCRYSTAL_LOCK_GUARD(c.mtx);
  c.cache.clear();
  static bool first = true;
  if (first) {
    first = false;
    registerCacheCleanupFunction( clearInfoCaches );
  }
  auto it = c.cache.find(uid);
  if ( it != c.cache.end() )
    return it->second.obj();
  auto iii = makeRC<const Info>( std::move(matinfo) );
  c.cache[uid] = iii;
  return iii.obj();
}

void NCL::clearInfoCaches()
{
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory - clearInfoCaches called."<<std::endl;
  {
    //Clear legacy cache:
    auto& c = legacyInfoCache();
    NCRYSTAL_LOCK_GUARD(c.mtx);
    c.cache.clear();
  }
  {
    //Clear new info caches (and, well everything):
    NC::clearCaches();
  }
}

void NCL::disableCaching()
{
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory - disableCaching called."<<std::endl;
  s_info_cache_enabled = false;
  FactImpl::setCachingEnabled(false);
}

void NCL::enableCaching()
{
  if (s_debug_factory)
    std::cout<<"NCrystal::Factory - enableCaching called."<<std::endl;
  s_info_cache_enabled = true;
  FactImpl::setCachingEnabled(true);
}

#include "NCrystal/NCDataSources.hh"
void NCL::registerInMemoryFileData( std::string virtualFileName,
                                    std::string&& data )
{
  DataSources::registerInMemoryFileData( std::move(virtualFileName),
                                         std::move(data) );
}

void NCL::registerInMemoryStaticFileData( std::string virtualFileName,
                                          const char* static_data )
{
  DataSources::registerInMemoryStaticFileData( std::move(virtualFileName),
                                               static_data );
}

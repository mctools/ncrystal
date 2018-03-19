////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

#include "NCrystal/NCBkgdPhonDebye.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCException.hh"
#include "NCBkgdPhonDebyeXS.hh"
#include <map>
#include <iostream>
#include <cassert>
#include <cstdlib>
#if __cplusplus >= 201103L
#  include <mutex>
#endif

namespace NCrystal {
  bool BkgdPhonDebye_checkinfo(const Info* ci, bool throwonerr )
  {
    nc_assert_always(ci);
    if ( ! ci->hasAtomInfo() ) {
      if (throwonerr) {
        NCRYSTAL_THROW(MissingInfo,"Passed info object lacks AtomInfo");
      } else {
        return false;
      }
    }
    if ( ! ci->hasTemperature() ) {
      if (throwonerr) {
        NCRYSTAL_THROW(MissingInfo,"Passed info object lacks Temperature information");
      } else {
        return false;
      }
    }
    if ( ! ci->hasXSectFree() ) {
      if (throwonerr) {
        NCRYSTAL_THROW(MissingInfo,"Passed info object lacks XSectFree information");
      } else {
        return false;
      }
    }
    if ( ! ( ci->hasPerElementDebyeTemperature() || ci->hasDebyeTemperature() ) ) {
      if (throwonerr) {
        NCRYSTAL_THROW(MissingInfo,"Passed info object lacks Debye temperature information");
      } else {
        return false;
      }
    }
    return true;
  }
  namespace BkgdPhonDebyeCache {
    typedef std::pair<uint64_t,uint64_t> Key;
    typedef std::map<Key,RCHolder<const BkgdPhonDebyeXS> > Map;
    static Map cache;
#if __cplusplus >= 201103L
    static std::mutex cache_mutex;
    //cache_mutex is used to protect cache access in case of multi-threading
    //(access should be rare, so a simple mutex approach is probably just fine).
#endif
  }
}

bool NCrystal::BkgdPhonDebye::hasSufficientInfo(const Info* ci)
{
  return ci && BkgdPhonDebye_checkinfo(ci,false);
}


NCrystal::BkgdPhonDebye::BkgdPhonDebye(const Info* ci, bool thermalise,
                                       unsigned nphonon, bool pzi, bool only_pzi, bool no_extrap )
  : ScatterXSCurve(ci,"BkgdPhonDebye",thermalise), m_xs(0)
{
  nc_assert_always(ci);
  BkgdPhonDebye_checkinfo(ci,true);

  RCHolder<const BkgdPhonDebyeXS> xs;
  const bool verbose = (std::getenv("NCRYSTAL_DEBUGSCATTER") ? true : false);
  //Construct unique cachekey for NCInfo instance and parameters above:
  BkgdPhonDebyeCache::Key cachekey(ci->getUniqueID(),0);
  if (pzi) cachekey.second += 1;
  if (only_pzi) cachekey.second += 2;
  if (no_extrap) cachekey.second += 4;
  cachekey.second += 8*nphonon;
  {
#if __cplusplus >= 201103L
    std::lock_guard<std::mutex> lock(BkgdPhonDebyeCache::cache_mutex);
#endif
    //Search cache:
    BkgdPhonDebyeCache::Map::iterator it = BkgdPhonDebyeCache::cache.find(cachekey);
    if (it==BkgdPhonDebyeCache::cache.end()) {
      //not in cache, must create and insert in cache:
      if (verbose)
        std::cout<<"NCBkgdPhonDebye failed to find relevant BkgdPhonDebyeXS object in cache (cachelength="
                 <<BkgdPhonDebyeCache::cache.size()<<"). Creating from scratch."<<std::endl;
      xs = createBkgdPhonDebyeXS(ci,nphonon,pzi,only_pzi,!no_extrap);
      BkgdPhonDebyeCache::cache[cachekey] = xs;
    } else {
      //found in cache!
      nc_assert_always(it->second.obj());
      if (verbose)
        std::cout<<"NCBkgdPhonDebye found relevant BkgdPhonDebyeXS object in cache (cachelength="
                 <<BkgdPhonDebyeCache::cache.size()<<")."<<std::endl;
      xs = it->second;
    }
    nc_assert_always(xs.obj()&&xs.obj()->refCount()>1);
    m_xs = xs.obj();
    m_xs->ref();
  }
  validate();
}

NCrystal::BkgdPhonDebye::~BkgdPhonDebye()
{
  if (m_xs && !std::getenv("NCRYSTAL_NEVERCLEARCACHES") ) {
#if __cplusplus >= 201103L
    std::lock_guard<std::mutex> lock(BkgdPhonDebyeCache::cache_mutex);
#endif
    //Remove from cache if we are the last user.
    BkgdPhonDebyeCache::Map::iterator it(BkgdPhonDebyeCache::cache.begin()),
                                      itE(BkgdPhonDebyeCache::cache.end());
    for (;it!=itE;++it) {
      if (it->second.obj()==m_xs) {
        assert(it->second.obj()->refCount()>=2);//normal assert, since we can't throw from destructors
        if (it->second.obj()->refCount()==2) {
          //cache itself holds 1 ref, we hold one in m_xs, and noone else does.
          BkgdPhonDebyeCache::cache.erase(it);
          if (std::getenv("NCRYSTAL_DEBUGSCATTER"))
            std::cout<<"NCBkgdPhonDebye destructor removed BkgdPhonDebyeXS object from cache (cachelength="
                     <<BkgdPhonDebyeCache::cache.size()<<")."<<std::endl;
        }
        break;
      }
    }
    m_xs->unref();
  }
}

double NCrystal::BkgdPhonDebye::crossSectionNonOriented(double ekin) const
{
  return m_xs->getXS(ekin);
}

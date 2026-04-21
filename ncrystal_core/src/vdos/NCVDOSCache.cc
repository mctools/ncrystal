////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

#include "NCrystal/internal/vdos/NCVDOSCache.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/fact_utils/NCFactoryUtils.hh"

namespace NC = NCrystal;

NC::DI_VDOSShPtr::~DI_VDOSShPtr() = default;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    static_assert( std::is_same<HashValue,std::size_t>::value, "" );
    HashValue vdosDataHash( const VDOSData& vd )
    {
      HashValue hash = calcHash( vd.elementMassAMU().dbl() );
      //Given the already ~0% clash probability and that our hash does not need
      //to be clash free, we don't hash the entire vd.vdos_density() vector.
      hash_combine(hash, vd.vdos_density().size() );
      hash_combine(hash, vd.vdos_density().front() );
      hash_combine(hash, vd.vdos_density().back() );
      hash_combine(hash, vd.vdos_egrid().first );
      hash_combine(hash, vd.vdos_egrid().second );
      hash_combine(hash, vd.temperature().dbl() );
      hash_combine(hash, vd.boundXS().dbl() );
      return hash;
    }

    struct VDOSCacheEntry {
      //Cache entry. In addition to the payload, we also keep track of when it
      //was last accessed. In case cache needs to be trimmed, this allows us to
      //remove those accessed furthest back in time first.
      VDOSCacheEntry(VDOSDataHashPtr pl,  std::size_t la )
        : payload( std::move(pl) ), last_accessed(la) {}
      VDOSDataHashPtr payload;
      std::size_t last_accessed = 0;
      //Default comparisons are on payload only
      bool operator<(const VDOSCacheEntry&o) const { return payload < o.payload; }
      bool operator<(const VDOSDataHashPtr&o) const { return payload < o; }
    };
    struct VDOSCacheDB {
      std::mutex mtx;
      std::size_t call_count = 0;
      std::vector<VDOSCacheEntry> cache;//keep sorted
    };
    VDOSCacheDB& getVDOSCacheDB()
    {
      static VDOSCacheDB db;
      static const bool dummy = []()
      {
        registerCacheCleanupFunction([](){
          auto& db2 = getVDOSCacheDB();
          NCRYSTAL_LOCK_GUARD(db2.mtx);
          db2.cache.clear();
        });
        return false;
       }();
      (void) dummy;
      return db;
    }

  }
}

NC::VDOSDataHashPtr::VDOSDataHashPtr( shared_obj<const VDOSData> data )
  : m_data( std::move(data) ),
    m_hash( vdosDataHash(m_data) )
{
}

NC::VDOSDataHashPtr NC::getCachedVDOSDataHashPtr( VDOSData&& data )
{
  VDOSDataHashPtr newobj( makeSO<VDOSData>(std::move(data)) );

  //Check if object identical to newobj is already in cache:
  auto& db = getVDOSCacheDB();
  NCRYSTAL_LOCK_GUARD(db.mtx);

  ++db.call_count;

  auto it = std::lower_bound(db.cache.begin(), db.cache.end(), newobj);
  if (it != db.cache.end() && newobj == it->payload ) {
    //Already found! Note down access and return:
    it->last_accessed = db.call_count;
    newobj = it->payload;
    return newobj;
  }
  //Not found, so we must insert and sort cache. Potentially we might also trim
  //the cache if it grows too wildly. For simplicity: if it would exceed 1000 we
  //discard the 100 entries that were accessed longest ago. We could potentially
  //improve this by checking via a weakptr if any of those are still alive
  //somewhere before discarding them, but it seems a bit of a needless
  //complication.
  if ( db.cache.size()+1 > 1000 ) {
    //Trim it down to, say, 900. We do it by sorting on .last_accessed,
    //discarding entries. The final sort below will take care of putting
    //everything back in order.
    std::sort( db.cache.begin(), db.cache.end(),
               [](const VDOSCacheEntry& a, const VDOSCacheEntry& b)
               { return a.last_accessed > b.last_accessed; });
    while ( db.cache.size() > 900 )
      db.cache.pop_back();
    //FIXME: Test code with very small limit, and also test cache clearing!
  }
  db.cache.emplace_back( newobj, db.call_count );
  std::sort( db.cache.begin(), db.cache.end() );
  return newobj;
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    struct EGCacheDB {
      std::mutex mtx;
      std::vector<EnergyGridHashPtr> cache;
    };
    EGCacheDB& getEGCacheDB()
    {
      static EGCacheDB db;
      static const bool dummy = []()
      {
        registerCacheCleanupFunction([](){
          auto& db2 = getEGCacheDB();
          NCRYSTAL_LOCK_GUARD(db2.mtx);
          db2.cache.clear();
        });
        return false;
      }();
      (void) dummy;
      return db;
    }

  }
  class EnergyGridHashPtrFact {
  public:
    static EnergyGridHashPtr create( std::shared_ptr<const VectD> v,
                                     std::size_t hash )
    {
      return { std::move(v), hash, UniqueID().getUniqueID() };
    }
  };
}

NC::EnergyGridHashPtr NC::getCachedEnergyGridHashPtr( VectD&& v )
{
  if ( v.empty() ) {
    //special case the empty (or null) grid (hash = 0). Always return same
    //instance with same UID:
    static auto egridnull = EnergyGridHashPtrFact::create( nullptr, 0 );
    return egridnull.clone();
  }
  //calc hash:
  nc_assert_always(v.size()>=3);
  HashValue hash = calcHash( v[0] );
  hash_combine( hash, v[1] );
  hash_combine( hash, v.back() );
  hash_combine( hash, v.size() );
  if (!hash)
    hash = 1;//make sure hash=0 <=> empty
  //Access DB under lock:
  auto& db = getEGCacheDB();
  NCRYSTAL_LOCK_GUARD(db.mtx);

  //FIXME: VERY Simplistic and unsorted linear search for now!
  for ( auto& e : db.cache ) {
    if ( hash == e.hash() && v == *e.dataShPtr() )
      return e.clone();
  }
  auto enew
    = EnergyGridHashPtrFact::create( std::make_shared<VectD>(std::move(v)),
                                     hash );
  db.cache.push_back(enew.clone());
  return enew;
}

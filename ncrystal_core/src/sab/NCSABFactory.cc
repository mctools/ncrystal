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

#include "NCrystal/internal/sab/NCSABFactory.hh"
#include "NCrystal/internal/sab/NCSABIntegrator.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/fact_utils/NCFactoryUtils.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace SAB {
    namespace {

      UniqueIDValue emptyEgridUIDVal()
      {
        static NC::UniqueID uid_empty;
        return uid_empty.getUniqueID();
      }

      //Cache key is (sabdata uid, egrid uid, sabdata ptr):
      //TODO: we should use new thin-key support instead of these
      //shared_obj<const NC::SABData>* pointers!
      typedef std::tuple<UniqueIDValue,
                         UniqueIDValue,
                         shared_obj<const NC::SABData>*> ScatHelperCacheKey;

      class ScatterHelperFactory
        : public NC::CachedFactoryBase<ScatHelperCacheKey,
                                       SABScatterHelper,
                                       20/*NStrongRefsKept*/> {
      public:
        const char* factoryName() const final { return "ScatterHelperFactory"; }
        std::string keyToString( const ScatHelperCacheKey& key ) const final
        {
          std::ostringstream ss;
          ss<<"(SABData id="<<std::get<0>(key).value
            <<";egrid id="<<std::get<1>(key).value<<")";
          return ss.str();
        }
      protected:
        virtual ShPtr actualCreate( const ScatHelperCacheKey& key ) const final
        {
          auto sabdata_shptr = *std::get<2>(key);
          nc_assert( sabdata_shptr->getUniqueID() == std::get<0>(key) );
          auto egrid_shptr = egridFromUniqueID(std::get<1>(key));
          return createScatterHelper(std::move(sabdata_shptr),
                                     std::move(egrid_shptr));
        }
      };

      ScatterHelperFactory& getScatterHelperFactory()
      {
        static ScatterHelperFactory s_scathelperfact;
        return s_scathelperfact;
      }

      //Egrid UID cache. NB: We never clean this particular cache in order to
      //preserve the id's, in case something somewhere is still hanging on to
      //one of them after a cache clearance.
      struct EgridUIDCacheDB {
        std::map< HashValue,
                  std::vector<std::pair<std::shared_ptr<const VectD>,
                                        UniqueID>>> egridHashCache;
        std::map< uint64_t, std::shared_ptr<const VectD>* > uid2egrid;
        std::mutex mtx;
      };

      EgridUIDCacheDB& getEgridUIDCacheDB()
      {
        static EgridUIDCacheDB db;
        return db;
      }
    }
  }
}

std::unique_ptr<const NC::SAB::SABScatterHelper>
NC::SAB::createScatterHelper( shared_obj<const NC::SABData> data,
                              std::shared_ptr<const VectD> energyGrid )
{
  nc_assert(!!data);
  SABIntegrator si(data,energyGrid.get());
  auto sh = si.createScatterHelper();
  return ncmake_unique<SABScatterHelper>(std::move(sh));
}

void NC::SAB::clearScatterHelperCache() {
  getScatterHelperFactory().cleanup();
}

NC::shared_obj<const NC::SAB::SABScatterHelper>
NC::SAB::createScatterHelperWithCache( shared_obj<const NC::SABData> dataptr,
                                       std::shared_ptr<const VectD> egrid )
{
  nc_assert_always(!!dataptr);

  ScatHelperCacheKey key( dataptr->getUniqueID(),
                          egridToUniqueID( egrid ),
                          &dataptr );

  return getScatterHelperFactory().create(key);
}

NC::UniqueIDValue NC::SAB::egridToUniqueID(const NC::VectD& egrid)
{
  if ( egrid.empty() ) {
    //empty => no need for hashing or locking
    return emptyEgridUIDVal();
  }

  //NB: code duplicated from here to following function
  auto hash = hashContainer(egrid);
  auto& db = getEgridUIDCacheDB();
  NCRYSTAL_LOCK_GUARD(db.mtx);
  auto& v = db.egridHashCache[hash];//In absence of hash collisions, v will have length 0 or 1.
  for (auto& e : v) {
    if ( *e.first == egrid )
      return e.second.getUniqueID();//exists in cache already
  }
  //Add new:
  v.emplace_back(std::make_shared<const VectD>(egrid),UniqueID() );
  auto uidval = v.back().second.getUniqueID();
  db.uid2egrid[uidval.value] = &v.back().first;
  return uidval;
}

NC::UniqueIDValue NC::SAB::egridToUniqueID(const std::shared_ptr<const NC::VectD>& egrid)
{
  if ( !egrid || egrid->empty() ) {
    //Treat nullptr as empty grid and handle specially.
    return emptyEgridUIDVal();
  }

  //code duplicated here from preceding function
  auto hash = hashContainer(*egrid);
  auto& db = getEgridUIDCacheDB();
  NCRYSTAL_LOCK_GUARD(db.mtx);
  auto& v = db.egridHashCache[hash];//In absence of hash collisions, v will have length 0 or 1.
  for (auto& e : v) {
    if ( *e.first == *egrid )
      return e.second.getUniqueID();//exists in cache already
  }
  //Add new:
  v.emplace_back( egrid, UniqueID() );
  auto uidval = v.back().second.getUniqueID();
  db.uid2egrid[uidval.value] = &v.back().first;
  return uidval;
}

std::shared_ptr<const NC::VectD>
NC::SAB::egridFromUniqueID( NC::UniqueIDValue uidval )
{
  if ( uidval == emptyEgridUIDVal() )
    return nullptr;
  auto& db = getEgridUIDCacheDB();
  NCRYSTAL_LOCK_GUARD(db.mtx);
  auto it = db.uid2egrid.find(uidval.value);
  if ( it == db.uid2egrid.end() )
    NCRYSTAL_THROW(LogicError,"egridFromUniqueID passed uid which was not"
                   " created by call to egridToUniqueID");
  return *it->second;
}

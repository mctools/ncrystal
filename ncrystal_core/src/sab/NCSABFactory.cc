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
#include "NCrystal/internal/sab/NCSABCfg.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace SAB {
    namespace {

      UniqueIDValue emptyEgridUIDVal()
      {
        static UniqueID uid_empty;
        return uid_empty.getUniqueID();
      }

      //Cache key is (sabdata uid, egrid uid, legacy opts). The thick key also
      //carries a reference to the SAB data.

      using EGridMargin = SABSampler::EGridMargin;

      struct LegacyOptsDecoded final {
        EGridMargin egridMargin = { EGridMargin::default_value };
        bool disable_betafix = false;
        LegacyOptsDecoded( int knllux )
          : LegacyOptsDecoded([knllux]()
          {
            nc_assert_always ( knllux
                               >= static_cast<int>(LegacySABAlgOpts::MIN)
                               && knllux
                               <= static_cast<int>(LegacySABAlgOpts::MAX) );
            return static_cast<LegacySABAlgOpts>(knllux);
          }())
        {
        }

        LegacyOptsDecoded( LegacySABAlgOpts opts )
        {
          static_assert( static_cast<int>(LegacySABAlgOpts::MAX) == -1, "" );
          static_assert( static_cast<int>(LegacySABAlgOpts::MIN) == -6, "" );
          switch( opts ) {
          case LegacySABAlgOpts::OVERSAMPLE10:
            egridMargin = EGridMargin{ 10.0 };
            break;
          case LegacySABAlgOpts::OVERSAMPLE50:
            egridMargin = EGridMargin{ 50.0 };
            break;
          case LegacySABAlgOpts::NOBETAFIX:
            disable_betafix = true;
            break;
          case LegacySABAlgOpts::NOBETAFIX_OVERSAMPLE10:
            egridMargin = EGridMargin{ 10.0 };
            disable_betafix = true;
            break;
          case LegacySABAlgOpts::NOBETAFIX_OVERSAMPLE50:
            egridMargin = EGridMargin{ 50.0 };
            disable_betafix = true;
            break;
          default:
          case LegacySABAlgOpts::DEFAULT:
            break;
          }
        }
      };

      using ScatHelperCacheKey_Thin = std::tuple<UniqueIDValue,
                                                 UniqueIDValue,
                                                 LegacySABAlgOpts>;

      struct ScatHelperCacheKey {
        ScatHelperCacheKey_Thin thin_key;
        std::shared_ptr<const SABData> sabdata_ptr;
      };

      struct ScatHelperCache_KeyThinner {
        using key_type = ScatHelperCacheKey;
        using thinned_key_type = ScatHelperCacheKey_Thin;
        template <class TMap>
        static typename TMap::mapped_type&
        cacheMapLookup( TMap& map, const key_type& key,
                        Optional<thinned_key_type>& tkey )
        {
          if ( !tkey.has_value() )
            tkey = key.thin_key;
          return map[tkey.value()];
        }
      };

      constexpr auto scathelperfact_nstrongrefskept = 20;

      class ScatterHelperFactory
        : public CachedFactoryBase<ScatHelperCacheKey,
                                   SABScatterHelper,
                                   scathelperfact_nstrongrefskept,
                                   ScatHelperCache_KeyThinner> {
      public:
        const char* factoryName() const final { return "ScatterHelperFactory"; }
        std::string keyToString( const ScatHelperCacheKey& key ) const final
        {
          std::ostringstream ss;
          ss<<"(SABData id="<<std::get<0>(key.thin_key).value
            <<";egrid id="<<std::get<1>(key.thin_key).value;
          if ( std::get<2>(key.thin_key) != LegacySABAlgOpts::DEFAULT ) {
            LegacyOptsDecoded o(std::get<2>(key.thin_key));
            LegacyOptsDecoded odef(LegacySABAlgOpts::DEFAULT);
            if ( o.egridMargin.value != odef.egridMargin.value )
              ss << ";egridmargin=" << o.egridMargin.value;
            if ( o.disable_betafix != odef.disable_betafix )
              ss << ";disablebetafix=" << int(o.disable_betafix);
          }
          ss<<")";
          return ss.str();
        }
      protected:
        virtual ShPtr actualCreate( const ScatHelperCacheKey& key ) const final
        {
          nc_assert( key.sabdata_ptr != nullptr );
          nc_assert( key.sabdata_ptr->getUniqueID() == std::get<0>(key.thin_key) );
          return createScatterHelper( key.sabdata_ptr,
                                      egridFromUniqueID(std::get<1>(key.thin_key)),
                                      std::get<2>(key.thin_key) );
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

      //New SABExtended (S.E.) cache, key is (sab uid, egrid uid, knllux):
      using SECacheKey_Thin = std::tuple<UniqueIDValue,UniqueIDValue,int>;
      struct SECacheKey {
        SECacheKey_Thin thin_key;
        std::shared_ptr<const SABData> sabdata_ptr;
      };

      struct SECache_KeyThinner {
        using key_type = SECacheKey;
        using thinned_key_type = SECacheKey_Thin;
        template <class TMap>
        static typename TMap::mapped_type&
        cacheMapLookup( TMap& map, const key_type& key,
                        Optional<thinned_key_type>& tkey )
        {
          if ( !tkey.has_value() )
            tkey = key.thin_key;
          return map[tkey.value()];
        }
      };

      constexpr auto factSE_nstrongrefskept = 20;

      class SEFactory final
        : public CachedFactoryBase<SECacheKey,
                                   SABUtils::SABExtended,
                                   factSE_nstrongrefskept,
                                   SECache_KeyThinner> {
      public:
        const char* factoryName() const override
        {
          return "SABExtendedFactory";
        }

        std::string keyToString( const SECacheKey& key ) const final
        {
          std::ostringstream ss;
          ss<<"(SABData id="<<std::get<0>(key.thin_key).value
            <<";egrid id="<<std::get<1>(key.thin_key).value
            <<";knllux="<<std::get<2>(key.thin_key)<<")";
          return ss.str();
        }
      protected:
        virtual ShPtr actualCreate( const SECacheKey& key ) const final
        {
          nc_assert( key.sabdata_ptr != nullptr );
          nc_assert( key.sabdata_ptr->getUniqueID()
                     == std::get<0>(key.thin_key) );
          int knllux = std::get<2>(key.thin_key);
          auto egrid = egridFromUniqueID(std::get<1>(key.thin_key));
          return createSABExtendedNoCache( knllux,
                                           key.sabdata_ptr,
                                           std::move(egrid) );
        }
      };

      SEFactory& getSEFactory()
      {
        static SEFactory s_fact;
        return s_fact;
      }

    }
  }
}

std::unique_ptr<const NC::SAB::SABScatterHelper>
NC::SAB::createScatterHelper( shared_obj<const SABData> data,
                              std::shared_ptr<const VectD> energyGrid,
                              LegacySABAlgOpts opts_raw)
{
  nc_assert(!!data);
  LegacyOptsDecoded opts(opts_raw);
  SABIntegrator si( data,
                    energyGrid.get(),
                    nullptr/*default extender*/,
                    opts.egridMargin,
                    opts.disable_betafix );
  auto sh = si.createScatterHelper();
  return ncmake_unique<SABScatterHelper>(std::move(sh));
}

NC::shared_obj<const NC::SAB::SABScatterHelper>
NC::SAB::createScatterHelperWithCache( shared_obj<const SABData> sabdataptr,
                                       std::shared_ptr<const VectD> egrid,
                                       LegacySABAlgOpts opts )
{
  ScatHelperCacheKey key;
  std::get<0>(key.thin_key) = sabdataptr->getUniqueID();
  std::get<1>(key.thin_key) = egridToUniqueID(egrid);
  std::get<2>(key.thin_key) = opts;
  key.sabdata_ptr = std::move(sabdataptr);
  return getScatterHelperFactory().create(key);
}

NC::UniqueIDValue NC::SAB::egridToUniqueID(const VectD& egrid)
{
  if ( egrid.empty() ) {
    //empty => no need for hashing or locking
    return emptyEgridUIDVal();
  }

  //NB: code duplicated from here to following function
  auto hash = hashContainer(egrid);
  auto& db = getEgridUIDCacheDB();
  NCRYSTAL_LOCK_GUARD(db.mtx);
  auto& v = db.egridHashCache[hash];//In absence of hash collisions,
                                    //v will have length 0 or 1.
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

NC::UniqueIDValue
NC::SAB::egridToUniqueID(const std::shared_ptr<const VectD>& egrid)
{
  if ( !egrid || egrid->empty() ) {
    //Treat nullptr as empty grid and handle specially.
    return emptyEgridUIDVal();
  }

  //code duplicated here from preceding function
  auto hash = hashContainer(*egrid);
  auto& db = getEgridUIDCacheDB();
  NCRYSTAL_LOCK_GUARD(db.mtx);
  auto& v = db.egridHashCache[hash];//In absence of hash collisions,
                                    //v will have length 0 or 1.
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
NC::SAB::egridFromUniqueID( UniqueIDValue uidval )
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

NC::shared_obj<const NC::SABUtils::SABExtended>
NC::SAB::createSABExtendedNoCache( int knllux,
                                   shared_obj<const SABData> sab,
                                   std::shared_ptr<const VectD> egrid )
{
  auto cfg = SABCfg::createConfig(knllux);
  return SABUtils::SABExtended::createWithFGExtender( cfg,
                                                      std::move(sab),
                                                      std::move(egrid) );
}

NC::shared_obj<const NC::SABUtils::SABExtended>
NC::SAB::createSABExtendedWithCache( int knllux,
                                     shared_obj<const SABData> sab,
                                     std::shared_ptr<const VectD> egrid )
{
  SECacheKey key;
  std::get<0>(key.thin_key) = sab->getUniqueID();
  std::get<1>(key.thin_key) = egridToUniqueID(egrid);
  std::get<2>(key.thin_key) = knllux;
  key.sabdata_ptr = std::move(sab);
  return getSEFactory().create(key);
}

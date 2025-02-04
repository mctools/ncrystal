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
#include "NCrystal/internal/utils/NCFileUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace FactImpl {

    template<std::size_t NCACHED>
    class TDProdDB {
    public:
      TextDataSP produceTextDataSP_PreferPreviousObject( TextData&& newtd )
      {
        nc_assert(newtd.dataUID().isUnset());
        uint64_t checkSum = newtd.rawData().calcCheckSum();

        auto it = m_db.begin();
        auto itE = m_db.end();

        //First check if we have a compatible object already:
        for ( ; it != itE; ++it ) {
          if ( it->first == checkSum
               && newtd.hasIdenticalMetaData(it->second)
               && newtd.rawData().hasSameContent( it->second->rawData() ) )
            {
              break;//Found!
            }
        }
        if ( it == itE ) {
          //Did not find existing. Finalize new TextData object construction by
          //actually assigning a UID, then insert in cache:
          auto newtdsp = makeSO<const TextData>(TextData::internal_consumeAndSetNewUID(std::move(newtd)));

          if ( m_db.size() == NCACHED ) {
            //Must first discard existing entry to make room:
            auto itLast = std::prev(m_db.end());
            for ( auto it2 = m_db.begin(); it2 != itLast; ++it2 )
              *it2 = std::move( *std::next(it2) );
            m_db.pop_back();
          }
          //add to back (which is the correct position indicating most recent
          //access):
          m_db.push_back(std::pair<uint64_t,TextDataSP>(checkSum,newtdsp));
          return newtdsp;
        } else {
          //Found existing! But before returning we need to reorder entries,
          //moving the requested entry to the back and sliding the other ones
          //down to make room:
          TextDataSP result(std::move(it->second));
          auto itLast = std::prev(m_db.end());
          for ( auto it2 = it; it2 != itLast; ++it2 )
            *it2 = std::move( *std::next(it2) );
          *itLast = std::pair<uint64_t,TextDataSP>(checkSum,result);
          return result;
        }
      }
      void clear() { m_db.clear(); }
    private:
      SmallVector<std::pair<uint64_t,TextDataSP>,NCACHED> m_db;
    };

    class TDProd {
      static constexpr std::size_t small_threshold_nBytes =   200000;//200kB
      static constexpr std::size_t large_threshold_nBytes = 10000000;//10MB
#ifndef NCRYSTAL_ALLOW_ULTRA_LARGE_FILES
      static constexpr std::size_t very_large_threshold_nBytes =  500000000;//500MB (max allowed size)
#else
      static constexpr std::size_t very_large_threshold_nBytes = 1000000000000000ull;//1000TB (!!)
#endif
      TDProdDB<200> m_db_small;//in hypothetical extreme case this can retain 100kB*200 = 20MB
      TDProdDB<10> m_db_large;//in hypothetical extreme case this can retain 10MB*10 = 100MB
      TDProdDB<3> m_db_veryLarge;//max kept is 500MB*3=1.5GB (but the user really asked for it!)
    public:
      void clear()
      {
        m_db_small.clear();
        m_db_large.clear();
        m_db_veryLarge.clear();
      }

      static TextData produceTextDataWithoutCache( const TextDataPath& textdatapath, TextDataSource&& tds )
      {
        std::string dataType = tds.dataType();
        Variant<std::string,RawStrData> data( std::move(tds).data() );
        Optional<TextData::LastKnownOnDiskAbsPath> lastKnownOnDiskPath;
        Optional<RawStrData> rawdata;

        //For consistency, simplicity, and to avoid spuriosly duplicated cached
        //entries, the DataSourceName name is always just the basename unless
        //explicitly chosen otherwise:
        std::string dataSourceName = tds.suggestedDataSourceName();
        if ( dataSourceName.empty() )
          dataSourceName = basename(textdatapath.path());

        if ( data.has_value<std::string>() ) {
          //On-disk path.
          std::string path = std::move(data.get<std::string>());

          //Convert to absolute path if relative:
          if (!path_is_absolute(path))
            path = path_join(ncgetcwd(),path);

          //Normalise if possible (canonical path, no symlinks, etc.):
          auto pn = tryRealPath( path );
          if ( !pn.empty() )
            path = std::move(pn);

          lastKnownOnDiskPath = TextData::LastKnownOnDiskAbsPath{path};
          Optional<std::string> content = readEntireFileToString( lastKnownOnDiskPath.value().value );
          if ( !content.has_value() )
            NCRYSTAL_THROW2(DataLoadError,"Missing or unreadable file: "<<path);
          rawdata = RawStrData(std::move(content.value()));

        } else {
          nc_assert( data.has_value<RawStrData>() );
          rawdata = std::move(data.get<RawStrData>());
        }

        nc_assert(rawdata.has_value());

        if ( dataType.empty() )
          dataType = FactImpl::guessDataType( rawdata.value(), textdatapath.path() );

        return TextData( TextData::internal_with_unset_textdatauid_t{},//no need to consume uid before we actually return to consumer
                         std::move(rawdata.value()),
                         TextData::DataType{dataType},
                         DataSourceName{dataSourceName},
                         std::move(lastKnownOnDiskPath) );
      }

      TextDataSP produceTextDataSP_PreferPreviousObject( TextData&& td )
      {
        std::size_t size = std::distance(td.rawData().begin(),td.rawData().end());
        if ( size <= small_threshold_nBytes ) {
          return m_db_small.produceTextDataSP_PreferPreviousObject(std::move(td));
        } else if ( size <= large_threshold_nBytes ) {
          return m_db_large.produceTextDataSP_PreferPreviousObject(std::move(td));
        } else if ( size <= very_large_threshold_nBytes ) {
          return m_db_veryLarge.produceTextDataSP_PreferPreviousObject(std::move(td));
        } else {
#ifndef NCRYSTAL_ALLOW_ULTRA_LARGE_FILES
          const char * extraguidance
            = " [NB: Recompile NCrystal with NCRYSTAL_ALLOW_ULTRA_LARGE_FILES to increase limit]";
#else
          const char * extraguidance = "";
#endif
          NCRYSTAL_THROW2(DataLoadError,"Input has unsupported data size ("<<(size*1e-6)
                          <<"MB, max allowed is "<<(very_large_threshold_nBytes*1e-6)
                          <<"MB): " << td.dataSourceName() << extraguidance );
          return OptionalTextDataSP{nullptr};//dummy
        }
      }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //Global DB:
    struct GlobalTDProd {
      TDProd db;
      std::mutex mtx;
    };
    GlobalTDProd& globalTDProd() { static GlobalTDProd db; return db; }
    void clearGlobalTDProdCache() {
      auto& db = globalTDProd();
      NCRYSTAL_LOCK_GUARD(db.mtx);
      db.db.clear();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    //Functions fwd declared and used in FactImpl.cc:
    TextDataSP produceTextDataSP_PreferPreviousObject( const TextDataPath& path, TextDataSource&& tds)
    {
      TextData td = TDProd::produceTextDataWithoutCache( path, std::move(tds) );
      nc_assert(td.dataUID().isUnset());//unset, if we don't actually need it we should ideally not consume an UID
      auto& db = globalTDProd();
      NCRYSTAL_LOCK_GUARD(db.mtx);
      static bool first = true;
      if ( first ) {
        first = false;
        registerCacheCleanupFunction(clearGlobalTDProdCache);
      }
      return db.db.produceTextDataSP_PreferPreviousObject( std::move(td) );
    }

  }
}

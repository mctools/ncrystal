#ifndef NCrystal_FactoryUtils_hh
#define NCrystal_FactoryUtils_hh

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

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/core/NCSmallVector.hh"
#include "NCrystal/internal/utils/NCMsg.hh"

namespace NCRYSTAL_NAMESPACE {

  constexpr unsigned CachedFactory_KeepAllStrongRefs = std::numeric_limits<unsigned>::max();

  template<class TKey>
  struct CFB_Unthinned_t {
    //Default key thinning strategy is to not actually do any thinning. If
    //implementing a cache which uses thinning, the Optional<thinned_key_type>&
    //arguments can be used as a place to cache the thinned key (to avoid
    //thinning again and again).
    using key_type = TKey;
    using thinned_key_type = TKey;
    template <class TMap>
    static typename TMap::mapped_type& cacheMapLookup( TMap& map,
                                                       const key_type& key,
                                                       Optional<thinned_key_type>& )
    {
      return map[key];
    }
  };

  template< class TKey,
            class TValue,
            unsigned NStrongRefsKept = 5,
            class TKeyThinner = CFB_Unthinned_t<TKey>>
  class CachedFactoryBase {
  public:

    /////////////////////////////////////////////////////////////////////////////////
    //
    // Base class for multi-thread safe factories. Client code must implement
    // the actualCreate function, plus a few cosmetic functions. Factories will
    // always keep weak references (weak_ptr) to the created objects, but can be
    // configured with the NStrongRefsKept parameter to also keep strong
    // references (shared_ptr) to them (thus ensuring the created objects are
    // kept alive and won't have to be recreated later). Strong refs are kept to
    // the last NStrongRefsKept objects accessed.
    //
    // The strong refs and any dangling weak refs can be cleaned up by an explicit
    // call to the cleanup function of a particular factory, or simply by calling the
    // global clearCaches function which will in turn call the cleanup function of
    // all factories.
    //
    /////////////////////////////////////////////////////////////////////////////////

    using key_type = TKey;
    using thinned_key_type = typename TKeyThinner::thinned_key_type;
    typedef TValue value_type;
    typedef std::shared_ptr<const value_type> ShPtr;
    typedef std::weak_ptr<const value_type> WeakPtr;

    ShPtr create(const key_type&);

    virtual std::string keyToString( const key_type& ) const = 0;
    virtual const char* factoryName() const = 0;

    //Cleanup cache, release all kept strong and weak references (this function
    //is automatically registered with and invoked by the global clearCaches
    //function). But note that it is NOT called from the factory destructor,
    //since that can trigger memory errors during programme shutdown.:
    void cleanup();

    //To automatically call a function whenever cleanup() is invoked:
    void registerCleanupCallback(voidfct_t);

    //Statistics about the current cache:
    struct Stats { std::size_t nstrongrefs, nweakrefs; };
    Stats currentStats();

    //NB: This next might seem sensible, but gives troubles since most
    //CacheFactoryBase instances are kept as global static objects, with
    //undefined destruction order: ~CachedFactoryBase() { cleanup(); }

  protected:
    virtual ShPtr actualCreate(const key_type&) const = 0;
  private:
    struct CacheEntry {
      WeakPtr weakPtr;
      unsigned nwork = 0;
      unsigned ncleanup = 0;
    };
    using CacheMap = std::map<thinned_key_type,CacheEntry>;
    CacheMap m_cache;
    std::mutex m_mutex;
    unsigned m_ncleanup = 0;
    class StrongRefKeeper;
    StrongRefKeeper m_strongRefs;
    bool m_cleanupNeedsRegistry = true;
    SmallVector<voidfct_t,1,SVMode::LOWFOOTPRINT> m_cleanupCallbacks;
  };

  ///////////////////////////////////////////////////////////////////////////
  // Global control affecting all CachedFactoryBase instances (the default //
  // settings of these can also be controlled by environment variables).   //
  ///////////////////////////////////////////////////////////////////////////

  void enableFactoryVerbosity( bool status = true );
  bool getFactoryVerbosity();

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace thread_details {
    std::string currentThreadIDForPrint();
  }

  template<class TKey, class TValue, unsigned NStrongRefsKept,class TKT>
  class CachedFactoryBase<TKey,TValue,NStrongRefsKept,TKT>::StrongRefKeeper {

    std::vector<ShPtr> m_v;

    static void reserveCapacity( std::vector<ShPtr>& v )
    {
      if (NStrongRefsKept>0)
        v.reserve(NStrongRefsKept > 512 ? 512 : NStrongRefsKept);
    }

  public:
    StrongRefKeeper()  { reserveCapacity(m_v); }
    void clear() { m_v.clear(); }
    bool empty() const { return m_v.empty(); }
    std::size_t size() const { return static_cast<std::size_t>(m_v.size()); }
    void wasAccessedAndIsNotInList( const ShPtr& sp ) {
      //was just created, so we know it is not already in the list.
      if (NStrongRefsKept==0)
        return;
      if (NStrongRefsKept==CachedFactory_KeepAllStrongRefs) {
        m_v.push_back( sp );
        return;
      }
      if ( m_v.size()==NStrongRefsKept ) {
        //make room, discard first entry (was accessed longest ago).
        for ( auto i : ncrange(typename decltype(m_v)::size_type(1),m_v.size() ) )
          m_v[i-1] = std::move(m_v[i]);
        m_v.pop_back();
      }
      m_v.push_back( sp );
    }

    void wasAccessed( const ShPtr& sp ) {
      //was accessed, might (or might not) already be in the list.
      //NStrongRefsKept==CachedFactory_KeepAllStrongRefs
      if (NStrongRefsKept==0)
        return;
      //Check if we already have it:
      auto it = m_v.begin();
      auto itE = m_v.end();
      for (;it!=itE;++it) {
        if ( *it == sp ) {
          //Already there! Shuffle it to the end:
          if ( std::next(it) == itE )
            return;//already at end
          auto itLast = std::prev(itE);
          for ( ;it!=itLast; ++it)
            *it = std::move(*std::next(it));
          *itLast = sp;
          return;
        }
      }
      //was not already in the list:
      wasAccessedAndIsNotInList(sp);
    }
  };

  template<class TKey,class TValue,unsigned N,class TKT>
  inline void CachedFactoryBase<TKey,TValue,N,TKT>::registerCleanupCallback(voidfct_t fn)
  {
    NCRYSTAL_LOCK_GUARD(m_mutex);
    m_cleanupCallbacks.push_back(fn);
  }

  template<class TKey,class TValue,unsigned N,class TKT>
  inline void CachedFactoryBase<TKey,TValue,N,TKT>::cleanup()
  {
    NCRYSTAL_LOCK_GUARD(m_mutex);
    ++m_ncleanup;
    m_strongRefs.clear();
    m_cache.clear();
    for ( const auto& fn : m_cleanupCallbacks )
      fn();
  }

  template<class TKey,class TValue,unsigned N,class TKT>
  inline typename CachedFactoryBase<TKey,TValue,N,TKT>::Stats CachedFactoryBase<TKey,TValue,N,TKT>::currentStats()
  {
    NCRYSTAL_LOCK_GUARD(m_mutex);
    Stats s;
    s.nstrongrefs = m_strongRefs.size();
    s.nweakrefs = static_cast<std::size_t>(m_cache.size());
    return s;
  }

  template<class TKey,class TValue,unsigned NStrongRefsKept,class TKT>
  inline std::shared_ptr<const TValue> CachedFactoryBase<TKey,TValue,NStrongRefsKept,TKT>::create(const TKey& key)
  {
    const bool verbose = getFactoryVerbosity();
    const std::string keystr = ( verbose ? keyToString(key) : std::string() );
    Optional<thinned_key_type> thinned_key;

    {
      NCRYSTAL_LOCK_GUARD(m_mutex);

      if ( m_cleanupNeedsRegistry ) {
        //Unrelated, but needs to be done while m_mutex is locked.
        m_cleanupNeedsRegistry = false;
        voidfct_t fct_cleanup = [this](){ this->cleanup(); };
        registerCacheCleanupFunction(fct_cleanup);
      }

      if ( verbose )
        NCRYSTAL_MSG(this->factoryName()
                     <<" (thread_"<<thread_details::currentThreadIDForPrint()<<")"
                     <<" : Request to provide object for key "<<keystr);

      auto& cache_entry = TKT::cacheMapLookup( m_cache, key, thinned_key );
      ShPtr res_existing = cache_entry.weakPtr.lock();
      if ( res_existing != nullptr ) {
        //Already there!
        if ( verbose )
          NCRYSTAL_MSG(this->factoryName()
                       <<" (thread_"<<thread_details::currentThreadIDForPrint()<<")"
                       <<" : Return pre-existing cached object for key "<<keystr);
        //Record access:
        m_strongRefs.wasAccessed( res_existing );
        return res_existing;
      }

      cache_entry.ncleanup = m_ncleanup;
      cache_entry.nwork += 1;

      //guard against cyclic dependencies:
      constexpr unsigned nwork_limit = 50;
      if ( cache_entry.nwork > nwork_limit ) {
        //Almost certainly a cyclic dependency. We could in principle get a
        //false positive in case a user would use a huge amount of threads to
        //simultaneously make the same request. For now we simply ignore that
        //possibility.
        //
        //NB: Tried nwork_limit=1000 but it resulted in segfaults. So lowered it
        //drastically.
        NCRYSTAL_THROW(BadInput,"Cyclic dependency in factory request"
                       " detected (check your input configurations"
                       " and data for cyclic references)!");
      }
    }//release m_mutex lock

    //Not in cache already, go ahead and construct it (without holding any
    //lock):

    if ( verbose )
      NCRYSTAL_MSG(this->factoryName()
                   <<" (thread_"<<thread_details::currentThreadIDForPrint()<<")"
                   << " : Creating (from scratch) object for key " << keystr);

    ShPtr res = actualCreate(key);

    {
      //reaquire lock, and populate result (or discard if another thread beat us
      //to it):
      NCRYSTAL_LOCK_GUARD(m_mutex);

      if ( verbose )
        NCRYSTAL_MSG(this->factoryName()
                     <<" (thread_"<<thread_details::currentThreadIDForPrint()<<")"
                     <<" : Finished construction");

      auto& cache_entry = TKT::cacheMapLookup( m_cache, key, thinned_key );
      --cache_entry.nwork;//no matter what, we no longer work on it.

      //Populate cache unless we were beaten to it:
      ShPtr res_existing = cache_entry.weakPtr.lock();
      if ( res_existing != nullptr ) {
        if ( verbose )
          NCRYSTAL_MSG(this->factoryName()
                       <<" (thread_"<<thread_details::currentThreadIDForPrint()<<")"
                       <<" : Finished construction but another thread beat us to it.");
        res.reset();//discard our own result, always return the first recorded
        m_strongRefs.wasAccessed( res_existing );
        return res_existing;
      }

      //Record out result, unless we got intercepted by a cleanup() call:
      if ( cache_entry.ncleanup == m_ncleanup ) {
        cache_entry.weakPtr = res;
        m_strongRefs.wasAccessedAndIsNotInList( res );
        return res;
      }
    }
    //We must have gotten intercepted by a cleanup() call, so let us try again:
    return this->create( std::move(key) );
  }
}

#endif

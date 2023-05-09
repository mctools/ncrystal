#ifndef NCrystal_FactoryUtils_hh
#define NCrystal_FactoryUtils_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2023 NCrystal developers                                   //
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

#include "NCrystal/NCDefs.hh"
#include "NCrystal/NCSmallVector.hh"
#include <chrono>
#include <iostream>
#ifndef NCRYSTAL_DISABLE_THREADS
#  include <thread>
#endif

namespace NCrystal {

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
    static typename TMap::mapped_type& cacheMapLookup( TMap& map, const key_type& key, Optional<thinned_key_type>& )
    {
      return map[key];
    }
  };

  template< class TKey, class TValue, unsigned NStrongRefsKept = 5, class TKeyThinner = CFB_Unthinned_t<TKey>>
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
    ShPtr createWithoutCache(const key_type&) const;

    virtual std::string keyToString( const key_type& ) const = 0;
    virtual const char* factoryName() const = 0;

    //Cleanup cache, release all kept strong and weak references (this function
    //automatically registered with and invoked by global clearCaches
    //function). But note that it is NOT called from the factory destructor,
    //since that can trigger memory errors during programme shutdown.:
    void cleanup();

    //To automatically call a function whenever cleanup() is invoked:
    void registerCleanupCallback(std::function<void()>);

    struct Stats { std::size_t nstrongrefs, nweakrefs; };
    Stats currentStats();

    //NB: This might seem sensible, but gives troubles since most
    //CacheFactoryBase instances are kept as global static objects, with
    //undefined destruction order: ~CachedFactoryBase() { cleanup(); }

  protected:
    virtual ShPtr actualCreate(const key_type&) const = 0;
  private:
    struct CacheEntry {
      bool underConstruction = false;
      bool wasInvalidatedDuringConstruction = false;
      WeakPtr weakPtr;
      void clear() { underConstruction = wasInvalidatedDuringConstruction = false; weakPtr.reset(); }
    };
    using CacheMap = std::map<thinned_key_type,CacheEntry>;
    CacheMap m_cache;
    std::mutex m_mutex;
    class StrongRefKeeper;
    StrongRefKeeper m_strongRefs;
    bool m_cleanupNeedsRegistry = true;
    SmallVector<std::function<void()>,1,SVMode::LOWFOOTPRINT> m_cleanupCallbacks;
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

namespace NCrystal {

  namespace thread_details {
#ifndef NCRYSTAL_DISABLE_THREADS
    inline std::thread::id currentThreadIDForPrint() { return std::this_thread::get_id(); }
#else
    inline constexpr const char * currentThreadIDForPrint() noexcept { return "<thread-id-unavailable>"; }
#endif
  }
  namespace detail {
#ifndef NCRYSTAL_DISABLE_THREADS
    void registerThreadWork(std::thread::id);
    void registerThreadWorkDone(std::thread::id);
    void registerThreadAsWaiting(std::thread::id);
    void registerThreadAsFinishedWaiting(std::thread::id);
#endif
  }

  template<class TKey, class TValue, unsigned NStrongRefsKept,class TKT>
  class CachedFactoryBase<TKey,TValue,NStrongRefsKept,TKT>::StrongRefKeeper {
    std::vector<ShPtr> m_v;
    void reserveCapacity( std::vector<ShPtr>& v )
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

    void releaseOne( const ShPtr& sp ) {
      if ( !sp || empty() )
        return;
      auto itE = m_v.end();
      auto it = std::find( m_v.begin(), itE, sp );
      if ( it == itE )
        return;//not there
      //Move it to the end and pop_back:
      auto itLast = std::prev(itE);
      for ( ;it!=itLast; ++it)
        *it = std::move(*std::next(it));
      m_v.pop_back();
    }

    void releaseMany( const std::set<ShPtr>& to_remove )
    {
      if ( to_remove.size() <= 1 ) {
        releaseOne( *to_remove.begin() );
        return;
      }
      decltype(m_v) new_v;
      reserveCapacity(new_v);
      for ( auto& sp : m_v ) {
        if ( !to_remove.count(sp) )
          new_v.push_back( std::move(sp) );
      }
      std::swap(new_v,m_v);
    }
  };

  template<class TKey,class TValue,unsigned N,class TKT>
  inline void CachedFactoryBase<TKey,TValue,N,TKT>::registerCleanupCallback(std::function<void()> fn)
  {
    NCRYSTAL_LOCK_GUARD(m_mutex); m_cleanupCallbacks.push_back(fn);
  }

  template<class TKey,class TValue,unsigned N,class TKT>
  inline void CachedFactoryBase<TKey,TValue,N,TKT>::cleanup()
  {
    NCRYSTAL_LOCK_GUARD(m_mutex);
    m_strongRefs.clear();
    auto it = m_cache.begin();
    auto itE = m_cache.end();
    while (it!=itE) {
      auto itNext = std::next(it);
      if ( it->second.underConstruction ) {
        it->second.wasInvalidatedDuringConstruction = true;
      } else {
        m_cache.erase(it);
      }
      it = itNext;
    }
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

  template<class TKey,class TValue,unsigned N,class TKT>
  inline std::shared_ptr<const TValue> CachedFactoryBase<TKey,TValue,N,TKT>::createWithoutCache( const TKey& key ) const
  {
    if ( getFactoryVerbosity() )
      std::cout<< this->factoryName()
               <<" (thread_"<<thread_details::currentThreadIDForPrint()<<")"
               <<" : Request to provide object for key "<<keyToString(key)<<" (without cache)"<<std::endl;
    return actualCreate(key);
  }

  template<class TKey,class TValue,unsigned NStrongRefsKept,class TKT>
  inline std::shared_ptr<const TValue> CachedFactoryBase<TKey,TValue,NStrongRefsKept,TKT>::create(const TKey& key)
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    class Guard {
      //Local guard class. Kind of like std::lock_guard, but can remove set
      //underConstruction flag as well.
      std::mutex& m_mutex;
      bool* m_constructFlag = nullptr;
      bool m_isLocked = false;
    public:
#ifdef NCRYSTAL_DEBUG_LOCKS
#  define NCRYSTAL_DEBUG_LOCKS_ARGS __FILE__,__LINE__
#else
#  define NCRYSTAL_DEBUG_LOCKS_ARGS
#endif
      constexpr bool isLocked() const noexcept { return m_isLocked; }

      bool weHoldConstructFlag() const { return m_constructFlag != nullptr; }
      Guard( std::mutex& mutex ) : m_mutex(mutex) {}
      void releaseConstructFlagWithoutSetting() {
        m_constructFlag = nullptr;
      }
      void setConstructFlagFalseAndRelease() {
        if ( m_constructFlag != nullptr ) {
          ensureLock(NCRYSTAL_DEBUG_LOCKS_ARGS);
          *m_constructFlag = false;
          m_constructFlag = nullptr;
        }
      }
      void setConstructFlagWithGuard( bool* constructFlag ) {
        nc_assert(constructFlag);
        ensureLock(NCRYSTAL_DEBUG_LOCKS_ARGS);
        m_constructFlag = constructFlag;
        *m_constructFlag = true;
      }
      ~Guard() {
        setConstructFlagFalseAndRelease();
        ensureUnlock(NCRYSTAL_DEBUG_LOCKS_ARGS);
      }
#ifdef NCRYSTAL_DEBUG_LOCKS
      void ensureLock(const char* file,unsigned lineno)
#else
      void ensureLock()
#endif
      {
        if (!m_isLocked) {
#ifdef NCRYSTAL_DEBUG_LOCKS
          std::cout<<"NCrystal::FactoryUtils:: About to lock mutex due to Guard::ensureLock() from "<<file<<" : "<<lineno<<std::endl;
#endif
          NCRYSTAL_LOCK_MUTEX(m_mutex);
          m_isLocked = true;
        }
      }
#ifdef NCRYSTAL_DEBUG_LOCKS
      void ensureUnlock(const char* file,unsigned lineno)
#else
      void ensureUnlock()
#endif
      {
        if (m_isLocked) {
#ifdef NCRYSTAL_DEBUG_LOCKS
          std::cout<<"NCrystal::FactoryUtils:: About to unlock mutex due to Guard::ensureUnlock() from "<<file<<" : "<<lineno<<std::endl;
#endif
          NCRYSTAL_UNLOCK_MUTEX(m_mutex);
          m_isLocked = false;
        }
      }
    };

    //////////////////////////////////////////////////////////////////////////////////////
    const bool verbose = getFactoryVerbosity();
    const std::string keystr = ( verbose ? keyToString(key) : std::string() );

    Guard guard(m_mutex);
    guard.ensureLock(NCRYSTAL_DEBUG_LOCKS_ARGS);
    if ( m_cleanupNeedsRegistry ) {
      m_cleanupNeedsRegistry = false;
      std::function<void()> fct_cleanup = [this](){ this->cleanup(); };
      registerCacheCleanupFunction(fct_cleanup);
    }

    if ( verbose )
      std::cout<< this->factoryName()
               <<" (thread_"<<thread_details::currentThreadIDForPrint()<<")"
               <<" : Request to provide object for key "<<keystr<<std::endl;


    Optional<thinned_key_type> thinned_key;

    auto& cache_entry = TKT::cacheMapLookup( m_cache, key, thinned_key );
    ShPtr res = cache_entry.weakPtr.lock();
    if (!!res) {
      if ( verbose )
        std::cout<< this->factoryName()
                 <<" (thread_"<<thread_details::currentThreadIDForPrint()<<")"
                 <<" : Return pre-existing cached object for key "<<keystr<<std::endl;
      nc_assert_always(!cache_entry.underConstruction);

      //Record access:
      nc_assert(guard.isLocked());
      m_strongRefs.wasAccessed( res );
      return res;//easy: already there
    }
    //Not there: check if already under construction or if we should construct:

    if (!cache_entry.underConstruction) {
      guard.setConstructFlagWithGuard(&cache_entry.underConstruction);
      nc_assert(guard.weHoldConstructFlag());
    }

    //Have to either construct ourselves or wait for other thread to do it. In
    //both cases, stop holding mutex lock:

    guard.ensureUnlock(NCRYSTAL_DEBUG_LOCKS_ARGS);

    if (guard.weHoldConstructFlag()) {
      if ( verbose )
        std::cout << this->factoryName()
                  <<" (thread_"<<thread_details::currentThreadIDForPrint()<<")"
                  << " : Creating (from scratch) object for key " << keystr << std::endl;
      //Invoke actual creation function without holding the mutex lock.
      {
#ifndef NCRYSTAL_DISABLE_THREADS
        struct IsWorkingGuard {
          std::thread::id m_id;
          IsWorkingGuard() : m_id(std::this_thread::get_id()) { detail::registerThreadWork(m_id); }
          ~IsWorkingGuard() { detail::registerThreadWorkDone(m_id); }
        } isworkingguard;
#endif
        res = actualCreate(key);
      }
      //Populate result while holding mutex lock:
      guard.ensureLock(NCRYSTAL_DEBUG_LOCKS_ARGS);
      cache_entry = TKT::cacheMapLookup( m_cache, key, thinned_key );//reacquire after getting lock back
      //no one else should have tried to create this:
      nc_assert_always(cache_entry.underConstruction);
      nc_assert_always(!cache_entry.weakPtr.lock());
      //Check if was invalidated while constructing:
      if ( cache_entry.wasInvalidatedDuringConstruction ) {
        //oups, we were invalidated. Throw away result and restart:
        if ( verbose )
          std::cout<< this->factoryName()
                   <<" (thread_"<<thread_details::currentThreadIDForPrint()<<")"
                   <<" : Throwing away constructed result due to invalidation from another thread"<<std::endl;
        guard.releaseConstructFlagWithoutSetting();
        cache_entry.clear();
        guard.ensureUnlock(NCRYSTAL_DEBUG_LOCKS_ARGS);
        return this->create(key);
      } else {
        //all ok, populate cache::
        if ( verbose )
          std::cout<< this->factoryName()
                   <<" (thread_"<<thread_details::currentThreadIDForPrint()<<")"
                   <<" : Finished construction"<<std::endl;
        cache_entry.weakPtr = res;
        m_strongRefs.wasAccessedAndIsNotInList( res );
        guard.setConstructFlagFalseAndRelease();
        guard.ensureUnlock(NCRYSTAL_DEBUG_LOCKS_ARGS);
        return res;
      }
    } else {
      //Wait for other thread to populate cache. Sleep and recheck periodically.
#ifndef NCRYSTAL_DISABLE_THREADS
        struct IsWaitingGuard {
          std::thread::id m_id;
          IsWaitingGuard() : m_id(std::this_thread::get_id()) { detail::registerThreadAsWaiting(m_id); }
          ~IsWaitingGuard() { detail::registerThreadAsFinishedWaiting(m_id); }
        } iswaitingguard;
#endif
      while (true) {
#ifndef NCRYSTAL_DISABLE_THREADS
        //Try to detect cyclic dependencies:
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
#else
        NCRYSTAL_THROW(LogicError,"Other thread seems to be doing work (or you have"
                       " a cyclical dependency!!) - but NCrystal was built with "
                       "NCRYSTAL_DISABLE_THREADS and can not support this.");
#endif
        guard.ensureLock(NCRYSTAL_DEBUG_LOCKS_ARGS);
        cache_entry = TKT::cacheMapLookup( m_cache, key, thinned_key );//reacquire after getting lock back
        if ( verbose )
          std::cout<< this->factoryName()
                   <<" (thread_"<<thread_details::currentThreadIDForPrint()<<")"
                   <<" : Waiting for other thread to create (from scratch) object for key "<<keystr<<std::endl;
        res = cache_entry.weakPtr.lock();
        if (!!res)
          return res;//success! [other thread just created it and put it in m_strongRefs, probably wasteful to do it again]

        //Not yet. Double-check other thread is still trying:
        if (!cache_entry.underConstruction) {
          //Not there and no other thread is currently trying to construct
          //it. Technically we can't know if this situation happened because the
          //other thread ended prematurely (e.g. exception thrown), or because
          //it succeeded, but the result was already used and cleaned up again
          //due to all shared pointers going out of scope (which can happen even
          //if keeping strong refs, since cleanup() might have been invoked). Or
          //it could be because it was invalidated and the other thread didnt
          //restart creation yet. The best we can do is to try again, and hope
          //this time there won't be any concurrent attempts to create. In case
          //of an underlying error condition (e.g. "file not found"), we will
          //likely trigger the exception again - this time in our own thread. In
          //case of no underlying error condition, we will end up with the
          //correct result in this thread.
          if ( verbose )
            std::cout<< this->factoryName()
                     <<" (thread_"<<thread_details::currentThreadIDForPrint()<<")"
                     <<" : Restarting since other thread did not as expected create (from scratch) object for key "<<keystr<<std::endl;

          guard.ensureUnlock(NCRYSTAL_DEBUG_LOCKS_ARGS);
          return this->create(key);
        }
        guard.ensureUnlock(NCRYSTAL_DEBUG_LOCKS_ARGS);
      }
    }
  }
}

#endif

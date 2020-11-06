#ifndef NCrystal_FactoryUtils_hh
#define NCrystal_FactoryUtils_hh

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

#include "NCrystal/NCDefs.hh"
#include <chrono>
#include <thread>
#include <iostream>

namespace NCrystal {

  template< class TKey, class TValue, bool factoryKeepsOwnRef = false >
  class CachedFactoryBase {
  public:

    /////////////////////////////////////////////////////////////////////////////////
    //
    // Base class for multi-thread safe factories. Client code must implement the
    // actualCreate function, plus a few cosmetic functions. Factories will
    // normally keep just weak references (weak_ptr) to the created objects, but
    // can be configured to instead keep also strong references (shared_ptr) to
    // them (thus ensuring the created objects are always kept alive and won't have
    // to be recreated later). This is done by either setting the
    // factoryKeepsOwnRef template parameter to true, by calling the global
    // enableAllFactoriesKeepStrongRefs function, or by setting the
    // NCRYSTAL_FACTORIES_KEEPS_STRONG_REFS environment variable.
    //
    // The strong refs and any dangling weak refs can be cleaned up by an explicit
    // call the cleanup function of a particular factory, or simply by calling the
    // global clearCaches function which will in turn call the cleanup function of
    // all factories.
    //
    /////////////////////////////////////////////////////////////////////////////////

    typedef TKey key_type;
    typedef TValue value_type;
    typedef std::shared_ptr<const value_type> ShPtr;
    typedef std::weak_ptr<const value_type> WeakPtr;

    ShPtr create(const key_type&);

    virtual std::string keyToString( const key_type& ) const = 0;
    virtual const char* factoryName() const = 0;

    //Release all kept strong references and all dangling weak references (this
    //function automatically registered with and invoked by global clearCaches
    //function):
    void cleanup();

  protected:
    virtual ShPtr actualCreate(const key_type&) = 0;

  private:
    struct CacheEntry { bool underConstruction = false; WeakPtr weakPtr; };
    std::map<key_type,CacheEntry> m_cache;
    std::mutex m_mutex;
    std::vector<ShPtr> m_strongRefs;
    bool m_cleanupNeedsRegistry = true;
  };

  ///////////////////////////////////////////////////////////////////////////
  // Global control affecting all CachedFactoryBase instances (the default //
  // settings of these can also be controlled by environment variables).   //
  ///////////////////////////////////////////////////////////////////////////

  void enableFactoryVerbosity( bool status = true );
  bool getFactoryVerbosity();

  void enableAllFactoriesKeepStrongRefs( bool status = true );
  bool getAllFactoriesKeepStrongRefs();

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCrystal {

  template<class TKey,class TValue,bool factoryKeepsOwnRef>
  inline void CachedFactoryBase<TKey,TValue,factoryKeepsOwnRef>::cleanup()
  {
    std::lock_guard<std::mutex> guard(m_mutex);
    m_strongRefs.clear();
    auto it = m_cache.begin();
    auto itE = m_cache.end();
    while (it!=itE) {
      auto itNext = std::next(it);
      if ( ! it->second.weakPtr.lock() && ! it->second.underConstruction )
        m_cache.erase(it);
      it = itNext;
    }
  }

  template<class TKey,class TValue,bool factoryKeepsOwnRef>
  inline std::shared_ptr<const TValue> CachedFactoryBase<TKey,TValue,factoryKeepsOwnRef>::create(const TKey& key)
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    class Guard {
      //Local guard class. Kind of like std::lock_guard, but can remove set
      //underConstruction flag as well.
      std::mutex& m_mutex;
      bool* m_constructFlag = nullptr;
      bool m_isLocked = false;
    public:
      bool weHoldConstructFlag() const { return m_constructFlag != nullptr; }
      Guard( std::mutex& mutex ) : m_mutex(mutex) {}
      void setConstructFlagWithGuard(bool* constructFlag) {
        nc_assert(constructFlag);
        ensureLock();
        m_constructFlag = constructFlag;
        *m_constructFlag = true;
      }
      ~Guard() {
        if (m_constructFlag) {
          ensureLock();
          *m_constructFlag = false;
        }
        ensureUnlock();
      }
      void ensureLock()
      {
        if (!m_isLocked) {
          m_mutex.lock();
          m_isLocked = true;
        }
      }
      void ensureUnlock()
      {
        if (m_isLocked) {
          m_mutex.unlock();
          m_isLocked = false;
        }
      }
    };

    //////////////////////////////////////////////////////////////////////////////////////
    const bool verbose = getFactoryVerbosity();
    std::string keystr = ( verbose ? keyToString(key) : std::string() );

    Guard guard(m_mutex);
    guard.ensureLock();

    if ( m_cleanupNeedsRegistry ) {
      m_cleanupNeedsRegistry = false;
      std::function<void()> fct_cleanup = [this](){ this->cleanup(); };
      registerCacheCleanupFunction(fct_cleanup);
    }

    if ( verbose )
      std::cout<< this->factoryName()
               <<" (thread_"<<std::this_thread::get_id()<<")"
               <<" : Request to provide object for key "<<keystr<<std::endl;

    auto& cache_entry = m_cache[key];
    ShPtr res = cache_entry.weakPtr.lock();
    if (!!res) {
      if ( verbose )
        std::cout<< this->factoryName()
                 <<" (thread_"<<std::this_thread::get_id()<<")"
                 <<" : Return pre-existing cached object for key "<<keystr<<std::endl;
      nc_assert_always(!cache_entry.underConstruction);
      return res;//easy: already there
    }
    //Not there: check if already under construction or if we should construct:
    if (!cache_entry.underConstruction)
      guard.setConstructFlagWithGuard(&cache_entry.underConstruction);

    //Have to either construct ourselves or wait for other thread to do it. In
    //both cases, stop holding mutex lock:
    guard.ensureUnlock();

    if (guard.weHoldConstructFlag()) {
      if ( verbose )
        std::cout << this->factoryName()
                  <<" (thread_"<<std::this_thread::get_id()<<")"
                  << " : Creating (from scratch) object for key " << keystr << std::endl;
      //Invoke actual creation function without holding the mutex lock.
      res = actualCreate(key);
      //Populate result while holding mutex lock:
      guard.ensureLock();
      cache_entry = m_cache[key];//reacquire after getting lock back
      nc_assert_always(!cache_entry.weakPtr.lock());//no one else should have tried to create this
      cache_entry.weakPtr = res;
      if ( factoryKeepsOwnRef || getAllFactoriesKeepStrongRefs() )
        m_strongRefs.push_back(res);
      return res;
    } else {
      //Wait for other thread to populate cache. Sleep and recheck periodically.
      while (true) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        guard.ensureLock();
        cache_entry = m_cache[key];//reacquire after getting lock back
        if ( verbose )
          std::cout<< this->factoryName()
                   <<" (thread_"<<std::this_thread::get_id()<<")"
                   <<" : Waiting for other thread to create (from scratch) object for key "<<keystr<<std::endl;
        res = cache_entry.weakPtr.lock();
        if (!!res)
          return res;//success!
        //Not yet. Double-check other thread is still trying:
        if (!cache_entry.underConstruction) {
          //Not there and no other thread is currently trying to construct
          //it. Technically we can't know if this situation happened because the
          //other thread ended prematurely (e.g. exception thrown), or because
          //it succeeded, but the result was already used and cleaned up again
          //due to all shared pointers going out of scope (which can happen even
          //if factoryKeepsOwnRef=true, since cleanup() might have been
          //invoked). The best we can do is to try again, and hope this time
          //there won't be any concurrent attempts to create. In case of an
          //underlying error condition (e.g. "file not found"), we will likely
          //trigger the exception again - this time in our own thread. In case
          //of no underlying error condition, we will end up with the correct
          //result in this thread.
          if ( verbose )
            std::cout<< this->factoryName()
                     <<" (thread_"<<std::this_thread::get_id()<<")"
                     <<" : Restarting since other thread did not as expected create (from scratch) object for key "<<keystr<<std::endl;

          guard.ensureUnlock();
          return this->create(key);
        }
        guard.ensureUnlock();
      }
    }
  }
}


#endif

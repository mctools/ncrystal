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

#include "NCrystal/internal/NCFactoryUtils.hh"
#include "NCrystal/internal/NCString.hh"

namespace NC = NCrystal;

namespace NCrystal {
  namespace {
    static std::atomic<bool> s_factoryVerbosity( ncgetenv_bool("DEBUG_FACTORY")
                                                 || ncgetenv_bool("DEBUGFACTORY")
                                                 || ncgetenv_bool("DEBUG_FACT")
                                                 || ncgetenv_bool("DEBUGFACT") );
  }
}

void NC::enableFactoryVerbosity( bool status )
{
  s_factoryVerbosity = status;
}

bool NC::getFactoryVerbosity()
{
  return s_factoryVerbosity;
}

#ifndef NCRYSTAL_DISABLE_THREADS
namespace NCrystal {
  namespace detail {
    struct ThreadDeadLockDetectDB {
      std::mutex mtx;
      struct ThreadStatus {
        std::thread::id thrid;
        unsigned nWorking = 0;
        bool isWaiting = false;
        bool unused() const noexcept { return !isWaiting&& nWorking == 0; }
        bool operator<(const ThreadStatus& o ) const { return thrid < o.thrid; }
        bool operator<(const std::thread::id& o ) const { return thrid < o; }
        ThreadStatus( std::thread::id tid ) : thrid(tid) {}
      };
      SmallVector<ThreadStatus,64> threadStates;
      ThreadStatus& getThreadStatus( const std::thread::id thrid )
      {
        auto it = std::lower_bound(threadStates.begin(),threadStates.end(),thrid);
        if ( it != threadStates.end() && thrid == it->thrid )
          return *it;
        //ThreadStatus for thread in question missing, add and sort:
        threadStates.emplace_back(thrid);
        std::sort(threadStates.begin(),threadStates.end());
        return getThreadStatus(thrid);
      }

      unsigned nPossibleUnusedItems = 0;
      void cleanupUnused()
      {
        nc_assert(nPossibleUnusedItems>4);
        nPossibleUnusedItems = 0;
        auto n_orig = threadStates.size();
        //Sort which puts all unused entries at the end, followed by truncation:
        std::sort(threadStates.begin(),threadStates.end(),
                  [](const ThreadStatus & a, const ThreadStatus & b) -> bool
                  {
                    return a.unused() == b.unused() ? a < b : b.unused();
                  });
        while ( !threadStates.empty() && threadStates.back().unused() )
          threadStates.pop_back();
        nc_assert(std::is_sorted(threadStates.begin(),threadStates.end()));
        if ( s_factoryVerbosity.load() )
          std::cout<< "FactoryUtils dead-lock protection cleanup triggered discarding "
                   <<n_orig-threadStates.size()<<" unused thread state entries ("
                   <<threadStates.size()<<" remains)"<<std::endl;
      }

    };
    ThreadDeadLockDetectDB& getDeadLockDB() { static ThreadDeadLockDetectDB db; return db; }
    void registerThreadWork( std::thread::id thrid)
    {
      auto& db = getDeadLockDB();
      NCRYSTAL_LOCK_GUARD(db.mtx);
      db.getThreadStatus( thrid ).nWorking += 1;
    }
    void registerThreadWorkDone( std::thread::id thrid )
    {
      auto& db = getDeadLockDB();
      NCRYSTAL_LOCK_GUARD(db.mtx);
      auto& ts = db.getThreadStatus( thrid );
      nc_assert(ts.nWorking>=1);
      if ( ts.nWorking == 1 ) {
        ts.nWorking = 0;
        ++db.nPossibleUnusedItems;
        if ( db.nPossibleUnusedItems > 512 && db.nPossibleUnusedItems % 64 == 0 )
          db.cleanupUnused();
      } else {
        ts.nWorking -= 1;
      }
    }
    void registerThreadAsWaiting( std::thread::id thrid )
    {
      auto& db = getDeadLockDB();
      NCRYSTAL_LOCK_GUARD(db.mtx);

      auto& ts = db.getThreadStatus( thrid );
      nc_assert(!ts.isWaiting);
      ts.isWaiting = true;
      //Check that at least 1 thread is still working without waiting, otherwise
      //throw exception about bad cyclic setup (this would otherwise be a deadlock):
      for ( const auto& e : db.threadStates )
        if ( !e.isWaiting && e.nWorking > 0 )
          return;

      //No threads are working, all are waiting!
      NCRYSTAL_THROW(BadInput,"Cyclic dependency in factory request detected "
                     "(check your input configurations and data for cyclic references)!");
    }
    void registerThreadAsFinishedWaiting( std::thread::id thrid )
    {
      std::cout<<"registerThreadAsFinished :"<<thrid<<std::endl;
      auto& db = getDeadLockDB();
      NCRYSTAL_LOCK_GUARD(db.mtx);
      auto& ts = db.getThreadStatus( thrid );
      nc_assert(ts.isWaiting);
      ts.isWaiting = false;
    }
  }
}
#endif

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

#include "NCrystal/NCMem.hh"
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cassert>
#include <vector>
#include <mutex>

namespace NCrystal {
  static long s_RCBase_nInstances = 0;
  static int s_RCBase_dbgmem = -1;
  static int RCBase_dbgmem()
  {
    if (s_RCBase_dbgmem == -1) {
      s_RCBase_dbgmem = 0;
      const char * envc = std::getenv("NCRYSTAL_DEBUGMEM");
      if (envc) {
        std::string env(envc);
        if (env=="1") s_RCBase_dbgmem = 1;
        else if (env=="2") s_RCBase_dbgmem = 2;
        else if (env!="0")
          NCRYSTAL_THROW(BadInput,"NCRYSTAL_DEBUGMEM environment variable must"
                         " equal either \"0\", \"1\" or \"2\" if set.");
      }
    }
    return s_RCBase_dbgmem;
  }
  void RCBase::enableMemDbg(int lvl)
  {
    if (lvl<0||lvl>2)
      NCRYSTAL_THROW(BadInput,"Debug lvl should be 0, 1 or 2.");
    s_RCBase_dbgmem = lvl;
  }
}

NCrystal::RCBase::RCBase() throw()
  : m_refCount(0)
{
  ++s_RCBase_nInstances;
  switch(RCBase_dbgmem()) {
  case 1:
    printf( "NCrystal::RCBase(). Number of active RCBase instances is now %li\n",
            s_RCBase_nInstances );
    break;
  case 2:
    printf( "NCrystal::RCBase() [%p]. Number of active RCBase instances is now %li\n",
            (void*)this, s_RCBase_nInstances );
    break;
  default:
    break;
  }
}

NCrystal::RCBase::~RCBase()
{
  //use assert rather than nc_assert in destructors
  assert(m_refCount==0&&"Inconsistent call to destructor of ref-counted class");
  assert(s_RCBase_nInstances>0&&"Inconsistent triggering of RCBase destructor");
  --s_RCBase_nInstances;
  switch(RCBase_dbgmem()) {
  case 1:
    printf( "NCrystal::~RCBase(). Number of active RCBase instances is now %li\n",
            s_RCBase_nInstances );
    break;
  case 2:
    printf( "NCrystal::~RCBase() [%p]. Number of active RCBase instances is now %li\n",
            (void*)this, s_RCBase_nInstances );
    break;
  default:
    break;
  }
}

long NCrystal::RCBase::nInstances()
{
  return s_RCBase_nInstances;
}

namespace NCrystal {
  static std::mutex s_cacheCleanerMutex;
  static std::vector<std::function<void()>> s_cacheCleanerMutexFcts;
}

void NCrystal::clearCaches()
{
  std::lock_guard<std::mutex> lock(s_cacheCleanerMutex);
  for (auto& f : s_cacheCleanerMutexFcts)
    f();
}


void NCrystal::registerCacheCleanupFunction( std::function<void()> f )
{
  std::lock_guard<std::mutex> lock(s_cacheCleanerMutex);
  s_cacheCleanerMutexFcts.emplace_back(f);
}

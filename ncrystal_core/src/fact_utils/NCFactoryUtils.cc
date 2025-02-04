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

#include "NCrystal/internal/fact_utils/NCFactoryUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
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
#  include <thread>
#  include <sstream>
std::string NC::thread_details::currentThreadIDForPrint()
{
  std::ostringstream ss;
  ss << std::this_thread::get_id();
  return ss.str();
}
#else
std::string NC::thread_details::currentThreadIDForPrint()
{
  return "<thread-id-unavailable>";
}
#endif

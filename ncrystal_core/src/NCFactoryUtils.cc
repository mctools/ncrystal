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

#include "NCrystal/internal/NCFactoryUtils.hh"
namespace NC = NCrystal;

namespace NCrystal {
  namespace {
    static std::atomic<bool> s_factoryVerbosity(bool(getenv("NCRYSTAL_DEBUG_FACTORY")));
    static std::atomic<bool> s_factoriesKeepStrongRefs(bool(getenv("NCRYSTAL_FACTORIES_KEEPS_STRONG_REFS")));
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

void NC::enableAllFactoriesKeepStrongRefs( bool status )
{
  s_factoriesKeepStrongRefs = status;
}

bool NC::getAllFactoriesKeepStrongRefs()
{
  return s_factoriesKeepStrongRefs;
}

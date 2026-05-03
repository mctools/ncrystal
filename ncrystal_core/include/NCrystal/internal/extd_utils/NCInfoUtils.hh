#ifndef NCrystal_InfoUtils_hh
#define NCrystal_InfoUtils_hh

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

// A few utilities related to Info objects.

#include "NCrystal/interfaces/NCInfo.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace InfoUtils {

    //Try to dig out a particular DynInfo fra an Info object. If a target
    //atomDisplayLabel is not provided, this will fail in case of multiple
    //DynInfos.
    const DynamicInfo* findDynInfo( const Info& info,
                                    Optional<std::string>
                                    atomDisplayLabel = NullOpt );
  }
}

#endif

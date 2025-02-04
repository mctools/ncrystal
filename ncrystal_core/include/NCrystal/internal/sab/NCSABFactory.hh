#ifndef NCrystal_SABFactory_hh
#define NCrystal_SABFactory_hh

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

#include "NCrystal/interfaces/NCSABData.hh"
#include "NCrystal/internal/sab/NCSABScatterHelper.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SAB {

    //Direct factory function with no caching:
    std::unique_ptr<const SABScatterHelper> createScatterHelper( shared_obj<const SABData>,
                                                                 std::shared_ptr<const VectD> energyGrid = nullptr );

    //Same with caching:
    void clearScatterHelperCache();
    shared_obj<const SABScatterHelper> createScatterHelperWithCache( shared_obj<const SABData>,
                                                                     std::shared_ptr<const VectD> energyGrid = nullptr );

    //For caching reasons, we keep a database of energy grid's and an associated
    //unique id. Note that it is expected that most energy grids specified will
    //either be "unspecified" (nullptr or empty) or just 3 entries long (emin
    //emax npts). Thus, this cache is not expected to actually become very big:
    UniqueIDValue egridToUniqueID(const VectD& egrid);
    UniqueIDValue egridToUniqueID(const std::shared_ptr<const VectD>& egrid);
    std::shared_ptr<const VectD> egridFromUniqueID(UniqueIDValue);
  }

}

#endif

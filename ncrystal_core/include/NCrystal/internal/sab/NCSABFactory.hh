#ifndef NCrystal_SABFactory_hh
#define NCrystal_SABFactory_hh

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

#include "NCrystal/interfaces/NCSABData.hh"
#include "NCrystal/internal/sab/NCSABScatterHelper.hh"
#include "NCrystal/internal/sab/NCSABExtended.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SAB {

    ////////////////////////
    // Create SABExtended //
    ////////////////////////

    shared_obj<const SABUtils::SABExtended>
    createSABExtendedNoCache( int knllux,
                              shared_obj<const SABData>,
                              std::shared_ptr<const VectD>
                              energyGrid = nullptr );

    shared_obj<const SABUtils::SABExtended>
    createSABExtendedWithCache( int knllux,
                                shared_obj<const SABData>,
                                std::shared_ptr<const VectD>
                                energyGrid = nullptr );

    ///////////////////////////////
    // Create shared energy grid //
    ///////////////////////////////

    //For caching reasons, we keep a database of energy grid's and an associated
    //unique id. Note that it is expected that most energy grids specified will
    //either be "unspecified" (nullptr or empty) or just 3 entries long (emin
    //emax npts). Thus, this cache is not expected to actually become very big:
    UniqueIDValue egridToUniqueID(const VectD& egrid);
    UniqueIDValue egridToUniqueID(const std::shared_ptr<const VectD>& egrid);
    std::shared_ptr<const VectD> egridFromUniqueID(UniqueIDValue);

    ///////////////////////
    // Legacy algorithms //
    ///////////////////////

    //For reference and validation work, it is possible to disable the
    //beta-endpoint fix (introduced in NCrystal 3.1.0), and we can reduce
    //sampling artifacts by sampling at higher energies (increasing
    //SABSampler::EGridMargin). The numerical value of the enum option is the
    //knllux value needed to select the legacy algorithm:

    enum class LegacySABAlgOpts : int {
      DEFAULT = -1,
      OVERSAMPLE10 = -2,
      OVERSAMPLE50 = -3,
      NOBETAFIX = -4,
      NOBETAFIX_OVERSAMPLE10 = -5,
      NOBETAFIX_OVERSAMPLE50 = -6,
      //For reference, provide access to the range of values covered:
      MIN = -6,
      MAX = -1
    };

    //Direct factory function with no caching:
    std::unique_ptr<const SABScatterHelper>
    createScatterHelper( shared_obj<const SABData>,
                         std::shared_ptr<const VectD> energyGrid = nullptr,
                         LegacySABAlgOpts = LegacySABAlgOpts::DEFAULT );

    //Same with caching:
    shared_obj<const SABScatterHelper>
    createScatterHelperWithCache( shared_obj<const SABData>,
                                  std::shared_ptr<const VectD> egrid = nullptr,
                                  LegacySABAlgOpts = LegacySABAlgOpts::DEFAULT );

  }

}

#endif

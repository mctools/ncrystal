#ifndef NCrystal_MMC_Query_hh
#define NCrystal_MMC_Query_hh

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

#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/minimc/NCMMC_CBMgr.hh"

namespace NCRYSTAL_NAMESPACE {

  class MatCfg;

  namespace MiniMC {

    namespace Query {

      //Generic JSON query for the MiniMC subsystem (format to be described
      //elsewhere, FIXME):
      using Query = SmallVector<StrView,8>;
      void JSONQuery( std::ostream&, const Query& );

      //Same, but allowing simulation callback. If a callback is provided, query
      //*must* be an ["mmc","run",...] query.
      void JSONQuery_flexmmcrun( std::ostream&, const Query&,
                                 const Optional<CB::CBMgrInput>& );

    }
  }
}

#endif

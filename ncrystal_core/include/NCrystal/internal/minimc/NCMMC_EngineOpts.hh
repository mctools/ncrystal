#ifndef NCrystal_MMC_EngineOpts_hh
#define NCrystal_MMC_EngineOpts_hh

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/utils/NCStrView.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    struct EngineOpts {
      //FIXME: Some docs here (or in nctool output / wiki).

      enum class IgnoreMiss : uint32_t { NO = 0,
                                         YES = 1,
                                         Default = NO };
      enum class TallyBreakdown : uint32_t { NO = 0,
                                             YES = 1,
                                             Default = YES };
      enum class IncludeAbsorption : uint32_t { NO = 0,
                                                YES = 1,
                                                Default = YES };
      IgnoreMiss ignoreMiss = IgnoreMiss::Default;
      ThreadCount nthreads = ThreadCount::auto_detect();
      TallyBreakdown tallyBreakdown = TallyBreakdown::Default;
      IncludeAbsorption includeAbsorption = IncludeAbsorption::Default;
    };

    //Parse to/from a string representation (can also be used to normalise a
    //string representation). The string representation will be compact, and not
    //include entries at default values:
    EngineOpts parseEngineOpts( StrView );
    std::string engineOptsToString( const EngineOpts& eopts );

    //Default streaming adds "EngineOpts(...)" around the string:
    std::ostream& operator<<( std::ostream& os, const EngineOpts& eopts );

    //Output as a JSON dictionary. This will always include all options,
    //including those at default values.
    void engineOptsToJSON(std::ostream&, const EngineOpts&);

  }
}

#endif

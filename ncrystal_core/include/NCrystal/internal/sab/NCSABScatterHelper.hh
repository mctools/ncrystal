#ifndef NCrystal_SABScatterHelper_hh
#define NCrystal_SABScatterHelper_hh

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

#include "NCrystal/internal/sab/NCSABSampler.hh"
#include "NCrystal/internal/sab/NCSABXSProvider.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SAB {

    //Convenience struct and helper class for when both sampler and xsprovider
    //are both needed and with similar lifetimes. Should normally be constructed
    //by using either the SABFactory functions or a SABIntegrator.
    //
    //For reference it can optionally hold a JSON descriptive dictionary (in a
    //format suitable for ProcImpl::specificJSONDescription).

    class SABScatterHelper : private MoveOnly {
    public:
      SABScatterHelper( SABXSProvider&& xp,
                        SABSampler&&sp,
                        Optional<std::string> json = NullOpt )
        : xsprovider(std::move(xp)),
          sampler(std::move(sp)),
          specificJSONDescription(std::move(json))
      {
      }
      SABScatterHelper() = default;//incomplete
      SABScatterHelper( SABScatterHelper&& ) = default;
      SABScatterHelper& operator=( SABScatterHelper&& ) = default;
      SABXSProvider xsprovider;
      SABSampler sampler;
      Optional<std::string> specificJSONDescription;

    };

  }

}

#endif

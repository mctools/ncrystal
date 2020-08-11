#ifndef NCrystal_SABScatterHelper_hh
#define NCrystal_SABScatterHelper_hh

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

#include "NCSABSampler.hh"
#include "NCSABXSProvider.hh"

namespace NCrystal {

  namespace SAB {

    //Convenience struct and helper class for when both sampler and xsprovider
    //are both needed and with similar lifetimes. Should normally be constructed
    //by using either the SABFactory functions or a SABIntegrator.

    class SABScatterHelper : private MoveOnly {
    public:
      SABScatterHelper(SABXSProvider&& xp,SABSampler&&sp) : xsprovider(std::move(xp)), sampler(std::move(sp)) {}
      SABScatterHelper() = default;//incomplete
      SABScatterHelper( SABScatterHelper&& ) = default;
      SABScatterHelper& operator=( SABScatterHelper&& ) = default;
      SABXSProvider xsprovider;
      SABSampler sampler;
    };

  }

}

#endif

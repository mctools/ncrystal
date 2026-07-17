#ifndef NCrystal_SABRefSampler_hh
#define NCrystal_SABRefSampler_hh

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

namespace NCRYSTAL_NAMESPACE {

  namespace SABRef {

    //////////////////////////////////////////////////////
    //                                                  //
    // Utility for reference sampling an S(alpha,beta)  //
    // table at a given energy point.                   //
    //                                                  //
    // Its primary purpose is to be used as a reference //
    // in tests, and not to be used for actual end-user //
    // sampling.                                        //
    //                                                  //
    //////////////////////////////////////////////////////

    //Sample multiple (alpha,beta) values at once:
    std::pair<VectD,VectD> refSampleAlphaBeta( RNG&,
                                               const SABData&,
                                               double E_div_kT,
                                               std::uint64_t nsample );

  }
}

#endif

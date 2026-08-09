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
#include "NCrystal/internal/sab/NCSABExtender.hh"

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

    //Sample multiple (alpha,beta) values at once. At energies beyond what is
    //covered by the sabdata ("Emax"), it is optionally possible to sample with
    //a SABExtender by providing a RefSampleExtension object:
    struct RefSampleExtension {
      //Extender (default is a free-gas extender).
      //NB: MUST HAVE SigmaBound=1!! (fixme enshrine in types)
      std::shared_ptr<SAB::SABExtender> extender = nullptr;
      //Emax (0eV means autodetect, with suggestedEMax from the SABData or an
      //expensive SABProcessor initialisation).
      NeutronEnergy emax = NeutronEnergy{0.0};
    };

    std::pair<VectD,VectD>
    refSampleAlphaBeta( RNG&,
                        shared_obj<const SABData>,
                        double E_div_kT,
                        std::uint64_t nsample,
                        Optional<RefSampleExtension> = NullOpt);

  }
}

#endif

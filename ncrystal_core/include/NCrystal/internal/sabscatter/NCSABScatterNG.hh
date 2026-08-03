#ifndef NCrystal_SABScatterNG_hh
#define NCrystal_SABScatterNG_hh

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

#include "NCrystal/interfaces/NCProcImpl.hh"

//fixme: rename NG<->Legacy? Or put all new stuff in new namespace?

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {
    class SABExtended;
  }

  class SABScatterNG final : public ProcImpl::ScatterIsotropicMat {
  public:

    //Provides cross-sections and samplings based on an S(alpha,beta) scattering
    //kernel.

    //Technically, a SABScatterNG instance essentially just wraps a SABExtended
    //instance as a Process, and applies a scale to it. The SABExtended instance
    //provides cross sections for a hypothetical SigmaBound=1barn, so the scale
    //should include the actual SigmaBound. Additionally it can include factors
    //related to composition, etc.

    const char * name() const noexcept override { return "SABScatterNG"; }

    using SABExtended = SABUtils::SABExtended;

    SABScatterNG( shared_obj<const SABExtended>, double scale );

    //Strongly typed SigmaBound (actual scale will become SigmaBound*scale):
    SABScatterNG( shared_obj<const SABExtended>,
                  SigmaBound, double scale = 1.0 );

    virtual ~SABScatterNG();

    CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ) const override;

    ScatterOutcomeIsotropic
    sampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const override;

    std::shared_ptr<Process> createMerged( const Process& other,
                                           double scale_self,
                                           double scale_other ) const override;

  protected:
    shared_obj<const SABExtended> m_sh;
    double m_scale = 1.0;
    Optional<std::string> specificJSONDescription() const override;
  };

}

#endif

#ifndef NCrystal_BkgdExtCurve_hh
#define NCrystal_BkgdExtCurve_hh

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

#include "NCrystal/interfaces/NCProcImpl.hh"
#include "NCrystal/interfaces/NCInfo.hh"

namespace NCRYSTAL_NAMESPACE {

  class BkgdExtCurve final : public ProcImpl::ScatterIsotropicMat {
  public:

    //Calculates background (non-Bragg) scattering in a crystal, based on
    //external functions for calculating cross-sections, using the XSectProvider
    //section in the passed Info object. Scatterings will be elastic and
    //isotropic.

    const char * name() const noexcept final { return "BkgdExtCurve"; }

    BkgdExtCurve( shared_obj<const Info> );
    virtual ~BkgdExtCurve();

    CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ) const final;
    ScatterOutcomeIsotropic sampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const final;
    ScatterOutcome sampleScatter(CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const final;

  protected:
    shared_obj<const Info> m_ci;
  };
}

#endif

#ifndef NCrystal_FreeGas_hh
#define NCrystal_FreeGas_hh

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
#include "NCrystal/interfaces/NCAtomData.hh"

namespace NCRYSTAL_NAMESPACE {

  class FreeGas final : public ProcImpl::ScatterIsotropicMat {
  public:

    ////////////////////////////////////////////////////////////////////////////////
    //                                                                            //
    // Thermal scattering with a free-gas model. Specifically this carefully      //
    // implements cross-sections and sampling of scatterings based on the neutron //
    // scattering function, S, given in eq. 19 of the NCrystal sampling paper     //
    // (10.1016/j.jcp.2018.11.043) and is thus a fully quantum mechanical model.  //
    //                                                                            //
    // See NCFreeGasUtils.hh for more details.                                    //
    //                                                                            //
    ////////////////////////////////////////////////////////////////////////////////

    const char * name() const noexcept final { return "FreeGas"; }

    //Explicitly provide target parameters or take parameters from AtomData object:
    FreeGas( Temperature, AtomMass, SigmaFree );
    FreeGas( Temperature, AtomMass, SigmaBound );
    FreeGas( Temperature, const AtomData& );

    CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ) const override;
    ScatterOutcomeIsotropic sampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const override;

    virtual ~FreeGas();

  protected:
    Optional<std::string> specificJSONDescription() const override;
    struct Impl;
    Pimpl<Impl> m_impl;
  };
}

#endif

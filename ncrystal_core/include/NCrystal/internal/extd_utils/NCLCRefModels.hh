#ifndef NCrystal_LCRefModels_hh
#define NCrystal_LCRefModels_hh

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
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"

namespace NCRYSTAL_NAMESPACE {

  class LCBraggRef final : public ProcImpl::ScatterAnisotropicMat {
  public:
    //Simple but very slow implementation of layered crystals. Mainly provided
    //as a reference (should give increasingly better result with higher
    //nsample).
    LCBraggRef(ProcImpl::ProcPtr scbragg, LCAxis lcaxis_lab, unsigned nsample = 1000);

    const char * name() const noexcept override { return "LCBraggRef"; }
    CrossSect crossSection(CachePtr&, NeutronEnergy, const NeutronDirection& ) const override;
    EnergyDomain domain() const noexcept override;
    ScatterOutcome sampleScatter(CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const override;
  private:
    ProcImpl::ProcPtr m_sc;
    Vector m_lcaxislab;
    unsigned m_nsample;
    unsigned m_nsampleprime;
  };

  class LCBraggRndmRot final : public ProcImpl::ScatterAnisotropicMat {
  public:
    //Like LCBraggRef, but using random crystallite rotations even for
    //crossSection calls - and reusing the same orientations in a subsequent
    //call to generateScattering. Again, this is mainly provided as a reference
    //and not really recommended for general usage.
    //
    //WARNING: The crossSection function of this class use the global default
    //RNG stream, and is therefore NOT MT-safe.
    LCBraggRndmRot(ProcImpl::ProcPtr scbragg, LCAxis lcaxis_lab, unsigned nsample = 1);
    const char * name() const noexcept override { return "LCBraggRndmRot"; }
    CrossSect crossSection(CachePtr&, NeutronEnergy, const NeutronDirection& ) const override;
    EnergyDomain domain() const noexcept override;
    ScatterOutcome sampleScatter(CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const override;
  private:
    ProcImpl::ProcPtr m_sc;
    Vector m_lcaxislab;
    unsigned m_nsample;
    class Cache : public CacheBase {
    public:
      std::vector<PhiRot> rotations;//rotations sampled
      VectD xscommul;//cross-sections at the sampled rotations.
      CachePtr sc_cacheptr;//for passing to m_sc
      std::pair<NeutronEnergy,Vector> neutron_state = {NeutronEnergy{-1.0},Vector{0,0,0}};
      void invalidateCache() override { neutron_state = {NeutronEnergy{-1.0},Vector{0,0,0}}; }
    };
    void updateCache(Cache&,NeutronEnergy,const Vector&) const;
  };

}

#endif


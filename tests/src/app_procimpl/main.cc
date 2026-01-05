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
namespace NC = NCrystal;

namespace {

  class MyIsoScat final : public NC::ProcImpl::ScatterIsotropicMat {
  public:
    const char* name() const noexcept final { return "MyIsoScat"; }
    NC::CrossSect crossSectionIsotropic(NC::CachePtr&, NC::NeutronEnergy ) const final
    {
      return NC::CrossSect{1.0};
    }
    NC::ScatterOutcomeIsotropic sampleScatterIsotropic(NC::CachePtr&, NC::RNG&, NC::NeutronEnergy ) const final
    {
      return { NC::NeutronEnergy{0.0}, NC::CosineScatAngle{0.5} };
    }
  };

  class MyAnisoScat final : public NC::ProcImpl::ScatterAnisotropicMat {
  public:
    const char* name() const noexcept final { return "MyAnisoScat"; }

    NC::CrossSect crossSection(NC::CachePtr&, NC::NeutronEnergy, const NC::NeutronDirection& ) const final
    {
      return NC::CrossSect{1.0};
    }

    NC::ScatterOutcome sampleScatter(NC::CachePtr&, NC::RNG&, NC::NeutronEnergy, const NC::NeutronDirection& ) const final
    {
      return { NC::NeutronEnergy(0.0), NC::NeutronDirection{0,0,1} };
    }

  };

  class MyIsoAbs final : public NC::ProcImpl::AbsorptionIsotropicMat {
  public:
    const char* name() const noexcept final { return "MyIsoAbs"; }
    NC::CrossSect crossSectionIsotropic(NC::CachePtr&, NC::NeutronEnergy ) const final
    {
      return NC::CrossSect{1.0};
    }
  };

}

int main() {

  //The testing is really just compile-time at this point.
  (void)std::make_shared<MyIsoScat>();
  (void)std::make_shared<MyAnisoScat>();
  (void)std::make_shared<MyIsoAbs>();
  (void)std::make_shared<NC::ProcImpl::ProcComposition>();
  (void)std::make_shared<NC::ProcImpl::NullScatter>();
  (void)std::make_shared<NC::ProcImpl::NullAbsorption>();

}

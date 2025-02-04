#ifndef NCrystal_AbsOOV_hh
#define NCrystal_AbsOOV_hh

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

  class AbsOOV : public ProcImpl::AbsorptionIsotropicMat {
  public:

    // Provide absorption cross section based on simple 1/velocity
    // (=OneOverVelocity=OOV) scaling. This is non-oriented.

    AbsOOV( SigmaAbsorption );

    const char * name() const noexcept final { return "AbsOOV"; }
    CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ) const final;

    EnergyDomain domain() const noexcept override { return m_domain; }

    //Simple additive merge:
    std::shared_ptr<Process> createMerged( const Process& other,
                                           double scale_self,
                                           double scale_other ) const override;

    //For initialising directly from Info objects:
    AbsOOV( const Info& info ) : AbsOOV( info.getXSectAbsorption() ) {}

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
    void evalManyXSIsotropic( CachePtr&, const double* ekin, std::size_t N,
                              double* out_xs ) const override;
#endif

  protected:
    Optional<std::string> specificJSONDescription() const override;
  private:
    double m_c;
    EnergyDomain m_domain;
  };
}

#endif

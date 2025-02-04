#ifndef NCrystal_SANSSphScat_hh
#define NCrystal_SANSSphScat_hh

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
#include "NCrystal/internal/extd_utils/NCSANSUtils.hh"

namespace NCRYSTAL_NAMESPACE {

  class SANSSphereScatter final : public ProcImpl::ScatterIsotropicMat {
  public:

    // Implements the simplest hard-sphere SANS model, assuming no
    // particle-particle interactions (i.e. S(Q)=1) and:
    //
    // P(Q) = [ 3*(sin(QR)-QR*cos(QR))/(QR)^3 ]^2
    //
    // Where R is the sphere radius. This is used in [eq.1] from NCSANSUtils.hh:
    //
    //   sigma_sans(k) = C_sans * (Vp/k^2) * integral_0^2k{Q*P(Q)*S(Q)}dQ
    //
    // Cross sections are evaluated at or around machine precision for any input
    // data, relying on an analytical evaluation of the cross section integral,
    // which is itself evaluated directly or via a Taylor expansion as
    // appropriate for numerical stability. See the corresponding .cc file for
    // details about formulas etc.

    const char * name() const noexcept override { return "SANSSphereScatter"; }

    //Constructor takes the sphere radius (in Angstrom), as well as SLD contrast
    //("delta-rho") and atomic number density:
    struct sphere_radius { double value; };
    SANSSphereScatter( SANSScaleFactor C_sans, sphere_radius );

    CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ) const override;
    ScatterOutcomeIsotropic sampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const override;

    //Simple additive merge:
    std::shared_ptr<Process> createMerged( const Process& other,
                                           double scale_self,
                                           double scale_other ) const override;


#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
    bool isPureElasticScatter() const override { return true; }
#endif

    struct detail_direct_t;
    SANSSphereScatter(detail_direct_t,double r, double scale);
  protected:
    Optional<std::string> specificJSONDescription() const override;
  private:
    double m_r, m_scale;
  };

}

#endif

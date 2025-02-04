#ifndef NCrystal_SANSUtils_hh
#define NCrystal_SANSUtils_hh

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

// SANS related utilities.

#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  class Info;

  //////////////////////////////////////////////////////////////////////////////////
  //
  // Given nano particles (e.g. spheres) with per-particle volume of Vp, form
  // factor P(Q), structure factor S(Q), SLD contrast of drho, and volume
  // fraction of the nano particles of phi, and finally an atomic number density
  // of nd, then the SANS contribution to the per-atom cross section value is
  // given by:
  //
  //   sigma_sans(k) = (2pi*phi*drho^2/nd) * (Vp/k^2) * integral_0^2k{ Q*P(Q)*S(Q) }dQ  [eq.1]
  //
  // See notes further down this file for references of this formula.
  //
  // The first paranthesis, C_sans=2pi*(phi/nd)*drho^2 captures a common SANS
  // scale not related to the specifics of a particular model, which gives us
  // the slightly more convenient form:
  //
  //   sigma_sans(k) = C_sans * (Vp/k^2) * integral_0^2k{Q*P(Q)*S(Q)}dQ
  //
  // If Vp has units of Aa^3 and we want sigma_sans(k) in barn, it is thus most
  // convenient to express C_sans in units of barn/Aa^3=1/cm.
  //
  // For convenience, safety, and in particular to avoid having to make the same
  // considerations of unit conversions in all SANS models, we define the
  // following strong type for C_sans. The constructor taking (drho,nd,phi)
  // makes it trivial to instantiate without having to fret over units.
  //
  //////////////////////////////////////////////////////////////////////////////////

  class NCRYSTAL_API SANSScaleFactor final : public EncapsulatedValue<SANSScaleFactor> {
  public:
    // This scale factor gives the SANS cross section if multiplied with
    // (Vp/k^2)*integral_0^2k{Q*P(Q)*S(Q)}dQ, where Vp, P(Q), and S(Q) depends
    // on the SANS model in question (assuming xsect units of barn and the unit
    // of Vp to be Aa^3).
    //
    // Specifically, the scale factor is defined as 2pi*(phi/nd)*drho^2
    // evaluated in units of barn/Aa^3 = cm^-1.
    //
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return "barn/Aa^3"; }
    SANSScaleFactor( SLDContrast drho, NumberDensity nd, double phi );
    ncconstexpr17 void validate() const noexcept {}
  };

  //////////////////////////////////////////////////////////////////////////////////
  //
  // Utility functions for SANS plugins which will examine a multiphase Info
  // object for occurances of a particular custom section
  // @CUSTOM_<customsectionname> in the phase list. It is assumed that SANS
  // physics will be present with the SLD contrast occuring between each such
  // phase ("phase A") and the phase following it in the list ("phase B"). It is
  // allowed for more than one such SANS setup (either directly in the list or
  // in daughter lists), which is why a vector of results is returned. Each
  // result will contain the relevant SANSScaleFactor and the actual content of
  // the associated custom data section.
  //
  // If any phase is itself multiphased, the search for the custom section in
  // question will be carried out recursively.
  //
  // In all cases, it is an error if any "B" phase contains the target custom
  // section.
  //
  //////////////////////////////////////////////////////////////////////////////////

  struct NCRYSTAL_API CustomSansPluginData {
    SANSScaleFactor scale;
    std::vector<VectS> customData;
  };
  NCRYSTAL_API std::vector<CustomSansPluginData> extractCustomDataForSANSPlugin( const Info& info,
                                                                                 const std::string& customsectionname );

  //To simply check whether or not the previous method would provide anything, it is cheaper to use:
  NCRYSTAL_API bool hasCustomDataForSANSPlugin( const Info& info, const std::string& customsectionname );
}

////////////////////////////////////////////////////////////////////////////////////////
//
// A few more details of how we got [eq.1] above in the form presented:
//
// Generally speaking, the macroscopic SANS cross section from scattering on
// certain nano particles can be described with (for instance see eq. 1 in "The
// SANS Toolbox", Hammouda 2012 available at
// https://www.ncnr.nist.gov/staff/hammouda/the_SANS_toolbox.pdf):
//
// MacroXS(Q) = (Np/V)*Vp^2*drho^2*P(Q)*S(Q)
//
// This has unit of inverse length and can be either interpreted as inverse mean
// free path or cross-section per volume. Note that N/V in Hammouda's eq. 1 is
// the number density of the nano particles, not of individual atoms, so we
// write Np/V here.
//
// To get xsext-per-atom we must divide by number density nd=N/V where N is the
// number of atoms:
//
// xs-per-atom(Q) = (V/N)*(Np/V)*Vp^2*drho^2*P(Q)*S(Q)
//
// Now, the total number of particles (e.g. spheres) must be the total
// volume occupied by the particles divided by the volume of a single particle, so Np =
// (phi*V)/Vp => Np/V = phi/Vp. Thus:
//
// xs-per-atom(Q) = (V/N)*(phi/Vp)*Vp^2*drho^2*P(Q)*S(Q)
//                = (V/N)*phi*Vp*drho^2*P(Q)*S(Q)
//                = (phi/nd)*Vp*drho^2*P(Q)*S(Q)
//
// Finally, this depends on Q, to get total xs we must integrate the above
// expression wrt. Q.
//
// Treating SANS as elastic scattering, integration over all solid angles
// amounts to (2pi/k^2)*integral_0^2k(...)*Q*dQ. Thus we arrive at the target formula:
//
// tot-xs-per-atom = (2pi*phi*drho^2/nd) * (Vp/k^2) * integral_0^2k{ Q*P(Q)*S(Q) }dQ
//                   C_sans * (Vp/k^2) * integral_0^2k{ Q*P(Q)*S(Q) }dQ
//
//////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::SANSScaleFactor::SANSScaleFactor( SLDContrast drho, NumberDensity nd, double phi )
  : SANSScaleFactor( [drho,nd,phi]()
  {
    drho.validate();
    nd.validate();
    if ( !(phi>0.0) || !(phi<0.5) )
      NCRYSTAL_THROW2(BadInput,"SANS phi value out of range: "<< phi <<" (should be >0.0 and <0.5)");
    //Units of drho^2/nd are (1e-6/Aa^2)^2 / (1/Aa^3) = 1e-12/Aa, so we must
    //multiply with 1e-4 to get barn/Aa^3=1e-8/Aa. And of course, we also must include 2pi*phi:
    constexpr double kkk = k2Pi * 1e-4;
    return kkk * phi * drho.dbl() * drho.dbl() / nd.dbl();
  }())
{
}


#endif

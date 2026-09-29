#ifndef NCrystal_SCTUtils_hh
#define NCrystal_SCTUtils_hh

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

#include "NCrystal/internal/phys_utils/NCFreeGasUtils.hh"
#include "NCrystal/internal/utils/NCSpline.hh"

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Short-collision-time (SCT) model: cross sections (SCTXSProvider) and       //
// (alpha,beta) sampling (SCTSampler) for neutrons scattering                 //
// inelastically on bound atoms at high energy/momentum transfer. It is       //
// the companion of the free-gas model in NCFreeGasUtils.hh, to which it      //
// reduces exactly when Teff==T.                                              //
//                                                                            //
// PHYSICS: At high momentum transfer a scattering event is over so           //
// quickly (a "short collision time") that the struck atom reacts as if       //
// free, moving with whatever velocity it had at that instant -- the         //
// impulse approximation, going back to Lamb [1]. A bound atom's velocity     //
// spread is however not the Maxwell spread at the material temperature       //
// T: each vibrational mode of energy hw carries a mean energy of            //
// (hw/2)*coth(hw/2kT), including its zero point motion, half of it           //
// kinetic. Weighting with the phonon density of states rho gives the         //
// "effective temperature",                                                   //
//                                                                            //
//    k*Teff = integral( rho(hw) * (hw/2) * coth(hw/2kT) )  >= k*T            //
//                                                                            //
// so that the mean atomic kinetic energy is (3/2)k*Teff (this is the         //
// same quantity computed by VDOSEval::calcEffectiveTemperature and           //
// tabulated in ENDF File 7 evaluations [3]). Thus at high energy and         //
// momentum transfer a bound atom scatters like a free gas at Teff, with      //
// the Doppler width of the recoil ridge set by Teff rather than T. A         //
// literal "free gas at Teff" model would however satisfy detailed            //
// balance at Teff instead of at the true temperature T, giving a grossly     //
// overestimated upscattering tail (see e.g. [6]). The standard fix,          //
// introduced for the ENDF thermal scattering files by the Gulf General       //
// Atomic ("GA") team [2] and in use ever since in e.g. NJOY's THERMR         //
// module [4,5], is to keep the free-gas-at-Teff form on the downscatter      //
// side only (where the impulse approximation is the correct physics),        //
// and to *define* the upscatter side by detailed balance at T. In the        //
// symmetric-S ENDF convention this reads (with alpha_E and beta the          //
// usual ENDF unit-less momentum/energy transfers, cf. below):                //
//                                                                            //
//   S_SCT(alpha_E,-|beta|)                                                   //
//     = exp( -(alpha_E-|beta|)^2*T/(4*alpha_E*Teff) - |beta|/2 )             //
//       / sqrt( 4*pi*alpha_E*Teff/T )                                        //
//                                                                            //
// with S_SCT(alpha_E,+|beta|) then following from symmetry (i.e. from        //
// detailed balance at T). NB: several published documents contain            //
// misprinted versions of this formula; reference [6] catalogues them and     //
// states the intended ("GA") form, which is also what the above is.          //
// Setting Teff=T recovers exactly the free-gas law, so the free-gas          //
// extension model is precisely the special case Teff=T of this model.        //
//                                                                            //
// CONVENTIONS: beta = (Efinal-Eini)/kT (beta<0 means neutron energy          //
// loss, "downscatter") and NCrystal's alpha convention is                    //
// alpha = ( Efinal+Eini-2*mu*sqrt(Eini*Efinal) ) / kT, i.e. WITHOUT the      //
// division by the target-to-neutron mass ratio A=M/m_n which the ENDF        //
// convention employs: alpha_E = alpha/A = hbar^2*q^2/(2*M*kT) (the           //
// ENDF/van-Hove convention, in which the free-gas and SCT formulas take      //
// their clean textbook forms, with the recoil ridge at beta=-alpha_E).       //
// We use the asymmetric S throughout (like the rest of NCrystal, cf.         //
// [8]), related to the symmetric one by S_asym=exp(-beta/2)*S_sym and        //
// satisfying detailed balance as S(alpha,beta)=exp(-beta)S(alpha,-beta).     //
// With c = Ekin/kT, cross sections follow from S via                         //
// sigma = sigma_bound*A/(4c) * integral( S dalpha_E dbeta ) over the         //
// kinematically accessible region (equivalently sigma_bound/(4c) times       //
// the integral in NCrystal's alpha units).                                   //
//                                                                            //
// IMPLEMENTATION (all formulas below were derived independently from        //
// the above model definition; the one nontrivial antiderivative is           //
// elementary and easily verified by differentiation):                        //
//                                                                            //
// Sampling: denote by (alpha*,beta*) the transfers of a given scattering     //
// event expressed in units of k*Teff instead of k*T, i.e. with r =           //
// Teff/T: (alpha,beta) = (r*alpha*, r*beta*). Substitution shows that on     //
// the downscatter side the SCT law is *identical*, event for event, to      //
// the free-gas law at Teff expressed in its own units, and that the          //
// kinematic domain maps onto itself under this scaling. The upscatter        //
// side differs only by the alpha-independent detailed-balance factor         //
// exp(-beta*(1-1/r)) = exp(-beta**(r-1)) <= 1. SCTSampler therefore          //
// simply runs a FreeGasSampler at Teff, multiplies the result by r, and      //
// accepts upscatter candidates with the above probability (rejection         //
// sampling; the acceptance rate is >=sigma_SCT/sigma_FG(Teff), which is      //
// never small at the energies where the model is used).                      //
//                                                                            //
// Cross section: by the same decomposition,                                  //
//                                                                            //
//   sigma_SCT(E) = sigma_FG(Teff)(E)                                         //
//     - (sigma_b*A/4c) * integral_{beta>0} I(beta)*(1-exp(-(r-1)beta))       //
//                                                                            //
// with c = E/kTeff and all transfers in units of kTeff, where I(beta) is     //
// the free-gas "beta-kernel", i.e. the asymmetric free-gas S integrated      //
// over the kinematically allowed alpha range. Substituting                   //
// u=sqrt(alpha_E), I has the closed form                                     //
//                                                                            //
//   I(beta) = 0.5*[ erf(y(u)) + exp(-beta)*erf(z(u)) ]_{u-}^{u+},            //
//   y(u) = (u^2+beta)/(2u),  z(u) = (u^2-beta)/(2u),                         //
//                                                                            //
// following from the elementary antiderivative                               //
//   int exp(-(u/2+beta/(2u))^2) du                                           //
//     = (sqrt(pi)/2)*[ erf(u/2+beta/(2u))                                    //
//                      + exp(-beta)*erf(u/2-beta/(2u)) ] + C.                //
//                                                                            //
// For min(c*A,c/A) large the erf's saturate, with corrections of order       //
// erfc(sqrt(c*A))+erfc(sqrt(c/A)), leaving exactly I(beta)=exp(-beta)        //
// (the familiar flat-in-E' downscatter spectrum plus exp(-beta)              //
// upscatter tail of a hot free gas), whence the correction integral          //
// becomes exactly (r-1)/r. Since c = E/kTeff >= O(50) whenever the model     //
// is used above a typical S(alpha,beta) table's Emax of several eV, this     //
// closed form covers most practical evaluations to double precision.         //
// The transition region at lower c (reachable for e.g. data with a low       //
// hardwired egrid Emax, or heavy targets where c/A saturates slowly) is      //
// covered by a construction-time cubic-spline lookup table of the           //
// correction integral versus ln(c) (node values via fixed-order Romberg;     //
// the spline is linear in the node data, which keeps smooth parameter        //
// scans smooth -- unlike shape-limited interpolations such as PCHIP),        //
// with direct adaptive quadrature of the closed form kernel as the           //
// ultimate fallback below the table range.                                   //
//                                                                            //
// The total cross section is only mildly affected by the difference          //
// between SCT and a free gas at T (an O(kTeff/E) effect); the important      //
// SCT improvements concern the secondary energy/angle distributions,         //
// whose Doppler widths differ by the full factor of r.                       //
//                                                                            //
// REFERENCES:                                                                //
//  [1] W. E. Lamb, "Capture of Neutrons by Atoms in a Crystal",              //
//      Phys. Rev. 55 (1939) 190.                                             //
//  [2] G. M. Borgonovi, "Neutron Scattering Kernels Calculations at          //
//      Epithermal Energies", Gulf General Atomic report GA-9950 (1970).      //
//  [3] "ENDF-6 Formats Manual", report ENDF-102 (BNL), Section 7.4           //
//      (File 7 thermal scattering, SCT recipe and Teff tabulation).          //
//  [4] R. E. MacFarlane, "New Thermal Neutron Scattering Files for           //
//      ENDF/B-VI Release 2", LA-12639-MS (1994).                             //
//  [5] "The NJOY Nuclear Data Processing System, Version 2016",              //
//      LA-UR-17-20093 (the THERMR chapter: SCT is applied whenever the       //
//      required (alpha,beta) fall outside the tabulated File 7 range).       //
//  [6] T. M. Sutton, T. H. Trumbull, C. R. Lubitz, "Comparison of Some       //
//      Monte Carlo Models for Bound Hydrogen Scattering", Proc. M&C 2009,    //
//      Saratoga Springs (2009); see also their CSEWG-Formats November        //
//      2008 presentation "The Short-Collision-Time Approximation in          //
//      ENDF" (definitive GA formula and misprint catalogue).                 //
//  [7] A. Sjolander, Arkiv foer Fysik 14 (1958) 315 (the harmonic            //
//      incoherent/Gaussian framework in which SCT is the short-time          //
//      limit, and the basis of the Teff definition above).                   //
//  [8] X.-X. Cai and T. Kittelmann, J. Comput. Phys. 380 (2019) 400,         //
//      doi:10.1016/j.jcp.2018.11.043 (NCrystal conventions and the           //
//      free-gas model implementation reused here).                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  class SCTXSProvider final : private MoveOnly {
  public:
    //Cross sections in the SCT model described above. Teff>=T is required
    //(a tiny numerical dip below T is clamped), and Teff==T gives results
    //identical to FreeGasXSProvider at T. table_npts is the node count
    //of the correction lookup table (worst-case relative accuracy of
    //the correction scales as npts^-4, 90 nodes giving ~1e-5;
    //normally provided via SABCfg::Cfg::sct_table_npts, i.e. tied to
    //the knllux/sablux luxury setting):
    SCTXSProvider( Temperature, Temperature teff, AtomMass, SigmaBound,
                   unsigned table_npts = 90 );
    ~SCTXSProvider();

    CrossSect crossSection( NeutronEnergy ekin ) const;

    Temperature effectiveTemperature() const noexcept { return m_teff; }

    //The free-gas beta-kernel I(c,beta) in closed erf form (beta>=0 only;
    //value is independent of alpha convention). Exposed for testing:
    static double evalFGBetaKernel( double c, double beta, double invA );

    //The upscatter correction integral, integral_{beta>0} of
    //I(c,beta)*(1-exp(-(r-1)*beta)), as used in production (asymptotic
    //form, splined lookup table, or quadrature, as appropriate) and
    //evaluated directly by quadrature. Exposed for testing:
    double evalUpscatterCorr( double c ) const;
    double evalUpscatterCorrIntegral( double c, double prec = 1e-9,
                                      unsigned maxlvl = 10 ) const;

  private:
    FreeGasXSProvider m_xsprovider;//free-gas at Teff
    Temperature m_teff;
    double m_r;//Teff/T (>=1.0, exactly 1.0 means pure free-gas behaviour)
    double m_invkTeff;
    double m_invA;
    double m_sb;
    double m_csat;//above this, corr = (r-1)/r to double precision
    double m_clow;//below this (rare), fall back to direct quadrature
    double m_lutlna = 0.0, m_lutlnb = 0.0;
    SplinedLookupTable m_corrlut;//corr vs ln(c), ln(c) in [lna,lnb]
    double corrIntegralFixedOrder( double c ) const;
  };

  class SCTSampler final : private MoveOnly {
  public:
    //Sample (alpha,beta) in the SCT model, in NCrystal's alpha convention
    //and units of the actual kT. For Teff==T the results are identical
    //(sample for sample) to FreeGasSampler::sampleAlphaBeta at T:
    SCTSampler( NeutronEnergy, Temperature, Temperature teff, AtomMass );
    ~SCTSampler();

    PairDD sampleAlphaBeta( RNG& ) const;

  private:
    FreeGasSampler m_fgs;//at Teff
    double m_r;//Teff/T
  };

}

#endif

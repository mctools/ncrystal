#ifndef NCrystal_FreeGasUtils_hh
#define NCrystal_FreeGasUtils_hh

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/phys_utils/NCKinUtils.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCMath.hh"

////////////////////////////////////////////////////////////////////////////////////
//                                                                                //
// Helper classes for getting scattering cross-sections (FreeGasXSProvider) or    //
// sampling scatterings (FreeGasSampler) in the free-gas model, where a neutron   //
// scatters on targets for which all internal energy states are ignored and which //
// are not interacting among themselves - except as required to maintain an       //
// overall thermal velocity balance of the system. The results are based on the   //
// free-gas scattering kernel, which has a close-form expression for              //
// S(alpha,beta) (see for instance Squires chapter 4.5 for the derivation of the  //
// S(q,omega) form, or equation 19 in the NCrystal sampling paper                 //
// (10.1016/j.jcp.2018.11.043) for the S(alpha,beta) form:                        //
//                                                                                //
// S(alpha,beta) = exp( -beta^2/(4alpha/A) - alpha/4A ) / sqrt(alpha)             //
//                                                                                //
// With kinematically allowed region for a given neutron beta>-ekin/kt and alpha  //
// in [alpha-(beta),alpha+(beta)] (see NCKinUtils.cc for the formula).            //
//                                                                                //
// While the integration of this expressions yields a relatively simple           //
// expression for the cross-section, the sampling of a given scattering is        //
// rather challenging - not the least since the code is expected to run with a    //
// different value of E in each sampling, precluding expensive E-dependent        //
// pre-calculations.                                                              //
//                                                                                //
// Firstly, to sample beta one obtains a very complicated expression for          //
// P(beta|E), which is also numerically unstable and expensive to                 //
// evaluate. Next, the shape of P(alpha|E,beta) follows directly from             //
// S(alpha,beta), but sampling it is also challenging due both to the form (in    //
// particular what is essentially a ~exp(-1/alpha) factor), and due to the        //
// varying limits [alpha-(beta),alpha+(beta)], meaning that different features    //
// of the shape will be important for different (E,beta).                         //
//                                                                                //
// The resulting implementation should provide correct results for any energy up  //
// to at least 1GeV neutrons (far beyond the validity of the free gas model       //
// anyway), and allows beta+alpha sampling at a rate roughly around 1MHz.         //
//                                                                                //
////////////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  class FreeGasXSProvider final {
  public:
    FreeGasXSProvider( Temperature, AtomMass, SigmaFree );
    FreeGasXSProvider( Temperature, AtomMass, SigmaBound );
    ~FreeGasXSProvider();

    //Get the cross-section:
    CrossSect crossSection( NeutronEnergy ekin ) const;

    //Evaluate (1+1/(2a^2))*erf(a)+exp(-a^2)/(sqrt(pi)*a) (used internally, but
    //exposed here for testing):
    static double evalXSShapeASq(double a_squared);

    SigmaFree sigmaFree() const { return SigmaFree{ m_sigmaFree }; }
  private:
    double m_sigmaFree, m_ca;
  };

  class FreeGasSampler final {
  public:

    //Evaluate quantum mechanical Free Gas model.

    FreeGasSampler(NeutronEnergy, Temperature, AtomMass );
    ~FreeGasSampler();

    //Beta/energy-transfer sampling:
    double sampleDeltaE(RNG&) const;
    double sampleBeta(RNG&) const;

    //Sample alpha (equivalent to sampling q^2 or scattering angle) for a given beta :
    double sampleAlpha(double beta, RNG&) const;

    //Combined sampling of (alpha,beta) or (delta_ekin,mu=cos(theta_scat)):
    PairDD sampleAlphaBeta( RNG& ) const;
    PairDD sampleDeltaEMu( RNG& ) const;

    //Exposed for testing purposes only:
    void testBetaDistEval ( double beta, double & f_exact, double & f_lb, double & f_ub );

  private:
    double m_c, m_kT, m_sqrtAc, m_invA, m_Adiv4, m_normfact, m_c_real;
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////


namespace NCRYSTAL_NAMESPACE {

  inline CrossSect FreeGasXSProvider::crossSection( NeutronEnergy ekin ) const
  {
    return CrossSect{ m_sigmaFree * evalXSShapeASq( m_ca * ekin.dbl() ) };
  }

  inline double FreeGasSampler::sampleDeltaE( RNG& rng ) const
  {
    return sampleBeta(rng)*m_kT;
  }

  inline PairDD FreeGasSampler::sampleAlphaBeta( RNG& rng ) const
  {
    double beta = sampleBeta(rng);
    if ( beta < -m_c || muIsotropicAtBeta(beta,m_c) ) {
      nc_assert( beta >= -m_c_real*1.001 );
      //close to kinematical end-point, or neutron has such an extreme energy
      //that it triggered various code-paths resulting in -m_c_real<beta<-m_c.
      auto alim = getAlphaLimits( m_c_real, beta );
      double alpha = alim.first+rng.generate()*(alim.second-alim.first);
      return std::make_pair(ncclamp(alpha,alim.first,alim.second), beta);
    }
    return std::make_pair(sampleAlpha(beta,rng),beta);
  }

  inline PairDD FreeGasSampler::sampleDeltaEMu( RNG& rng ) const
  {
    double beta = sampleBeta(rng);
    if ( beta <= -m_c || muIsotropicAtBeta(beta,m_c) ) {
      nc_assert( beta >= -m_c_real*1.001 );
      //close to kinematical end-point, or neutron has such an extreme energy
      //that it triggered various code-paths resulting in -m_c_real<beta<-m_c.
      return std::make_pair(beta*m_kT,rng.generate()*2.0-1.0);
    }
    auto res = convertAlphaBetaToDeltaEMu( sampleAlpha(beta,rng),
                                           beta,
                                           NeutronEnergy{m_c*m_kT},
                                           m_kT );
    return { res.deltaE, res.mu };
  }

}

#endif

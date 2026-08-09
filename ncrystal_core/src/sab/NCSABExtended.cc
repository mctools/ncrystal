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

#include "NCrystal/internal/sab/NCSABExtended.hh"
namespace NC = NCrystal;
namespace NCS = NCrystal::SABUtils;

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {
    namespace {
      inline double sIntegralFromXS( CrossSect xs,
                                     NeutronEnergy ekin,
                                     double invkT )
      {
        return xs.dbl() * ekin.dbl() * 4.0 * invkT;
      }

    }
  }
}


NCS::SABExtended::SABExtended( shared_obj<const SABProcessor> p,
                               shared_obj<const SAB::SABExtender> e )
  : m_p( std::move(p) ),
    m_e( std::move(e) ),
    m_kT( m_p->kT() ),
    m_invkT( 1.0 / m_kT )
{
  nc_assert(!m_p->getEDivKTGrid().empty());
  nc_assert(!m_p->getPhaseSpaceIntegralAtGrid().empty());
  auto emaxInfo = m_p->getEMaxInfo();
  m_emax = emaxInfo.ekin;
  //Initialise constants needed for high-E extrapolation (see comments below
  //where they are being used):
  const auto extenderXS_emax = m_e->crossSection(m_emax);
  const auto tableXS_emax = emaxInfo.crossSectionUnitSigmaBound;
  m_kExtension = ( tableXS_emax.dbl() - extenderXS_emax.dbl() ) * m_emax.dbl();
  m_extSAtEmax = sIntegralFromXS( extenderXS_emax, m_emax, m_invkT );
  m_tableSAtEmax = emaxInfo.phaseSpaceIntegral;
}

NC::CrossSect
NCS::SABExtended::xsHighE( NeutronEnergy ekin ) const
{
  nc_assert( ekin > m_emax );
  //  integral_E(S) = ( tableintegral_Emax(S)
  //                    - extenderintegral_Emax(S) )
  //                  + extenderintegral_E(S)
  //  Now, in general XS(E) = [C/E] * integral_E(S),   C=sigmaB*kT/4. So:
  //    XS_E = [C/E] * integral_E(S)
  //            = [Emax/E]*( [C/Emax]*tableintegral_Emax(S)
  //                         - [C/Emax]*extenderintegral_Emax(S) )
  //               + [C/E]*extenderintegral_E(S)
  //            = [Emax/E] *(tableXS_Emax-extenderXS_Emax) + extenderXS_E
  //            = k / E + extenderXS_E
  return CrossSect{ m_kExtension / ekin.dbl()
                    + m_e->crossSection( ekin ).dbl() };
}

NC::PairDD NCS::SABExtended::scatABHighE( RNG& rng, NeutronEnergy ekin ) const
{
  //Sample (alpha,beta) for ekin>m_emax. Inside the kinematic curve
  //corresponding to E=m_emax, S(alpha,beta) is given by the tabulated kernel
  //(captured in the SABProcessor), and outside the curve it is modelled by the
  //S(alpha,beta) values represented by the SABExtender.

  nc_assert( ekin >= m_emax );

  //////////////////////////////////////////////////////////////////////////
  // First a little safeguard to prevent long tails in sampling times of  //
  // events just above emax.                                              //
  //////////////////////////////////////////////////////////////////////////

  //Safe-guard against evens at emax+epsilion needing O(1/epsilon) free gas
  //sampling attempts in the rare case where the event should fall outside the
  //emax boundary. This is on average a small O(10%) slowdown for events with
  //ekin in [emax,1*1,emax], but this is preferable to having rare events with
  //truly abysmal efficiency.

  const NeutronEnergy near_emax{ m_emax.dbl() * 1.1};//1.1 => slowdown is
                                                     //O(10%), but tail events
                                                     //are kept at
                                                     //O(1/(1.1-1))=O(10)
                                                     //samplings.
  if ( ekin < near_emax ) {
    const double four_E_div_kT = 4.0 * ekin.dbl() * m_invkT;
#  ifndef NDEBUG
    int nloop = 0;
#  endif
    while ( true ) {
      nc_assert( nloop++ < 1000 );
      auto ab = scatABHighE( rng, near_emax );
      if ( ncsquare( ab.first - ab.second ) <= four_E_div_kT*ab.first )
        return ab;
    }
  }

  //Ok, we are well above emax:
  nc_assert( ekin.dbl() >= near_emax.dbl() );

  /////////////////////////////////////////////////////////////////////////
  // Figure out and do a random choice against the probability of        //
  // sampling (alpha,beta) point inside the E=m_emax kinematic boundary. //
  // This probability is given by comparing the relative S-integrals     //
  // (which are proportional to sigma(E)*E).                             //
  /////////////////////////////////////////////////////////////////////////

  {
    //xs = sigma_bound * sintegral / 4E, so with sigma_bound==1barn (as we
    //require here), we get sintegral = 4*E*xs
    const double extS = sIntegralFromXS(m_e->crossSection(ekin), ekin, m_invkT);
    nc_assert( extS >= m_extSAtEmax );
    const double extContrib = (extS - m_extSAtEmax);
    if ( rng() * ( m_tableSAtEmax+extContrib) <= m_tableSAtEmax ) {
      //results end up in central region, covered by the SABProcessor:
      auto ab = m_p->sampleScatterAlphaBeta( rng, m_emax );
      return { ab.alpha, ab.beta };
    }
  }

  //////////////////////////////////////////////////////////////////////////
  // Ok, we need to sample with the SABExtender, and only accepts results //
  // not within the m_emax phasespace.                                    //
  //////////////////////////////////////////////////////////////////////////

#ifndef NDEBUG
  int nloop2(0);
#endif
  const double fouremax_div_kT = 4*m_emax.dbl()*m_invkT;
  while ( true ) {
    nc_assert(nloop2++<1000);
    auto ab = m_e->sampleAlphaBeta(rng,ekin);
    if ( !( ncsquare( ab.first - ab.second ) < fouremax_div_kT*ab.first ) )
      return ab;
  }
}

NC::ScatterOutcomeIsotropic
NCS::SABExtended::sampleScatter( RNG& rng, NeutronEnergy ekin ) const
{
  auto ab = sampleScatterAlphaBeta( rng, ekin );

  //Fixme: this code is duplicated in too many places (and is it even
  //technically correct that we end up isotropic??)
  double dE, mu;
  if ( muIsotropicAtBeta(ab.second,ekin.get()*m_invkT) ) {
    //close to kinematical end-point, numerically safe fall-back:
    dE = ab.second*m_kT;
    mu = rng.generate()*2.0 - 1.0;
  } else {
    auto dEmu = convertAlphaBetaToDeltaEMu( ab.first, ab.second, ekin, m_kT );
    dE = dEmu.deltaE;
    mu = dEmu.mu;
  }
  return { NeutronEnergy{ncmax(0.0,ekin.dbl()+dE)}, CosineScatAngle{mu} };
 }

NC::shared_obj<const NC::SABUtils::SABExtended>
NC::SABUtils::
SABExtended::createWithFGExtender( shared_obj<const SABProcessor> processor )
{
  //fixme: inconsistent namespaces
  auto sab = processor->sabDataPtr();
  auto ext = makeSO<SAB::SABFGExtender>( sab->temperature(),
                                         sab->elementMassAMU(),
                                         SigmaBound{1.0} );
  return makeSO<SABUtils::SABExtended>( std::move(processor), std::move(ext) );
}


NC::shared_obj<const NC::SABUtils::SABExtended>
NC::SABUtils::
SABExtended::createWithFGExtender( const SABCfg::Cfg& cfg,
                                   shared_obj<const SABData> sab,
                                   std::shared_ptr<const VectD> egrid )
{
  return
    createWithFGExtender( makeSO<SABUtils::SABProcessor>( cfg,
                                                          std::move(sab),
                                                          std::move(egrid) ) );
}

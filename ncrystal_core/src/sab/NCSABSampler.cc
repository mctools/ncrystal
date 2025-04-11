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

#include "NCrystal/internal/sab/NCSABSampler.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;

NC::SABSampler::~SABSampler() = default;

NC::SABSampler::SABSampler( Temperature temperature,
                            VectD&& egrid,
                            SABSamplerAtEList&& samplers,
                            std::shared_ptr<const SAB::SABExtender> extender,
                            double xsAtEmax,
                            EGridMargin egridMargin )
{
  setData( temperature, std::move(egrid), std::move(samplers), std::move(extender), xsAtEmax, egridMargin );
}


void NC::SABSampler::setData( Temperature temperature,
                              VectD&& egrid,
                              SABSamplerAtEList&& samplers,
                              std::shared_ptr<const SAB::SABExtender> extender,
                              double xsAtEmax,
                              EGridMargin egridMargin )
{
  m_egrid = std::move(egrid);
  m_samplers = std::move(samplers);
  m_kT = temperature.kT();
  m_extender = std::move(extender);
  m_xsAtEmax = xsAtEmax;
  m_k1 = m_xsAtEmax * m_egrid.back();
  m_k2 = m_extender->crossSection( NeutronEnergy{m_egrid.back()} ).dbl() * m_egrid.back();
  nc_assert_always( m_egridMargin.value >= 1.0 && m_egridMargin.value < 1e3 );
  m_egridMargin = egridMargin;
}

NC::PairDD NC::SABSampler::sampleHighE(NeutronEnergy ekin, RNG& rng) const
{
  const double emax = m_egrid.back();
  nc_assert( ekin.get() >= emax );
  //Sample (alpha,beta) using provided SABExtender. A returned alpha value of
  //-1.0 indicates that the usual code should sample (alpha,beta) from the
  //tabulated kernel with ekin=m_egrid.back().

  //Principle: Inside the kinematic curve corresponding to
  //E=Emax=m_egrid.back(), S(alpha,beta) is given by the tabulated kernel, and
  //outside the curve it is modelled by the S(alpha,beta) values represented by
  //m_extender.

  //Simply speaking, the strategy is to go ahead and sample via the extender
  //(using E=ekin>Emax), and then resample any S(alpha,beta) values that
  //happen to fall inside the reach of the tabulated data, i.e. inside the
  //kinematic curve associated with Emax, with the tabulated data (using now
  //ekin=Emax). However, in general the cross-sections (and therefore the
  //S-integrals) at Emax for the tabulated data and m_extender are *not*
  //exactly identical, and therefore one must correct for this to avoid
  //biasing the distributions.

  //At a given E>Emax, the
  //probability for sampling a point inside the Emax curve should be given by:
  //
  //  P_inside = table_Sintegral_inside_Emax_curve / ( table_Sintegral_inside_Emax_curve + extender_Sintegral_between_Emax_curve_and_E_curve)
  //
  //Using sigma(E) = [C/E]*Sintegral_inside_E_curve, we can evaluate this in terms of cross-sections:
  //
  //  S_integral_inside_E_curve \propto E*sigma(E)
  //
  //Thus:
  //
  //  P_inside = table_XS(Emax)*Emax / ( table_XS(Emax) * Emax + extenderXS(E)*E-extenderXS(Emax)*Emax)
  //           = table_XS(Emax)*Emax / ( [table_XS(Emax)-extenderXS(Emax)]*Emax + extenderXS(E)*E)

  const double extenderXSMultE = ekin.get() * m_extender->crossSection(ekin).dbl();
  const double P_inside = m_k1 / ( (m_k1-m_k2) + extenderXSMultE );

  // But sampling from extender only will give (alpha,beta) values inside the
  // Emax-curve with different probability:
  //
  //   P_extender_inside = extender_Sintegral_inside_Emax_curve / extender_Sintegral_inside_E_curve
  //
  // Which once again can be evaluated with cross sections:
  //
  //   P_extender_inside = extender_Sintegral_inside_Emax_curve / extender_Sintegral_inside_E_curve
  //                       [ extenderXS(Emax)*Emax ] / [ extenderXS(E)*E ]

  const double P_extender_inside = m_k2 / extenderXSMultE;

  // If P_extender_inside>=P_inside, we sample from extender and throw away sampled
  // values inside the Emax kinematics curve P_discardinside of the time. When a value inside
  // the Emax curve is not thrown away, it is resambled from the table with ekin=Emax.
  //
  const double P_discardinside = ( P_extender_inside >= P_inside ? (1.0-P_inside/P_extender_inside) : 0.0 );

  if ( P_discardinside > 0.95 )
    NCRYSTAL_THROW(BadInput,"Scattering Kernel does not appear to match up very well with the chosen extrapolation model at Emax.");

  // If on the other hand, P_extender_inside<P_inside, then we must boost the
  // number of values sampled inside emax. P_extrainside of the time we can go
  // directly and sample the table with ekin=Emax and return the value to the
  // caller. Otherwise we do as when P_extender_inside>=P_inside, but with
  // P_discardinside=0.

  if (P_extender_inside<P_inside) {
    // P_extrainside can be determined via the requirement:
    //
    //      P_extrainside + (1-P_extrainside)*P_extender_inside = P_inside
    //  <=> P_extrainside * (1-P_extender_inside) = (P_inside-P_extender_inside)
    //
    // If 1-P_extender_inside is not extremely close to zero, this is easy to
    // evaluate. Otherwise (and this should normally only happen when E~=Emax), it
    // is probably fine to put P_extrainside=1.0.
    double aa = 1.0-P_extender_inside;
    const double P_extrainside = aa>1e-10 ? (P_inside-P_extender_inside)/aa : 1.0;
    if (rng.generate()<P_extrainside)
      return {-1.0,0.0};//must sample tabulated data with ekin=emax
  }

  const double emax_div_kt = emax/m_kT;
  while (true) {
    //sample with extender:
    auto alphabeta = m_extender->sampleAlphaBeta(rng,ekin);
    //if outside emax curve, always return immediately:
    if ( alphabeta.second <= -emax_div_kt )
      return alphabeta;
    auto alims =  getAlphaLimits( emax_div_kt,alphabeta.second );
    if (!valueInInterval(alims.first,alims.second,alphabeta.first))
      return alphabeta;
    //inside emax curve. Check P_discardinside:
    if ( P_discardinside && rng.generate()<P_discardinside )
      continue;//discard to avoid bias
    //Sample with the table using ekin=emax!
    return {-1.0,0.0};
  }
}

NC::PairDD NC::SABSampler::sampleAlphaBeta(NeutronEnergy ekin, RNG& rng) const
{
  nc_assert( m_egrid.size()>1 && m_egrid.size()==m_samplers.size() );
  double alpha,beta;

  decltype(m_samplers.begin()) itSampler;

  auto itEkinUpper = std::upper_bound ( m_egrid.begin(), m_egrid.end(), ekin.dbl() );

  bool ultra_small_ekin_mode = false;
  const double ultra_small_ekin = m_egrid.front();

  if ( itEkinUpper == m_egrid.end() ) {

    //High-E extrapolation via m_extender.
    auto alphabeta =  sampleHighE(ekin, rng);
    if (alphabeta.first>=0.0)
      return alphabeta;
    //HighE code decided that we must sample the kernel with ekin=emax:
    ekin = NeutronEnergy{m_egrid.back()};
    itSampler = std::prev(m_samplers.end());

  } else if ( itEkinUpper == m_egrid.begin() ) {

    //Low-E extrapolation. Beta-distribution is essentially unchanged at this
    //energy, but must treat alpha-sampling specially.
    itSampler = m_samplers.begin();
    ultra_small_ekin_mode = (ekin.get()<ultra_small_ekin);

  } else {
    //Inside range of energy grid. Apply egridMargin if relevant+possible and pick overlay sampler.
    if ( m_egridMargin.value > 1.0 ) {
      while ( std::next(itEkinUpper) != m_egrid.end() && ekin.dbl()*m_egridMargin.value > *itEkinUpper )
        ++itEkinUpper;
    }
    itSampler = m_samplers.begin()+std::distance(m_egrid.begin(), itEkinUpper);
  }

  //Sample with *itSampler, using the rejection method to make the results
  //correct at the given ekin value:
  const double ekin_div_kT = ekin.get()/m_kT;
  const double sampling_ekin_div_kT = (ultra_small_ekin_mode ? ultra_small_ekin/m_kT : ekin_div_kT);
  int loopmax(100);
  while (loopmax--) {
    std::tie(alpha,beta) = (*itSampler)->sampleAlphaBeta(sampling_ekin_div_kT, rng);
    if (beta<-ekin_div_kT)
      continue;
    double alow,aupp;
    {
      auto alims = getAlphaLimits( ekin_div_kT, beta );
      alow = alims.first;
      aupp = alims.second;
    }
    if (valueInInterval(alow,aupp,alpha))
      return { alpha, beta };
    if (ultra_small_ekin_mode) {
      //energy is very low, so |alpha+ - alpha-| is also very low, and
      //S(alpha+,beta)~=S(alpha-,beta) (or at least we could approximate S with
      //a straight line over the interval. To avoid wasting time on potentially
      //abysmal acceptance rates, we simply generate an isotropic alpha.
      //
      //TODO: the auto-emin-determination code finds the point where S has
      //linear behaviour, not necessarily flat, so it would be consistent to
      //sample with linear interpolation of S from (alpha-,beta) to
      //(alpha+,beta).
      return { alow + rng.generate()*(aupp-alow), beta };
    }
  }
  NCRYSTAL_THROW2(CalcError,"Infinite looping in sampleAlphaBeta(ekin="<<ekin<<")");
}

NC::PairDD NC::SABSampler::sampleDeltaEMu(NeutronEnergy ekin, RNG& rng) const
{
  auto alphabeta = sampleAlphaBeta(ekin,rng);
  if ( NC::muIsotropicAtBeta(alphabeta.second,ekin.get()/m_kT) )
    return std::make_pair( alphabeta.second*m_kT, rng.generate()*2.0 - 1.0 );
  auto res = convertAlphaBetaToDeltaEMu(alphabeta,ekin,m_kT);
  return { res.deltaE, res.mu };
}

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

#include "NCrystal/internal/sab/NCSABXSProvider.hh"

namespace NC = NCrystal;

NC::SABXSProvider::~SABXSProvider() = default;

NC::SABXSProvider::SABXSProvider( VectD&& egrid,
                                  VectD&& xsvals,
                                  std::shared_ptr<const SAB::SABExtender> extender )
  : m_kExtension(0.0)
{
  setData(std::move(egrid),std::move(xsvals),std::move(extender));
}

void NC::SABXSProvider::setData( VectD&& egrid,
                                 VectD&& xsvals,
                                 std::shared_ptr<const SAB::SABExtender> extender )
{
  m_egrid = std::move(egrid);
  m_xs = std::move(xsvals);
  m_extender = std::move(extender);
  nc_assert_always(!!m_extender);
  nc_assert_always(!m_egrid.empty());
  nc_assert_always(!m_xs.empty());

  const double emax = m_egrid.back();
  const double extenderXS_emax = m_extender->crossSection(NeutronEnergy{emax}).dbl();
  const double tableXS_emax = m_xs.back();
  //constant needed for high-E extrapolation (see comments below where it is
  //being used):
  m_kExtension = ( tableXS_emax - extenderXS_emax ) * emax;
}

NC::CrossSect NC::SABXSProvider::crossSection( NeutronEnergy ekin ) const
{
  nc_assert( ! m_xs.empty() && m_xs.size() == m_egrid.size() );

  auto itEkinUpper = std::upper_bound ( m_egrid.begin(), m_egrid.end(), ekin.dbl() );
  if ( itEkinUpper == m_egrid.end()) {
    //  integral_E(S) = (tableintegral_Emax(S)-extenderintegral_Emax(S))+extenderintegral_E(S)
    //  Now, in general XS(E) = [C/E] * integral_E(S),   C=sigmaB*kT/4. So:
    //    XS_E = [C/E] * integral_E(S)
    //            = [Emax/E]*([C/Emax]*tableintegral_Emax(S)-[C/Emax]*extenderintegral_Emax(S))+[C/E]*extenderintegral_E(S)
    //            = [Emax/E] *(tableXS_Emax-extenderXS_Emax) + extenderXS_E = k / E + extenderXS_E
    return CrossSect{ m_kExtension / ekin.dbl() + m_extender->crossSection( ekin ).dbl() };
  } else if ( itEkinUpper == m_egrid.begin() ) {

    //Energy is below lowest tabulated energy. At very small energies, the
    //kinematically allowed region of (alpha,beta) space becomes ever thinner,
    //with |alpha+(beta)-alpha-(beta)| ~= 4*sqrt(E/kT)*sqrt(beta). As E->0, it
    //will therefore eventually be true that S(alpha|beta) will be approximately
    //flat (even just a linear function of alpha) within [alpha-,alpha+]. Thus
    //the entire integral of S(alpha,beta) will eventually decrease as sqrt(E)
    //(the beta integral will run from [-E/kT,inf] which is essentially the same
    //[0,inf] for all "very small" energies). In addition to the integral over
    //S, sigma(E) also includes a factor proportional to 1/E. Thus, sigma(E)
    //will decrease as 1/sqrt(E) for small energies (we have thus essentially
    //derived, or at least argued for, the "1/v law").

    return CrossSect { ekin.dbl() > 0.0 ? std::sqrt( m_egrid.front() / ekin.dbl() ) * m_xs.front() : kInfinity };
  } else {
    //linear interpolation in grid
    auto itEkinLower = std::prev(itEkinUpper);
    auto itXSLower = m_xs.begin() + std::distance(m_egrid.begin(), itEkinLower);
    auto itXSUpper = std::next(itXSLower);
    const double dXS = *itXSUpper - *itXSLower;
    const double dEkin = *itEkinUpper - *itEkinLower;
    nc_assert(dEkin>0.0);
    double xs = *itXSLower + dXS * ( ekin.dbl() - *itEkinLower ) / dEkin;
    nc_assert(xs>=std::min<double>(*itXSLower,*itXSUpper));
    nc_assert(xs<=std::max<double>(*itXSLower,*itXSUpper));
    return CrossSect{ xs };
  }
}


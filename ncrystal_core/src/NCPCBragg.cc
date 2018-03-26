////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

#include "NCrystal/NCPCBragg.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCException.hh"
#include "NCMath.hh"
#include <functional>//std::greater
#include <algorithm>//std::upper_bound, std::lower_bound

NCrystal::PCBragg::PCBragg(const Info*ci)
  : ScatterIsotropic("PCBragg"),
    m_threshold_wl(0.0),
    m_threshold_ekin(wl2ekin(0.0)),
    m_xsectfact(-1.0)
{
  nc_assert_always(ci);

  //We always require HKL info, since it is needed for both crossSection(..),
  //threshold(), and shootScatterAngle().
  if (!ci->hasHKLInfo())
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks HKL information.");

  m_threshold_wl = ci->nHKL() ? ci->hklBegin()->dspacing * 2.0 : 0.0;
  m_threshold_ekin = wl2ekin(m_threshold_wl);

  if (ci->hasStructureInfo()) {
    nc_assert_always(ci->getStructureInfo().n_atoms>0);
    m_xsectfact = 0.5/(ci->getStructureInfo().volume * ci->getStructureInfo().n_atoms);
    nc_assert_always(m_xsectfact>0);
  }

  m_2d.reserve(ci->nHKL());
  m_fdm_commul.reserve(ci->nHKL());
  HKLList::const_iterator it = ci->hklBegin();
  HKLList::const_iterator itE = ci->hklEnd();
  double fdmsum(0.0);
  for (;it!=itE;++it) {
    if (it->dspacing<0)
      NCRYSTAL_THROW(CalcError,"Inconsistent HKL data implies negative d_spacing.");
    double c = it->fsquared * it->multiplicity * it->dspacing;
    if (c<0)
      NCRYSTAL_THROW(CalcError,"Inconsistent HKL data implies negative |F|^2*multiplicity.");
    m_fdm_commul.push_back(fdmsum += c);
    m_2d.push_back(2.0*it->dspacing);
  }
  if (m_fdm_commul.empty()||m_fdm_commul.back()<=0.0) {
    m_threshold_wl = 0.0;
    m_threshold_ekin = wl2ekin(m_threshold_wl);
  }

  if ( m_fdm_commul.size()>1) {
    //Merge entries with identical d-spacing for efficiency (could be skipped if
    //we need to track which original HKL entry triggered a scattering).
    std::vector<double> newfdm, new2d;
    size_t n = m_fdm_commul.size();
    nc_assert(m_2d.size()==n);
    for (size_t i = 0; i<n; ++i) {
      double dd = m_2d.at(i);
      //look at next entries for identical d-spacings:
      while ( i+1<n && ncabs(m_2d.at(i+1)-dd)<1e-10) {
        ++i;
      }
      newfdm.push_back(m_fdm_commul.at(i));
      new2d.push_back(dd);
    }
    std::vector<double>(newfdm.begin(),newfdm.end()).swap(m_fdm_commul);
    std::vector<double>(new2d.begin(),new2d.end()).swap(m_2d);
  }
  validate();
}

NCrystal::PCBragg::~PCBragg()
{
}

double NCrystal::PCBragg::crossSectionNonOriented(double ekin) const
{
  if (ekin<m_threshold_ekin)
    return 0.0;

  double wl = ekin2wl(ekin);
  //TODO for NC2: we could do this directly on energy, no need to convert energy
  //to wl!  Meaning that m_2d should actually store the energy equivalent of
  //2*d, and genSinThetaBragg should use energy when looking up there.

  if (m_xsectfact<0)
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks Structure information needed for crossSection(..).");

  //Use binary search to find iterator to first hkl plane which has 2*dspacing <
  //neutron_wavelength, and use the resulting index of the iterator to look up
  //the summed cross-section contribution of all hkl planes before it:
  std::vector<double>::const_iterator it = std::upper_bound(m_2d.begin(),m_2d.end(),wl,std::greater<double>());
  unsigned idx=(unsigned)(it-m_2d.begin());
  if (idx==0)
    return 0.0;//extra check, but should only really happen below m_threshold

  nc_assert(idx<=m_fdm_commul.size());
  return m_fdm_commul[idx-1] * m_xsectfact * wl * wl;
}

void NCrystal::PCBragg::genSinThetaBragg(double wl,double& sinthetabragg) const
{
  //generates sin(theta_Bragg). A value of 2 signals that Bragg scattering is not possible.

  unsigned idx(0);
  if (wl>m_threshold_wl) {
    idx = 0;//no plane
  } else {
    //Use binary search to find iterator to first hkl plane which has 2*dspacing <
    //neutron_wavelength, and use the resulting index of the iterator to look up
    //the summed cross-section contribution of all hkl planes before it:
    std::vector<double>::const_iterator it = std::upper_bound(m_2d.begin(),m_2d.end(),wl,std::greater<double>());
    idx=(unsigned)(it-m_2d.begin());
  }

  if (idx==0) {
    //Bragg condition not satisfied => no scattering from any plane.
    //This should not happen except if the wavelength is above threshold_wl, or
    //perhaps in case of errors due to floating point precision when near the
    //threshold. So for lack of a better option, we indicate to that calling
    //routine that it should pick an isotropic scattering.
    sinthetabragg = 2;
    return;
  }

  nc_assert(idx<=m_fdm_commul.size());

  //select a plane at random:
  double rand_contrib = this->rand() * m_fdm_commul[idx-1];

  //Use binary search again, this time to find which plane was selected at random:
  std::vector<double>::const_iterator it2 = std::lower_bound(m_fdm_commul.begin(),m_fdm_commul.begin()+(idx-1),rand_contrib);
  unsigned idx_rand = (unsigned)(it2-m_fdm_commul.begin());
  if (idx_rand>=m_2d.size())
    NCRYSTAL_THROW(CalcError,"Inconsistency encountered while selecting scattering plane.");

  sinthetabragg = wl / m_2d[idx_rand];
}

void NCrystal::PCBragg::generateScatteringNonOriented( double ekin, double& angle, double& dekin ) const
{
  dekin = 0;//strictly elastic
  double x;
  genSinThetaBragg(ekin2wl(ekin),x);
  if (x==2) {
    angle = randIsotropicScatterAngle();
    return;
  }

#if 0
  //if asin is faster than acos (gcc 4.9.2 on fedora 21 indicates acos is faster):
  angle = 2.0 * asin( x );
#else
  //if acos is faster than asin, using the identity cos(2*asin(x)) = 1-2*x^2
  //(useful also if we ever decide to return cosangle rather than angle):
  angle = acos( 1.0 - 2.0 * x * x );
#endif
}

void NCrystal::PCBragg::generateScattering( double ekin, const double (&indir)[3],
                                            double (&outdir)[3], double& dekin ) const
{
  //Reimplement generateScattering since we can avoid expensive trigonometric
  //function calls entirely.

  dekin = 0;//strictly elastic
  double x;
  genSinThetaBragg(ekin2wl(ekin),x);
  if (x==2) {
    randIsotropicDirection(outdir);
    return;
  }

  //cheap cos+sin of angle with a few handy trigonometric identities:
  double x2 = x*x;
  double cosangle = 1.0 - 2.0 * x2;
  double sinangle = 2.0 * x * sqrt( ncabs(1.0 - x2) );//abs to protect against rounding errors
  randDirectionGivenScatterAngle(cosangle,sinangle,indir,outdir);
}

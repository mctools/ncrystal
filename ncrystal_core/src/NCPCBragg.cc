////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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

#include "NCrystal/internal/NCPCBragg.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/internal/NCPlaneProvider.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include <functional>//std::greater

void NCrystal::PCBragg::init( const StructureInfo& si,
                              std::vector<PairDD >& data )//(dspacing,fsquared_sum)
{
  nc_assert_always(si.n_atoms>0);
  nc_assert_always(si.volume>0);
  if (!(si.volume>0) || !(si.n_atoms>=1) )
    NCRYSTAL_THROW(BadInput,"Passed structure info object has invalid volume or n_atoms fields.");
  init(si.volume * si.n_atoms, data);
}

void NCrystal::PCBragg::init( double v0_times_natoms,
                              std::vector<PairDD >& data )//(dspacing,fsquared_sum)
{
  if (!(v0_times_natoms>0) )
    NCRYSTAL_THROW(BadInput,"v0_times_natoms is not a positive number.");
  double xsectfact = 0.5/v0_times_natoms;
  xsectfact *= wl2ekin(1.0);//Adjust units so we can get cross sections through
                            //multiplication with 1/ekin instead of wl^2.
  std::sort(data.begin(),data.end(),std::greater<PairDD >());
  VectD v2dE;
  v2dE.reserve(data.size());
  VectD fdm_commul;
  fdm_commul.reserve(data.size());
  StableSum fdmsum2;
  std::vector<PairDD >::const_iterator it(data.begin()),itE(data.end());
  double prev_dsp = -kInfinity;
  for (;it!=itE;++it) {
    if (!(it->first>0.0))
      NCRYSTAL_THROW(CalcError,"Inconsistent plane data implies non-positive (or NaN) d_spacing.");
    if (ncabs(prev_dsp-it->first)<1e-11) {
      double c = it->first * it->second * xsectfact;
      fdmsum2.add(c);
      fdm_commul.back() = fdmsum2.sum();
    } else {
      prev_dsp = it->first;
      double c = it->first * it->second * xsectfact;
      fdmsum2.add(c);
      fdm_commul.push_back(fdmsum2.sum());
      v2dE.push_back(wl2ekin(2.0*it->first));
    }
  }
  if (fdm_commul.empty()||fdm_commul.back()<=0.0) {
    m_threshold = kInfinity;
    fdm_commul.clear();
    v2dE.clear();
  } else {
    m_threshold = *(v2dE.begin());
  }
  //Transfer while squeezing memory:
  VectD(fdm_commul.begin(),fdm_commul.end()).swap(m_fdm_commul);
  VectD(v2dE.begin(),v2dE.end()).swap(m_2dE);
  nc_assert(m_threshold>0);
  validate();
}

NCrystal::PCBragg::PCBragg(const StructureInfo& si, PlaneProvider * pp)
  : ScatterIsotropic("PCBragg"),
    m_threshold(kInfinity)
{

  std::vector<std::pair<double, double> > data;
  data.reserve(4096);
  double d,f;
  Vector demi_normal;
  pp->prepareLoop();
  while (pp->getNextPlane(d,f,demi_normal)) {
    if (f<0)
      NCRYSTAL_THROW(CalcError,"Inconsistent data implies negative |F|^2.");
    f*=2;//getNextPlane provides demi-normals, e.g. only half of the normals.
    if (data.empty()||data.back().first!=d) {
      data.emplace_back(d,f);
    } else {
      data.back().second += f;
    }
  }
  init(si,data);
}

NCrystal::PCBragg::PCBragg( const StructureInfo& si, std::vector<PairDD >&  data)
  : ScatterIsotropic("PCBragg"),
    m_threshold(kInfinity)
{
  init(si,data);
}

NCrystal::PCBragg::PCBragg( double v0_times_natoms, std::vector<PairDD >&  data)
  : ScatterIsotropic("PCBragg"),
    m_threshold(kInfinity)
{
  init(v0_times_natoms,data);
}

NCrystal::PCBragg::PCBragg(const Info*ci)
  : ScatterIsotropic("PCBragg"),
    m_threshold(kInfinity)
{
  nc_assert_always(ci);
  if (!ci->hasHKLInfo())
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks HKL information.");
  if (!ci->hasStructureInfo())
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks Structure information.");
  std::vector<std::pair<double, double> > data;
  data.reserve(ci->nHKL());
  HKLList::const_iterator it = ci->hklBegin();
  HKLList::const_iterator itE = ci->hklEnd();
  for (;it!=itE;++it) {
    double f = it->fsquared * it->multiplicity;
    if (f<0)
      NCRYSTAL_THROW(CalcError,"Inconsistent data implies negative |F|^2*multiplicity.");
    if (data.empty()||data.back().first!=it->dspacing) {
      data.emplace_back(it->dspacing,f);
    } else {
      data.back().second += f;
    }
  }
  init(ci->getStructureInfo(),data);
}

NCrystal::PCBragg::~PCBragg()
{
}

void NCrystal::PCBragg::domain(double& ekin_low, double& ekin_high) const
{
  ekin_low = m_threshold;
  ekin_high = kInfinity;
}

std::size_t NCrystal::PCBragg::findLastValidPlaneIdx(double ekin) const {
  //Quick binary search to find index of the plane with the smallest d-spacing
  //satisfying wl<=2d, but in energy-space: Finding the index of the plane with
  //the largest value of ekin2wl(2d) satisfying ekin>=ekin2wl(2d).  We already
  //know that ekin>=m_2dE[0], so we search from one past this entry:
  nc_assert( ekin >= m_threshold );
  return (std::upper_bound(m_2dE.begin() + 1,m_2dE.end(),ekin) - m_2dE.begin()) - 1;
}

double NCrystal::PCBragg::crossSectionNonOriented(double ekin) const
{
  if (ekin<m_threshold)
    return 0.0;
  std::size_t idx = findLastValidPlaneIdx(ekin);
  nc_assert(idx<m_fdm_commul.size());
  return m_fdm_commul[idx] / ekin;
}

double NCrystal::PCBragg::genScatterMu(RandomBase* rng, double ekin) const
{
  nc_assert(ekin>=m_threshold);

  std::size_t idx = findLastValidPlaneIdx(ekin);
  nc_assert(idx<m_fdm_commul.size());

  //randomly select one plane by contribution:
  VectD::const_iterator itFCUpper = m_fdm_commul.begin()+idx;
  VectD::const_iterator itFC = std::lower_bound( m_fdm_commul.begin(),
                                                               itFCUpper,
                                                               rng->generate() * (*itFCUpper) );
  std::size_t idx_rand = (std::size_t)( itFC - m_fdm_commul.begin() );
  nc_assert(idx_rand<m_2dE.size());
  double sin_theta_bragg_squared = m_2dE[idx_rand] / ekin;

  //scatter angle A=2*theta_bragg, so with x=sin^2(theta_bragg), we have:
  //   x = sin^2(A/2)= (1-cosA)/2 => 1-2x = cosA = mu
  const double mu = 1.0 - 2.0 * sin_theta_bragg_squared;
  nc_assert(ncabs(mu)<=1.0);
  return mu;
}

void NCrystal::PCBragg::generateScatteringNonOriented( double ekin, double& angle, double& dekin ) const
{
  dekin = 0;//strictly elastic

  if (ekin<m_threshold) {
    //scatterings not possible here
    angle = 0.0;
  } else {
    angle = std::acos(genScatterMu(getRNG(),ekin));
  }
}

void NCrystal::PCBragg::generateScattering( double ekin, const double (&indir)[3],
                                            double (&outdir)[3], double& dekin ) const
{
  //Reimplement generateScattering to avoid expensive trigonometric function
  //calls.

  dekin = 0;//strictly elastic

  if (ekin<m_threshold) {
    //scatterings not possible here
    outdir[0] = indir[0];
    outdir[1] = indir[1];
    outdir[2] = indir[2];
  } else {
    RandomBase * rng = getRNG();
    randDirectionGivenScatterMu(rng,genScatterMu(rng,ekin),indir,outdir);
  }
}

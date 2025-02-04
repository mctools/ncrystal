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

#include "NCrystal/internal/phys_utils/NCGaussMos.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC=NCrystal;

namespace NCRYSTAL_NAMESPACE {
  inline double GaussMos_cacheRound(double x) {
    //Cut off input at 15 decimals, which should be a negligible effect on any
    //realistic value of wavelength or 1/2dspacing, but ensures that we don't
    //get call-order irreproducibilities. The ncmax fct ensures that results
    //will be >0.0, so the code won't have to worry about FPE's when e.g. 1e-200
    //gets rounded to 0.0.
    return std::floor(ncmax(x,1e-15)*1e15+0.5)*1e-15;
  }
}

NC::GaussMos::GaussMos( MosaicityFWHM mosaicity, double prec, double ntrunc )
  : m_mos_truncN( ntrunc==0 ? GaussOnSphere::estimateNTruncFromPrec(prec) : ntrunc),
    m_prec(prec)
{
  double override_ntrunc = ncgetenv_dbl("GAUSSMOS_OVERRIDE_NTRUNC");
  if (override_ntrunc)
    m_mos_truncN = override_ntrunc;
  nc_assert(prec>=0);
  nc_assert(m_mos_truncN>=0);
  //Set mosaicity and trigger one call to updateDerivedValues:
  setMosaicity(mosaicity);
  nc_assert( m_mos_fwhm.dbl() != -99 );
  nc_assert( m_mos_sigma.dbl() != -99 );
}

NC::GaussMos::GaussMos( MosaicitySigma mosaicity, double prec, double ntrunc )
  : GaussMos( mosaicity.fwhm(), prec, ntrunc )
{
  //Same, but keeping exact MosaicitySigma value:
  m_mos_sigma = mosaicity;
  m_mos_sigma.validate();
}

NC::GaussMos::~GaussMos() = default;

void NC::GaussMos::updateDerivedValues()
{
  double truncangle = m_mos_truncN*m_mos_sigma.dbl();
  if ( ! (truncangle < kPiHalf) )
    NCRYSTAL_THROW(BadInput,"Mosaicity too large, truncation angle (sigma*Ntrunc) must be less than pi/2");
  m_gos.set(m_mos_sigma.dbl(), truncangle, m_prec );
}

void NC::GaussMos::setMosaicity( MosaicityFWHM mosaicity )
{
  mosaicity.validate();
  nc_assert_always(mosaicity.get()>0);
  m_mos_fwhm = mosaicity;
  m_mos_sigma = mosaicity.sigma();
  updateDerivedValues();
}

void NC::GaussMos::setMosaicity( MosaicitySigma mosaicity )
{
  //Same, but keeping exact MosaicitySigma value:
  auto mf = mosaicity.fwhm();
  setMosaicity(mf);
  m_mos_sigma = mosaicity;
}

void NC::GaussMos::setTruncationN(double N)
{
  nc_assert(N>0);
  if (m_mos_truncN != N) {
    m_mos_truncN = N;
    updateDerivedValues();
  }
}

void NC::GaussMos::setPrecision(double p)
{
  nc_assert(p>=0);
  if (m_prec != p) {
    m_prec = p;
    updateDerivedValues();
  }
}

void NC::GaussMos::setDSpacingSpread(double dd)
{
  if (dd==m_delta_d)
    return;
  NCRYSTAL_THROW(LogicError,"GaussMos::setDSpacingSpread not actually implemented and debugged fully yet");//TODO: implement+test this?
  nc_assert(dd>=0&&dd<0.99);
  m_delta_d = dd;
}

double NC::GaussMos::calcRawCrossSectionValueInit(InteractionPars& ip, double cos_angle_indir_normal ) const
{
  nc_assert(ip.m_Q<=0.);
  if ( ip.m_Qprime == -1 ) {
    //First call after getting new wavelength or d-spacing.
    ip.m_cos_perfect_theta = std::sqrt(ip.m_cos_perfect_theta_sq);// = sinalpha
    double tmp2 = ip.m_cos_perfect_theta*ip.m_sin_perfect_theta;
    if (tmp2>0) {
      //usual case:
      ip.m_Qprime = ip.m_wl3 / tmp2;
    } else {
      //Avoid a zero division when wavelength is 0 or 2*dspacing. The actual
      //cross-section will be 0 or infinity at the two limits (special value -2.0
      //designates infinity).
      ip.m_Qprime = (ip.m_sin_perfect_theta>0.5&&ip.m_xsfact) ? -2.0/*wl~=2dsp*/ : 0.0/*wl~=0*/;
    }
  }
  if (ip.m_Qprime>0.) {
    ip.m_Q = ip.m_Qprime * ip.m_xsfact;
    nc_assert(ip.m_Q>0.);
    return calcRawCrossSectionValue(ip,cos_angle_indir_normal);
  }

  if (ip.m_Qprime) {
    nc_assert(ip.m_Qprime==-2);
    return kInfinity;//Q=inf, W factor irrelevant (assuming called within truncation radius)
  } else {
    return 0.0;//Q=0, W factor irrelevant
  }
}

double NC::GaussMos::calcCrossSections( InteractionPars& ip,
                                        const NC::Vector& indir,
                                        const std::vector<NC::Vector>& deminormals,
                                        std::vector<NC::GaussMos::ScatCache>& cache,
                                        VectD& xs_commul ) const
{
  nc_assert(ip.isValid()&&ip.m_wl>0);
  nc_assert(indir.isUnitVector());
  std::vector<Vector>::const_iterator it(deminormals.begin()), itE(deminormals.end());
  double xsoffset = xs_commul.empty() ? 0.0 : xs_commul.back();
  double xssum(0.0);
  const double cptsq = ip.m_cos_perfect_theta_sq;
  const double cta = m_gos.getCosTruncangle();
  for(;it!=itE;++it) {
    const Vector& normal = *it;
    const double dot = normal.dot(indir);
    double sdotcptsq = (1.0 - dot * dot)*cptsq;
    double ds = dot * ip.m_sin_perfect_theta;

    //First a combined check, which usually allows us to skip both normal and
    //anti-normal:
    double A0 = ncmax( 0.0, cta - ncabs(ds) );
    if ( sdotcptsq <= A0*A0 )
      continue;

    //At least one of the two normals should contribute, so deal with them:
    double Am = ncmax( 0.0, cta - ds );
    if ( sdotcptsq > Am*Am ) {
      //anti-normal is within truncated Gauss
      double xs = calcRawCrossSectionValue(ip, dot );
      if (xs) {
        xs_commul.push_back(xsoffset + (xssum += xs));
        cache.emplace_back(-normal, ip.m_inv2dsp);
      }
    }
    double Ap = ncmax( 0.0, cta + ds );
    if ( sdotcptsq > Ap*Ap ) {
      //normal is within truncated Gauss
      double xs = calcRawCrossSectionValue(ip, -dot );
      if (xs) {
        xs_commul.push_back(xsoffset + (xssum += xs));
        cache.emplace_back(normal, ip.m_inv2dsp);
      }
    }

  }
  return xssum;
}

void NC::GaussMos::genScat( RNG& rng, const ScatCache& cache, double wl_raw, const NC::Vector& indir, NC::Vector& outdir) const
{
  nc_assert(wl_raw>0.);
  nc_assert(cache.plane_inv2d()>0.);
  double wl = GaussMos_cacheRound(wl_raw);
  double inv2d = GaussMos_cacheRound(cache.plane_inv2d());
  nc_assert(wl>0);
  nc_assert(inv2d>0);

  //NB: indir has opposite sign convention from that used in derivation of equations. Hence a few extra minus signs in this method.
  nc_assert(cache.isValid());
  nc_assert(indir.isUnitVector());
  double sinthetabragg = wl * inv2d;
  nc_assert(ncabs(sinthetabragg)<=1.0);

  if ( sinthetabragg==0. ) {
    //completely forward scattering
    outdir = indir;
    return;
  }

  //NB: In principle we could cache sa and sg on the cache object, for increased
  //efficiency when scattering from the same neutron state constantly. However,
  //this is a small improvement, coming with mem/cpu costs, and would not help
  //the usual stepping-MC case.
  double ca = sinthetabragg;
  double sa = std::sqrt(1.0-ca*ca);
  double cg = ncclamp(-(indir.dot(cache.plane_normal())),-1.0,1.0);
  double sg = std::sqrt(1.0-cg*cg);
  double ct,st;
  if (!m_gos.genPointOnCircle( rng, cg,sg,ca,sa,ct,st) ) {
    //something went wrong, or numerical imprecision caused us to be called at
    //vanishing cross-section. Assume the latter.
    outdir = indir;
    return;
  }

  //genPointOnCircle selected an actual position of the normal for us (e.g. it
  //selected an actual crystallite). It did not do that in the lab frame, but in
  //a frame where the neutron sits at (0,0,1) rather than -indir, while the
  //nominal normal sits at (sg,0,cg) rather than cache.plane_normal, and the
  //selected normal sits at (sa*ct,sa*st,ca). Perform the reflection of the
  //neutron in that frame as well (using double-angle formulas to get
  //sin(2*alpha) and cos(2*alpha)):
  double s2a = 2*sa*ca;
  double c2a = ca*ca-sa*sa;
  outdir.set(s2a*ct,s2a*st,c2a);

  //Now rotate back to the lab frame:
  rotateToFrame( sg, cg, cache.plane_normal(), -indir, outdir, &rng );

  //For numerical safety:
  outdir.normalise();
  return;
}

void NC::GaussMos::InteractionPars::set(double wl_raw, double inv2dsp_raw, double xsfact) {
  nc_assert(wl_raw>0);
  nc_assert(inv2dsp_raw>0);
  nc_assert(xsfact>0.);
  m_xsfact = xsfact * 0.5;//absorb factor 1/2 here, to avoid a multiplication later
  double wl = GaussMos_cacheRound(wl_raw);
  double inv2dsp = GaussMos_cacheRound(inv2dsp_raw);
  nc_assert(wl>0);
  nc_assert(inv2dsp>0);
  if (wl==m_wl) {
    if (inv2dsp==m_inv2dsp) {
      //great, we don't have to invalidate the Qprime/alpha values, but m_Q might have changed:
      m_Q = (m_Qprime > 0.0 ? m_Qprime * m_xsfact : -1/*will trigger recalc later*/);
      return;
    }
  } else {
    nc_assert(wl>0);
    m_wl = wl;
    m_wl3 = m_wl*m_wl*m_wl;
  }
  //wl or inv2dsp changed => new thetabragg and invalid Q/Qprime/alpha/cos_perfect_theta values:
  nc_assert(inv2dsp>0.);
  m_inv2dsp = inv2dsp;
  m_sin_perfect_theta = wl * inv2dsp;
  nc_assert(valueInInterval(0.,1.,m_sin_perfect_theta));
  m_cos_perfect_theta_sq = 1 - m_sin_perfect_theta*m_sin_perfect_theta;
  m_Q = m_Qprime = m_cos_perfect_theta = -1;//invalidate
  m_alpha = -99;//invalidate
}

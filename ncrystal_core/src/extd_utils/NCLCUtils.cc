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

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/internal/extd_utils/NCLCUtils.hh"
#include "NCrystal/internal/extd_utils/NCPlaneProvider.hh"
#include "NCrystal/internal/phys_utils/NCGaussMos.hh"
#include "NCrystal/internal/utils/NCRomberg.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include <functional>//std::greater

namespace NC = NCrystal;

//uncomment to generate special data files for debugging the overlay functions:
//#define NCRYSTAL_LCHELPER_WRITE_OVERLAYS
#ifdef NCRYSTAL_LCHELPER_WRITE_OVERLAYS
#  include <set>
#  include <fstream>
#  include <sstream>
#endif

//uncomment for debugging purposes, to exclude anti-normals or include only anti-normals:
//#define NCRYSTAL_LCUTILS_ANTINORMALS_EXCLUDED
//#define NCRYSTAL_LCUTILS_ANTINORMALS_ONLY

#define NCRYSTAL_LCUTILS_DISCRFACT (1099511627776.0) // 2^40 ~= 1.1e12
namespace NCRYSTAL_NAMESPACE {
  uint64_t LCdiscretizeValue(double value) {
    nc_assert_always(value>=0.0&&value<1e7);//range limited by uint64_t bits
    return static_cast<uint64_t>(value*NCRYSTAL_LCUTILS_DISCRFACT+0.5);
  }
  double LCdediscretizeValue(uint64_t dvalue) {
    return dvalue * ( 1.0 / NCRYSTAL_LCUTILS_DISCRFACT );
  }
  typedef std::pair<uint64_t,uint64_t> LCInitKey;//discretised (dspacing,alpha)
  typedef std::map<LCInitKey,LCPlaneSet,std::greater<LCInitKey> > LCInitMap;
}

NC::LCHelper::LCHelper( NC::LCAxis lcaxis,
                        NC::LCAxis lcaxis_labframe,
                        MosaicityFWHM mosaicity_fwhm,
                        double unitcell_volume_times_natoms,
                        PlaneProvider* pp,
                        double prec,
                        double ntrunc )
  : m_lcaxislab(lcaxis_labframe.as<Vector>().unit()),
    m_lcstdframe(mosaicity_fwhm,prec,ntrunc),
    m_xsfact( 1.0 / unitcell_volume_times_natoms )
{
  nc_assert(pp);
  nc_assert_always(unitcell_volume_times_natoms>0);
  nc_assert_always(m_xsfact>0);
  nc_assert_always(lcaxis.as<Vector>().isUnitVector());
  nc_assert_always(lcaxis_labframe.as<Vector>().isUnitVector());
  lcaxis.as<Vector>().normalise();//clear tiny errors

  //Collect planes into temporary map, in order to merge those with similar
  //(angle2lcaxis,dspacing) values:
  LCInitMap initmap;
  pp->prepareLoop();
  {
    Optional<PlaneProvider::Plane> opt_plane;
    while ( (opt_plane = pp->getNextPlane() ).has_value() ) {
      auto& plane = opt_plane.value();

      nc_assert(plane.demi_normal.isUnitVector());
      nc_assert(plane.fsq>0.0);
      nc_assert(plane.dspacing>0.0);

      //The only thing that matters is the angle between lcaxis (in crystal frame)
      //and that of the normals. By using the absolute value, we pick the one of the
      //two normals pointing into the same hemisphere as lcaxis:
      double alpha;
      {
        double cosalpha = ncabs ( lcaxis.as<Vector>().dot(plane.demi_normal) );
        if (ncabs(cosalpha)>0.9999999) {
          //deminormal is parallel to lcaxis!
          alpha = 0;
        } else {
          alpha = std::acos(cosalpha);
        }
      }
      nc_assert( alpha>=0.0 && alpha <= kPiHalf );

      //avoid floating point keys + merge entries withing 1/DISCRFACT ~= 1e-12:
      nc_assert_always(plane.dspacing<1e7);//range limited by uint64_t bits
      uint64_t ui_dsp = LCdiscretizeValue(plane.dspacing);
      uint64_t ui_alpha = LCdiscretizeValue(alpha);

      LCInitKey key(ui_dsp,ui_alpha);
      LCInitMap::iterator it = initmap.find(key);
      if ( it == initmap.end() ) {
        //For consistency, alpha/dspacing values should be derived from the key:
        double dsp_of_key = LCdediscretizeValue(ui_dsp);
        double alpha_of_key = LCdediscretizeValue(ui_alpha);
        //rounding can give alpha_of_key slightly above pi/2, but should not be much!
        nc_assert(alpha_of_key<kPiHalf*(1.0+1e-10));
        nc_assert(alpha_of_key>=0);
        nc_assert(dsp_of_key>0);
        LCPlaneSet ps( dsp_of_key, ncmin(alpha_of_key,kPiHalf), m_lcstdframe.gaussMos().mosaicityTruncationAngle(), plane.fsq );
        initmap.insert(std::make_pair(key,ps));
      } else {
        it->second.addFsq(plane.fsq);
      }
    }
  }
  //Finish up initialisation by transferring planes from the temporary map to m_planes:
  m_planes.reserve(initmap.size());
  LCInitMap::iterator it(initmap.begin()),itE(initmap.end());
  for (;it!=itE;++it)
    m_planes.push_back(it->second);
}

NC::LCPlaneSet::LCPlaneSet(double dspacing, double thealpha,
                           double truncangle, double fsquared)
  : twodsp(2.0*dspacing),
    inv_twodsp(0.5/dspacing),
    cosalpha( thealpha ? cos_mpi2pi2(thealpha) : 1.0 ),
    sinalpha( thealpha ? sin_mpi2pi2(thealpha) : 0.0 ),
    cosalphaminus( thealpha>truncangle ? cos_mpi2pi2(thealpha-truncangle) : 1.0),
    cosalphaplus( cos_mpipi(thealpha+truncangle) ),
    fsq(fsquared)
{
  nc_assert(fsquared>=0.0);
  nc_assert(thealpha<=kPiHalf);
  nc_assert(thealpha>=0.0);
  nc_assert(dspacing>0.0);
  nc_assert(truncangle>0.0);
  nc_assert(truncangle<kPiHalf*0.999999);
  nc_assert(thealpha + truncangle<kPi);
  nc_assert(ncmax(0.0,thealpha - truncangle)<kPiHalf);
  nc_assert(cosalphaplus<cosalphaminus);
}

double NC::LCHelper::braggThreshold() const
{
  return m_planes.empty() ? 0.0 : m_planes.begin()->twodsp;
}

NC::LCROIFinder::LCROIFinder(double wl, double c3, double cta, double sta)
  : m_wl(wl),
    m_c3(ncabs(c3)),//ncabs, to ensure alpha_neutron < pi/2 (rotation symmetry guarantees same results)
    m_s3(std::sqrt(ncabs(1.0-c3*c3))),//sinus always positive on [0,pi]
    m_cta(cta),
    m_sta(sta),
    m_prev2d(-1.),
    m_c2(-99),
    m_c23(-99),
    m_s2(-99),
    m_s23(-99),
    m_s2approximated(true)
{
  //The following notation is used for angles in this class:
  //
  //alpha1: angle between demi-normal and lcaxis [provided on the planesets]  [NB: This is denoted theta_n in the paper!]
  //alpha2: pi/2-theta_Bragg where sin(theta_Bragg) = lambda/2*dspacing.      [NB: This is denoted alpha in the paper!]
  //alpha3: angle between neutron and lcaxis.                                 [NB: This is denoted theta_k in the paper!]
  //
  //We also shorten cos(alphai)=ci and sin(alphai)=si, e.g. c2 is cos(alpha2)
  nc_assert(m_c3<=1.0&&m_c3>=0);
  nc_assert(wl>0.0);
  nc_assert(ncabs(cta*cta+sta*sta-1.0)<1e-6);
}

void NC::LCROIFinder::findROIs(const LCPlaneSet * plane,std::vector<NC::LCROI>& roilist)
{
  const double c1minus = plane->cosalphaminus;
  const double c1plus = plane->cosalphaplus;

  nc_assert( m_wl <= plane->twodsp );
  if ( m_prev2d != plane->twodsp ) {
    nc_assert(!ncisnan(m_prev2d)&&!ncisnan(plane->twodsp));
    m_prev2d = plane->twodsp;
    m_c2 = plane->inv_twodsp * m_wl;// cos(alpha2)=cos(pi/2-theta_bragg)=sin(theta_bragg)=wl/2*d
    nc_assert( m_c2>0 && m_c2<=1.0 );
    m_c23 = m_c3*m_c2;
    //m_s2=sqrt(1-c2^2), but we delay the sqrt call since we can reject a lot of
    //planes using just an approximation (Taylor expansion around c2=0,
    //neglecting O(c2^6) terms). The approximated value is strictly larger than
    //sqrt(1-c2^2), which is crucial to the implementation here.
    m_s2approximated = true;
    double c2sq = m_c2*m_c2;
    m_s2 = 1.0-(0.5-0.125*c2sq)*c2sq;//approximation
    m_s23 = m_s3*m_s2;
    nc_assert(m_s23>=0);
  }

  //We first check whether the z-range of the alpha2 cone around the neutron
  //overlaps with the mosaicity band around the normals in the planeset. If it
  //does not, we are done with both normal and anti-normal already, since it is
  //not possible to overlap with the anti-normal's mosaicity band unless there
  //is an overlap with the normal's mosaicity band.

  //NB: We are using the identity cos(alpha2-+alpha3) = m_c23+-m_s23.

  if (intervalsDisjoint(m_c23-m_s23,m_c23+m_s23,c1plus,c1minus)) {
#ifndef NDEBUG
    //Using exact s2 value, verify that neither normal or anti-normal can
    //contribute:
    if (m_s2approximated) {
      const double exacts23 = m_s3*std::sqrt(ncabs(1.0-m_c2*m_c2));
      nc_assert(exacts23>=0);
      nc_assert(intervalsDisjoint(m_c23-exacts23,m_c23+exacts23,c1plus,c1minus));
      nc_assert(intervalsDisjoint(m_c23-exacts23,m_c23+exacts23,-c1minus,-c1plus));
    } else {
      //sanity check that anti-normal does not contribute
      nc_assert(intervalsDisjoint(m_c23-m_s23,m_c23+m_s23,-c1minus,-c1plus));
    }
#endif
    return;
  }

  //It would appear that the normals in the plane-set contribute to the
  //cross-section. However, if the check was made using the approximated value
  //of m_s2, we have to redo it using the precise value, since it could be that
  //it did not actually contribute:
  if (m_s2approximated) {
    m_s2approximated = false;
    m_s2 = std::sqrt( 1.0 - m_c2 * m_c2 );
    m_s23 = m_s3*m_s2;
    nc_assert(m_s23>=0);
    if (intervalsDisjoint(m_c23-m_s23,m_c23+m_s23,c1plus,c1minus)) {
      //sanity check that anti-normal does not contribute
      nc_assert(intervalsDisjoint(m_c23-m_s23,m_c23+m_s23,-c1minus,-c1plus));
      return;
    }
  }

  //Ok, the mosaicity band around the normals does overlap with the alpha2 cone
  //around the neutron. Check if this is also the case for the anti-normals:

  const bool antinormal_contributes = intervalsOverlap(m_c23-m_s23,m_c23+m_s23,-c1minus,-c1plus);

  //Time to add ROI's. First handle the two degenerate cases:
  if (plane->isOnAxis()) {
    //normals are parallel to lcaxis, add special on-axis ROI(s):
    roilist.emplace_back(plane,1.0/*normal*/);
    if (antinormal_contributes)
      roilist.emplace_back(plane,-1.0/*anti-normal*/);
    return;
  }

  if (ncabs(m_s3)<1e-10) {
    //Neutron is aligned with lcaxis, all crystallite rotations give the same
    //result. Add special ROI(s) indicating this.
    roilist.emplace_back(plane,1.0/*normal*/);
    if (antinormal_contributes)
      roilist.emplace_back(plane,-1.0/*anti-normal*/);
    return;
  }

  //Usual non-degenerate case, must determine the appropriate phi-range of
  //crystallite rotations which leads to contributions:

  const double c2ta   = m_c2 * m_cta;//>0
  const double s2ta   = m_s2 * m_sta;//>0

  const double c2low  = c2ta - s2ta;
  const double c2high = ( m_s2<m_sta ? 1.0 : c2ta+s2ta );
  nc_assert(ncabs(c2high-std::cos(ncmax(0.0,std::acos(m_c2)-std::acos(m_cta))))<1e-10);

  const double c1 = plane->cosalpha;
  const double s1 = plane->sinalpha;
  const double a = c1 * m_c3;//cosalpha*normsign*m_c3 [NB: is this an obsolete comment?]
  const double invb = s1 * m_s3;
  nc_assert(invb>0.0);//s1~0,m_s3~=0 degenerate cases handled above
  const double b = 1.0 / invb;
  double mab = -a*b;
  double c2lowb = c2low*b;
  double c2highb = c2high*b;
  //NB: Can in principle avoid expensive acos calls when cos values are +-1. But
  //might be good to approximate acos anyway?
  double cosphi1 = ncmax(-1.0,ncmin(1.0, mab + c2lowb ));
  double cosphi2 = ncmax(-1.0,ncmin(1.0, mab + c2highb ));;
  const double mindist = 1e-10;//do not add ROIs with too short phi-range.
  if ( ncabs(cosphi1-cosphi2) > mindist )
    roilist.emplace_back(std::acos(ncmax(cosphi1,cosphi2)),std::acos(ncmin(cosphi1,cosphi2)),plane,1.0/*normal*/);
  if (antinormal_contributes) {
    cosphi1 = ncmax(-1.0,ncmin(1.0, mab - c2lowb ));
    cosphi2 = ncmax(-1.0,ncmin(1.0, mab - c2highb ));;
    if ( ncabs(cosphi1-cosphi2) > mindist )
      roilist.emplace_back(std::acos(ncmax(cosphi1,cosphi2)),std::acos(ncmin(cosphi1,cosphi2)),plane,-1.0/*anti-normal*/);
  }
}

double NC::LCHelper::crossSectionNoCache( double wl, const Vector& indir ) const
{
  Cache cache;
  return crossSection( cache, wl, indir );
}

void NC::LCHelper::genScatterNoCache( RNG& rng, double wl, const Vector& indir, Vector& outdir ) const
{
  Cache cache;
  genScatter( cache, rng, wl, indir, outdir );
}

bool NC::LCHelper::isValid(NC::LCHelper::Cache& cache, double wl, double c3 ) const
{
  //ncabs, to ensure alpha_neutron <= pi/2 (rotation symmetry guarantees same results)
  nc_assert(wl>=0&&wl<1e7&&c3>=-1.0&&c3<=1.0);
  return cache.m_signature.first == LCdiscretizeValue(wl)
    && cache.m_signature.second == LCdiscretizeValue(ncabs(c3));
}

bool NC::LCHelper::isValid(NC::LCHelper::Cache& cache, double wl, const NC::Vector& indir) const
{
  nc_assert(wl>0.0);
  nc_assert(indir.isUnitVector());
  return isValid(cache,wl,m_lcaxislab.dot(indir));
}

void NC::LCHelper::ensureValid(NC::LCHelper::Cache& cache, double wl, const NC::Vector& indir) const
{
  nc_assert(wl>0.0);
  nc_assert(indir.isUnitVector());
  double c3 = m_lcaxislab.dot(indir);
  nc_assert(wl>=0&&wl<1e7&&c3>=-1.0&&c3<=1.0);
  uint64_t discrwl = LCdiscretizeValue(wl);
  uint64_t discrc3 = LCdiscretizeValue(ncabs(c3));
  if ( cache.m_signature.first == discrwl && cache.m_signature.second == discrc3 )
    return;
  forceUpdateCache(cache,discrwl,discrc3);
}


void NC::LCHelper::forceUpdateCache( NC::LCHelper::Cache& cache, uint64_t discr_wl, uint64_t discr_c3 ) const
{
  cache.m_signature.first = discr_wl;
  cache.m_signature.second = discr_c3;
  cache.m_wl = LCdediscretizeValue(discr_wl);
  cache.m_c3 = LCdediscretizeValue(discr_c3);
  nc_assert(cache.m_c3>=0&&cache.m_c3<=1.0+1e-6);
  cache.m_c3 = ncmin(cache.m_c3,1.0);
  cache.m_s3 = std::sqrt(ncabs(1.0-cache.m_c3*cache.m_c3));
  cache.m_roilist.clear();
  cache.m_roixs_commul.clear();
  cache.m_roi_overlays.clear();
  const double wl(cache.m_wl), c3(cache.m_c3), s3(cache.m_s3);

  const double cta = m_lcstdframe.gaussMos().mosaicityCosTruncationAngle();
  const double sta = m_lcstdframe.gaussMos().mosaicitySinTruncationAngle();
  LCROIFinder roifinder(wl,c3,cta,sta);
  std::vector<LCPlaneSet>::const_iterator it(m_planes.begin()),itE(m_planes.end());
  for (;it!=itE;++it) {
    if ( wl > it->twodsp )
      break;//done, no other planes can contribute, since m_planes is sorted by dspacing
#ifndef NDEBUG
    std::size_t nold = cache.m_roilist.size();
#endif
    roifinder.findROIs(&(*it),cache.m_roilist);
#ifndef NDEBUG
    {
      //sanity check of ROI ranges:
      Vector vneutron(s3,0.,c3);
      for (std::size_t ii = nold; ii<cache.m_roilist.size(); ++ii) {
        LCROI & roi = cache.m_roilist.at(ii);
        if (roi.isDegenerate())
          continue;
        LCStdFrame::NormalPars normal(roi.planeset,roi.normal_sign);
        Vector v1 = LCStdFrame::normalInStdFrame(normal,std::cos(roi.rotmin),std::sin(roi.rotmin));
        Vector v2 = LCStdFrame::normalInStdFrame(normal,std::cos(roi.rotmax),std::sin(roi.rotmax));
        const double idealangle = std::acos( wl * roi.planeset->inv_twodsp );
        nc_assert(ncabs(vneutron.angle(v1)-idealangle)<m_lcstdframe.gaussMos().mosaicityTruncationAngle()*1.0001);
        nc_assert(ncabs(vneutron.angle(v2)-idealangle)<m_lcstdframe.gaussMos().mosaicityTruncationAngle()*1.0001);
      }
    }
#endif
#if defined(NCRYSTAL_LCUTILS_ANTINORMALS_ONLY) || defined(NCRYSTAL_LCUTILS_ANTINORMALS_EXCLUDED)
    {
      decltype(cache.m_roilist) modified_roilist;
      //    double normal_sign;//1.0 for normal(s), -1.0 for anti-normal(s).
#  ifdef NCRYSTAL_LCUTILS_ANTINORMALS_ONLY
      constexpr double modify_target_normsign = -1.0;
#  else
      constexpr double modify_target_normsign = 1.0;
#  endif
      for (const auto& e: cache.m_roilist) {
        if (e.normal_sign == modify_target_normsign)
          modified_roilist.emplace_back(e);
      }
      cache.m_roilist.swap(modified_roilist);
    }
#endif
  }

  if (cache.m_roilist.empty())
    return;

  cache.m_roixs_commul.reserve(cache.m_roilist.size());
  double sumxs = 0.0;
  std::vector<LCROI>::const_iterator itROI(cache.m_roilist.begin()),itROIE(cache.m_roilist.end());

  LCStdFrame::NeutronPars neutron(wl,c3,s3);

  for ( ; itROI!=itROIE; ++itROI ) {
    LCStdFrame::NormalPars normal(itROI->planeset,itROI->normal_sign);
    double roi_xs;
    if (itROI->isDegenerate()) {
      if (itROI->normalIsOnAxis()) {
        roi_xs = m_lcstdframe.calcXS_OnAxis(neutron,normal);
      } else {
        nc_assert(itROI->neutronIsOnAxis());
        //xs is independent of particular crystallite rotation around lcaxis. We
        //evaluate at phi=90deg, as an average value in case of infinitesimal
        //variations when is only aaalmost on axis..
        roi_xs = m_lcstdframe.calcXS(neutron,normal,0.0/*cos(phi=90deg)*/);
      }
    } else {
      //Result of integration needs division by interval length (pi) to get the
      //average, which we calculated over [0,pi] rather than [-pi,pi]. Due to
      //the symmetry xs(phi)=xs(-phi), this gives the same result:
      roi_xs = m_lcstdframe.calcXSIntegral(neutron,normal,itROI->rotmin,itROI->rotmax) * kInvPi;
    }
    nc_assert(roi_xs>0);//otherwise it should not have been a ROI!
    cache.m_roixs_commul.push_back(sumxs += roi_xs);
  }

  nc_assert(cache.m_roixs_commul.size()==cache.m_roilist.size());
}

double NC::LCHelper::crossSection( NC::LCHelper::Cache& cache, double wl, const NC::Vector& indir ) const
{
  ensureValid(cache,wl,indir);
  return cache.m_roixs_commul.empty() ? 0.0 : (m_xsfact * cache.m_roixs_commul.back());
}

void NC::LCHelper::Cache::reset()
{
  //same result as Cache() constructor
  m_signature.first = m_signature.second = std::numeric_limits<uint64_t>::max();
  m_wl = m_c3 = m_s3 = -99.0;
  m_roilist.clear();
  m_roixs_commul.clear();
}

namespace NCRYSTAL_NAMESPACE {
  class LCStdFrameIntegrator : public Romberg {
  public:
    LCStdFrameIntegrator(const GaussMos* gm, const LCStdFrame::NormalPars& normal, const LCStdFrame::NeutronPars& neutron)
      : Romberg(),
        m_ip(neutron.wl, normal.planeset->inv_twodsp, normal.planeset->fsq),
        m_gm(gm),
        m_sinnormalmults3(normal.planeset->sinalpha*neutron.s3*normal.sign),//include normal.sign here, since it should multiply cos(phi) in evalFunc
        m_cosnormalmultc3(normal.planeset->cosalpha*normal.sign*neutron.c3),
        m_acc(ncclamp(gm->precision(),1e-7,1e-2))
    {
      nc_assert(normal.planeset&&!normal.planeset->isOnAxis());
      nc_assert(gm&&neutron.wl>0);
    }

    virtual bool accept(unsigned /*level*/, double prev_estimate, double estimate,double,double) const
    {
      return ncabs(estimate-prev_estimate) <= m_acc*ncabs(estimate);
    }

    virtual double evalFunc(double) const {
      //Must be implemented, but should never be called since we provide
      //evalFuncMany+evalFuncManySum.
      nc_assert_always(false);
    }

    virtual void evalFuncMany(double* fvals, unsigned n, double offset, double delta) const
    {
      nc_assert(offset>=0&&offset<kPi*1.00001);
      nc_assert(delta>0&&delta*(n-1)<=kPi*1.00001);
      CosSinGridGen grid(n,offset,delta);
      unsigned i = 0;
      do {
        double cosgamma = m_sinnormalmults3 * grid.current_cosval() + m_cosnormalmultc3;
        nc_assert(NC::ncabs(cosgamma)<1.000000001);
        fvals[i++] = m_gm->calcRawCrossSectionValue(m_ip,cosgamma);
      } while (grid.step());
    }

    virtual double evalFuncManySum(unsigned n, double offset, double delta) const
    {
      nc_assert(offset>=0&&offset<kPi*1.00001);
      nc_assert(delta>0&&delta*n<=kPi*1.00001);
      CosSinGridGen grid(n,offset,delta);
      double sum(0.);
      do {
        double cosgamma = m_sinnormalmults3 * grid.current_cosval() + m_cosnormalmultc3;
        nc_assert(NC::ncabs(cosgamma)<1.000000001);
        sum += m_gm->calcRawCrossSectionValue(m_ip,cosgamma);
      } while (grid.step());
      return sum;
    }

  private:
    mutable GaussMos::InteractionPars m_ip;//NB: Mutable here is OK, since we use local per-thread LCStdFrameIntegrator objects.
    const GaussMos * m_gm;
    const double m_sinnormalmults3;
    const double m_cosnormalmultc3;
    const double m_acc;
  };
}

void NC::LCHelper::genPhiVal(RNG& rng, const LCROI& roi, const Overlay& overlay, double& phi, double& overlay_at_phi)
{
  const float* it = std::lower_bound( overlay.data, overlay.data+Overlay::ndata, overlay.data[Overlay::ndata-1] * rng.generate() );
  unsigned ichoice = std::min<unsigned>((unsigned)(it - overlay.data),Overlay::ndata-1);
  overlay_at_phi = overlay.nonCommulVal(ichoice);
  double rel_phi_pos = (ichoice + rng.generate())/Overlay::ndata;
  phi = roi.rotmin + rel_phi_pos*roi.length();
}

void NC::LCHelper::genScatter( LCHelper::Cache& cache, RNG& rng, double wl, const Vector& indir, Vector& outdir ) const
{

  ensureValid(cache,wl,indir);

  double roixssum = cache.m_roixs_commul.empty() ? 0.0 : cache.m_roixs_commul.back();
  if (!roixssum) {
    //scattering not possible here.
    outdir = indir;
    return;
  }

  //Choose ROI, according to cross-section of each ROI:
  std::size_t idx = pickRandIdxByWeight(rng,cache.m_roixs_commul);
  nc_assert(idx<cache.m_roilist.size());
  const LCROI& roi = cache.m_roilist[idx];

  //Now, generate the scattering in the chosen ROI. In case of on-Axis ROI, this
  //can go straight ahead. In case of off-Axis, one must first decide upon the
  //exact rotation of the crystallite in which the scattering occurs.


  LCStdFrame::NeutronPars neutron(cache.m_wl,cache.m_c3,cache.m_s3);
  LCStdFrame::NormalPars normal(roi.planeset,roi.normal_sign);

  if ( roi.normalIsOnAxis() ) {
    //Normal parallel to lcaxis.
    m_lcstdframe.genScat_OnAxis(rng,neutron,normal,outdir);
  } else {
    //Off-axis normal. Use rejection-method to pick phi-value in ROI (unless
    //neutron is on axis, in which case we just pick one at random):
    double phi,cosphi;
    if (roi.neutronIsOnAxis()) {
      phi = rng.generate()*kPi;
      cosphi = cos_mpipi(phi);
    } else {

      //Find overlay object:
      if (cache.m_roi_overlays.empty())
        cache.m_roi_overlays.resize(cache.m_roilist.size());
      nc_assert(idx<cache.m_roi_overlays.size());
      Overlay& overlay = cache.m_roi_overlays[idx];

      if (!overlay.data) {

        //Didn't scatter on this normal before, prepare overlay by sampling xs
        //values at edges of overlay histogram bins (For convenience and
        //consistency, use the integrator class to do this):
        double tmp[Overlay::ndata+1];
        LCStdFrameIntegrator integrator(&m_lcstdframe.gaussMos(), normal,neutron);
        integrator.evalFuncMany(&tmp[0], Overlay::ndata+1, roi.rotmin, roi.length()/Overlay::ndata);

        //Adding 2% of maxval to all bins significantly increases safety
        //of non-central bins, with low impact on the acceptance rate:
        double * it(&tmp[0]);
        double * itLast(it+Overlay::ndata);
        double * itE(itLast+1);
        double maxval = 0.0;
        for (;it!=itE;++it)
          maxval = ncmax(maxval,*it);
        double safety_offset = 0.02 * maxval;

        //And a multiplicative factor ensures that the overlay function will
        //never be too small in central bins:
        const double safety_factor = 1.7;

        //Finally, put into overlay.data as commulative array:
        overlay.prepareNullArray();
        float * itData = overlay.data;
        float sum(0.0);
        for (it = &tmp[0]; it!=itLast ; ++it, ++itData )
          *itData = ( sum += (ncmax(*it,*(it+1)) * safety_factor+safety_offset) );

      }
      const int maxtries = 1000;
      int triesleft = maxtries;
#ifdef NCRYSTAL_LCHELPER_WRITE_OVERLAYS
      std::string written_name;
      {
        static std::set<std::tuple<const LCPlaneSet *,double,double,double> > written;
        auto key = std::make_tuple(roi.planeset,cache.m_wl,roi.rotmin,roi.rotmax);
        if (!written.count(key)) {
          written.insert(key);
          std::stringstream s; s<<"phi_overlay_roi"<<written.size()<<".txt";
          std::ofstream ofs (s.str(), std::ofstream::out);
          written_name = s.str();
          std::vector<PairDD > sampleoverlay;
          unsigned n=1000;
          sampleoverlay.reserve(n);
          for (unsigned i=0; i<n; ++i) {
            double ph,overlay_at_phi;
            genPhiVal(rng,roi,overlay,ph,overlay_at_phi);
            sampleoverlay.push_back(std::make_pair(ph,overlay_at_phi));
          }
          std::sort(sampleoverlay.begin(),sampleoverlay.end());
          double ddd = (roi.rotmax-roi.rotmin)/(n-1);
          for (size_t i=0;i<n;++i) {
            double ph = roi.rotmin + i*ddd;
            double overlay_relphi= i*(1.0/(n-1));//maps [rotmin,rotmax] to [0,1]
            nc_assert(overlay_relphi>=0.&&overlay_relphi<=1.+1e-13);
            unsigned overlay_bin = std::min<unsigned>(Overlay::ndata-1,(unsigned)(overlay_relphi*Overlay::ndata));
            ofs << sampleoverlay.at(i).first << " " << sampleoverlay.at(i).second << " "
                << ph << " " << m_lcstdframe.calcXS(neutron,normal,std::cos(ph)) << " "
                << overlay.nonCommulVal(overlay_bin) << "\n";
          }
        }
      }
#endif
      while(triesleft--) {
        double overlay_at_phi;
        genPhiVal(rng,roi,overlay,phi,overlay_at_phi);
        cosphi = cos_mpipi(phi);
        double xsphi = m_lcstdframe.calcXS(neutron,normal,cosphi);
        if ( xsphi > overlay_at_phi ) {
          static bool first = true;
          if (first) {
            first = false;
            std::ostringstream ss;
            ss<<"Problems sampling with rejection method during LCHelper::genScatter "
              "invocation. Overlay function was not larger than actual cross-section value at sampled point "
              "(overshot by factor of "<<(overlay_at_phi?xsphi/overlay_at_phi:kInfinity)<<"). Further warnings"
              " of this type will not be emitted.";
#ifdef NCRYSTAL_LCHELPER_WRITE_OVERLAYS
            ss<<" The associated overlay file was written to: "<<written_name;
            Msg::outputMsg(ss.str(),MsgType::Warning);
#endif
          }
        }
        if ( xsphi > overlay_at_phi * rng.generate() )
          break;
      }
      if (triesleft<=0) {
        static bool first = true;
        if (first) {
          first = false;
          std::ostringstream ss;
          ss<<"NCrystal WARNING: Problems sampling with rejection method during LCHelper::genScatter "
            "invocation. Did not accept sampled value after "<<maxtries<<" attempts. Further warnings"
            " of this type will not be emitted.";
#ifdef NCRYSTAL_LCHELPER_WRITE_OVERLAYS
          ss<<" The associated overlay file was written to: "<<written_name;
          Msg::outputMsg(ss.str(),MsgType::Warning);
#endif
        }
      }
    }
    nc_assert(phi>=0&&phi<=kPi);
    double sinphisign = (rng.coinflip()?1.0:-1.0);//pick normals in [-pi,pi], not just in [0,pi]
    m_lcstdframe.genScat(rng,neutron,normal,cosphi,sinphisign*std::sqrt(1.0-cosphi*cosphi),outdir);
  }

  //The scattering captured in outdir took place in the standard frame (see
  //description of LCStdFrame). Here, the lcaxis was aligned with the z-axis
  //while indir=(s3,0,c3). We rotate outdir from this frame to the lab frame:
  double indirsign = m_lcaxislab.dot(indir) >= 0.0 ? 1.0 : -1.0;
  rotateToFrame( cache.m_s3, cache.m_c3, indir, m_lcaxislab*indirsign, outdir, &rng );

  //Throughout LCUtils internals we have for (misguided attempts at?) clarity,
  //employed a sign convention for the neutron direction vector which is
  //opposite to its momentum. Thus, we must perform a final sign flip on the
  //resulting direction:
  outdir *= -1.0;
}

NC::LCStdFrame::LCStdFrame(MosaicityFWHM mosaicity, double prec, double ntrunc)
  : m_gm(mosaicity,prec,ntrunc)
{
}

double NC::LCStdFrame::calcXS_OnAxis( const NC::LCStdFrame::NeutronPars& neutron,
                                      const NC::LCStdFrame::NormalPars& normal ) const
{
  nc_assert(normal.planeset->isOnAxis());
  double cosval = normal.sign * neutron.c3;
  GaussMos::InteractionPars ip(neutron.wl, normal.planeset->inv_twodsp, normal.planeset->fsq);
  return m_gm.calcRawCrossSectionValue( ip, cosval );
}

void NC::LCStdFrame::genScat_OnAxis( RNG& rng,
                                     const NC::LCStdFrame::NeutronPars& neutron,
                                     const NC::LCStdFrame::NormalPars& normal,
                                     NC::Vector& outdir ) const
{
  GaussMos::ScatCache scatcache(Vector(0.,0.,normal.sign), normal.planeset->inv_twodsp );

  //Invoke GaussMos helper to carry out the actual  scattering. Note that the
  //sign convention of the neutron indir vector in the GaussMos class differs
  //from the one used throughout LCUtils, hence the minus signs in the
  //definition of indir_stdframe:
  Vector indir_stdframe(-neutron.s3,0.,-neutron.c3);
  m_gm.genScat( rng, scatcache, neutron.wl, indir_stdframe, outdir );
}

NC::Vector NC::LCStdFrame::normalInStdFrame( const NC::LCStdFrame::NormalPars& normal, double cosphi, double sinphi )
{
  nc_assert( ncabs ( cosphi*cosphi+sinphi*sinphi - 1.0 ) < 1.0e-10 );
  const double c1 = normal.planeset->cosalpha;
  const double s1 = normal.planeset->sinalpha;
  const double ns1 = normal.sign*s1;
  return Vector(ns1*cosphi, ns1*sinphi,normal.sign*c1);
}

double NC::LCStdFrame::calcXS( const NC::LCStdFrame::NeutronPars& neutron,
                               const NC::LCStdFrame::NormalPars& normal,
                               double cosphi ) const
{
  nc_assert(!normal.planeset->isOnAxis());
  nc_assert(ncabs(cosphi)<=1);
  //gamma is angle between neutron and plane normal:
  const double c1 = normal.planeset->cosalpha;
  const double s1 = normal.planeset->sinalpha;
  const double cosgamma = (s1 * neutron.s3 * cosphi + c1 * neutron.c3)*normal.sign;//normal.sign multiplies both cosphi and c1
  GaussMos::InteractionPars ip(neutron.wl, normal.planeset->inv_twodsp, normal.planeset->fsq);
  return m_gm.calcRawCrossSectionValue( ip, cosgamma );
}


double NC::LCStdFrame::calcXSIntegral( const NC::LCStdFrame::NeutronPars& neutron,
                                       const NC::LCStdFrame::NormalPars& normal,
                                       double phimin, double phimax ) const
{
  nc_assert(phimin>=0.0&&phimin<kPi);
  nc_assert(phimax>0.0&&phimax<=kPi);
  nc_assert(phimax>phimin);
  LCStdFrameIntegrator integrator(&m_gm, normal,neutron);
  return integrator.integrate(phimin,phimax);
}

void NC::LCStdFrame::genScat( RNG& rng,
                              const NC::LCStdFrame::NeutronPars& neutron,
                              const NC::LCStdFrame::NormalPars& normal,
                              double cosphi,
                              double sinphi,
                              NC::Vector& outdir ) const
{
  nc_assert(ncabs(cosphi*cosphi+sinphi*sinphi-1.0)<=1e-10);

  GaussMos::ScatCache scatcache(normalInStdFrame( normal, cosphi, sinphi ), normal.planeset->inv_twodsp);

  //Invoke GaussMos helper to carry out the actually scattering. Note that the
  //sign convention of the neutron indir vector in the GaussMos class differs
  //from the one used throughout LCUtils, hence the minus signs in the
  //definition of indir_stdframe:
  Vector indir_stdframe(-neutron.s3,0.,-neutron.c3);
  m_gm.genScat( rng, scatcache, neutron.wl, indir_stdframe, outdir );
}

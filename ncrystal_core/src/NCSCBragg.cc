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

#include "NCrystal/NCSCBragg.hh"
#include "NCrystal/NCSCOrientation.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCException.hh"
#include "NCMath.hh"
#include "NCEqRefl.hh"
#include "NCSCGeoComputation.hh"
#include <algorithm>
#include <vector>
#include <iomanip>

#define NCSCBragg_LACKSORIENTATION (-2.0)
#define NCSCBragg_INVALIDATECACHE (-1.0)
#define NCSCBragg_NSIGMA (5)
#define NCSCBragg_CONST1 (0.424660900144009534) // = 1 / (2sqrt(2log(2)))
#define NCSCBragg_CONST2 (0.939437278699651435) // = 1 / ( sqrt(2pi)*CONST1 )

struct NCrystal::SCBragg::pimpl {
  pimpl(SCBragg*scb,const NCrystal::Info* cinfo, double mosaicity, double dd);
  ~pimpl();

  void setOrientationVectors(const Vector&,const Vector&,const Vector&,const Vector&);

  void genScat( double wl, const Vector& indir, Vector& outdir, double& de ) const;
  bool channelXSInLab( unsigned, double,const Vector &,double & ) const ;
  void sampleConeAndReflect(double, double,double,const Vector&,
                            const Vector&, Vector& ) const;
  void updateCache(double wl, const Vector& ) const;
  double calcWavelengthThreshold(double contributionthreshold_times_V0_times_numAtom,
                                 double dspacing, double fsquared);

  struct Cache {
    double wl;
    Vector dir;
    std::vector<double> contribs;
    struct ContribInfo {
      ContribInfo(unsigned f,unsigned c, double n)
        : famidx(f), chanidx(c), normsign(n) {}
      unsigned famidx;
      unsigned chanidx;
      double normsign;
    };
    std::vector<ContribInfo> famchan;
  };

  SCBragg * m_scb;
  const Info * m_cinfo;
  Reflections * m_Reflections;
  SCGeoComputation * m_geo;
  double m_delta_d;     //d-spacing deviation
  double m_gauss_sigma; //mosaic spread
  double m_sqrt2_div_sigma;
  double m_neginv_2_gauss_sigma_sq;
  double m_gauss_norm;
  double m_max_distance_mag2;
  double m_threshold_ekin;
  double m_threshold_wl;
  bool m_enable_backscattering;

  //Cache - should be in tread-local storage if calling this in a
  //multi-threaded application:
  mutable Cache m_cache;
};


NCrystal::SCBragg::pimpl::pimpl(SCBragg*scb,const NCrystal::Info* cinfo, double mosaicity, double dd)
  : m_scb(scb),
    m_cinfo(0),
    m_delta_d(dd),
    m_gauss_sigma( NCSCBragg_CONST1 * mosaicity ),//width of truncated gaussian
    m_gauss_norm( NCSCBragg_CONST2 / ( mosaicity * ncerf( M_SQRT1_2 * NCSCBragg_NSIGMA ) ) ),//normalisation for truncated gaussian
    m_threshold_ekin(std::numeric_limits<double>::infinity()),
    m_threshold_wl(0),
    m_enable_backscattering(true)//hardcode to true for now, pending discussions
{
  nc_assert_always(cinfo);

  m_neginv_2_gauss_sigma_sq = -1.0 / (2.0 * m_gauss_sigma * m_gauss_sigma);
  m_sqrt2_div_sigma = sqrt(2.0) / m_gauss_sigma;

  //TODO for NC2 (question for XX): Is the following comment still correct?
  //With the reflection normal and incident directions as unit vectors, the
  //angle between them is a^2=2*cos(angle) therefore, (N+D).mag2()==2*cos(angle)
  //is true the W function is a Gaussian of angle which peaks at angle=0.  we
  //set the cut-off value of the angle to
  //angle_cut=m_gauss_sigma*NCSCBragg_NSIGMA so, the reflection is only
  //considered when (N+D).mag2() is smaller than 2*cos(angle_cut)

  m_max_distance_mag2 = 2.0 - 2.0 * std::cos( NCSCBragg_NSIGMA * m_gauss_sigma );

  //Start with invalid cache:
  m_cache.wl = NCSCBragg_LACKSORIENTATION;

  //Always need both hkl and structure info:
  if (!cinfo->hasHKLInfo())
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks HKL information.");
  if (!cinfo->hasStructureInfo())
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks Structure information.");

  //Threshold doesn't depend on orientation:
  m_threshold_wl = cinfo->nHKL() ? cinfo->hklBegin()->dspacing * 2.0 : 0.0;
  m_threshold_ekin = wl2ekin(m_threshold_wl);

  //Setup based on structure info:
  m_Reflections = new Reflections();
  m_Reflections->numAtom=cinfo->getStructureInfo().n_atoms;
  m_Reflections->V0=cinfo->getStructureInfo().volume;//the unit for V0 is aa^3

  const StructureInfo & si = cinfo->getStructureInfo();
  m_geo = new SCGeoComputation(si.lattice_a,si.lattice_b,si.lattice_c,
                               si.alpha,si.beta,si.gamma);
  m_cinfo = cinfo;
  m_cinfo->ref();
}

NCrystal::SCBragg::pimpl::~pimpl()
{
  if (m_cinfo)
    m_cinfo->unref();
  delete m_Reflections;
  delete m_geo;
}

NCrystal::SCBragg::SCBragg( const NCrystal::Info* cinfo,
                            const SCOrientation& sco,
                            double mosaicity,
                            double dd )
  : Scatter("SCBragg"),
    m_pimpl(new pimpl(this,cinfo,mosaicity,dd))
{
  if (!sco.isComplete())
    NCRYSTAL_THROW(BadInput,"Incomplete SCOrientation object - must set both primary and secondary directions.");

  Vector dirc[2];
  for (size_t i=0; i < 2; ++i) {
    if (sco.m_crystal_is_hkl[i]) {
      dirc[i] = m_pimpl->m_geo->getReciDir(sco.m_crystal[i][0],sco.m_crystal[i][1],sco.m_crystal[i][2]);
    } else {
      dirc[i] = asVect(sco.m_crystal[i]);
    }
  }

  Vector dirl[2];
  dirl[0] = asVect(sco.m_lab[0]);
  dirl[1] = asVect(sco.m_lab[1]);

  if (dirc[0].isParallel(dirc[1],1.0e-6))
    NCRYSTAL_THROW(BadInput,"Chosen SCOrientation directions in the crystal reference frame are too parallel.");

  if (dirl[0].isParallel(dirl[1],1.0e-6))
    NCRYSTAL_THROW(BadInput,"Chosen SCOrientation directions in the laboratory frame are too parallel.");

  const double anglec = dirc[0].angle(dirc[1]);
  const double anglel = dirl[0].angle(dirl[1]);
  if ( ncabs(anglec-anglel)>sco.m_tolerance ) {
    NCRYSTAL_THROW2(BadInput,"Chosen SCOrientation directions in the lab frame are "<<std::setprecision(8)
                    <<anglel*180/M_PI<<" deg apart, while the chosen directions in the crystal frame"
                    " are "<<anglec*180/M_PI<<" deg apart. This is not within the specified"
                    " tolerance of "<<sco.m_tolerance<<" rad. = "<<sco.m_tolerance*180/M_PI<<" deg.");
  }
  //We are within the tolerance, but now ensure exact anglec==anglel by removing
  //components of secondary direction parallel to the primary direction:

  for (size_t i=0; i < 2; ++i) {
    dirc[i].normalise();
    dirl[i].normalise();
  }
  dirc[1] -= dirc[0] * dirc[1].dot(dirc[0]);
  dirl[1] -= dirl[0] * dirl[1].dot(dirl[0]);
  dirc[1].normalise();
  dirl[1].normalise();

  m_pimpl->setOrientationVectors(dirl[0],dirl[1],dirc[0],dirc[1]);
  validate();
}

NCrystal::SCBragg::~SCBragg()
{
  delete m_pimpl;
}

double NCrystal::SCBragg::pimpl::calcWavelengthThreshold(double contributionthreshold_times_V0_times_numAtom,
                                                         double dspacing, double fsquared)
{
  //This function is essentially calculating the threshold wavelength at which
  //Q(wavelength)==contribution_threshold.
  double CR = contributionthreshold_times_V0_times_numAtom;
  nc_assert(fsquared>0);//otherwise return infinity
  double k = CR / fsquared;
  double d2 = dspacing*dspacing;
  double a = k*k+64*d2;
  nc_assert(a>0&&d2>0);
  return sqrt(k*(sqrt(a)-k)/(8*d2));
}


void NCrystal::SCBragg::pimpl::setOrientationVectors( const NCrystal::Vector& vecInLab1,
                                                      const NCrystal::Vector& vecInLab2,
                                                      const NCrystal::Vector& planeNor1,
                                                      const NCrystal::Vector& planeNor2 )
{
  m_cache.wl = NCSCBragg_INVALIDATECACHE;//Invalidate cache and note that setOrientationVectors was called

  //construct equivalent reflection calculator
  if(vecInLab1.isParallel(vecInLab2) )
    NCRYSTAL_THROW(LogicError,"Vectors in the lab reference frame cannot be parallel.");
  if(planeNor1.isParallel(planeNor2))
    NCRYSTAL_THROW(LogicError,"Vectors in the crystal reference frame cannot be parallel.");
  if(ncabs(vecInLab1.angle(vecInLab2) - planeNor1.angle(planeNor2)) > 1e-8)
    NCRYSTAL_THROW(LogicError,"Vector angles in the two reference frames are different.");

  m_geo->calcTransform(vecInLab1.unit(),vecInLab2.unit(),
                       planeNor1.unit(),planeNor2.unit());

  //expand crystal info
  nc_assert_always(m_cinfo->hasHKLInfo()==true);

  m_Reflections->family_list.clear();
  m_Reflections->family_list.reserve(m_cinfo->nHKL());
  HKLList::const_iterator it_hklE=m_cinfo->hklEnd();

  const double contribution_threshold = 1.0e-5;//barn [NB: hardcoded below as well!]
  const double CR = contribution_threshold*m_Reflections->V0 * m_Reflections->numAtom;

  if(!m_cinfo->hasHKLDemiNormals())
  {
    if (!m_cinfo->getStructureInfo().spacegroup)
      NCRYSTAL_THROW(MissingInfo,"Passed Info object has neither HKL normals nor spacegroup number available.");
    EqRefl equReflCalculator(m_cinfo->getStructureInfo().spacegroup);
    for(HKLList::const_iterator it_hkl=m_cinfo->hklBegin();it_hkl!=it_hklE;++it_hkl)
    {
      const std::set<EqRefl::HKL>& eqPlanes = equReflCalculator.getEquivalentReflections(it_hkl->h, it_hkl->k, it_hkl->l);
      double wlthr = calcWavelengthThreshold(CR,it_hkl->dspacing,it_hkl->fsquared);
#if __cplusplus >= 201103L
      m_Reflections->family_list.emplace_back(it_hkl->h,it_hkl->k,it_hkl->l,it_hkl->fsquared,it_hkl->dspacing,wlthr);
#else
      m_Reflections->family_list.push_back(ReflectionFamily(it_hkl->h,it_hkl->k,it_hkl->l,it_hkl->fsquared,it_hkl->dspacing,wlthr));
#endif
      ReflectionFamily& aFamily = m_Reflections->family_list.back();
      aFamily.eqhkl_normals.reserve(eqPlanes.size());
      for (std::set<EqRefl::HKL>::const_iterator it = eqPlanes.begin(); it != eqPlanes.end(); ++it)
      {
        Vector cry = m_geo->getReciDir(it->h,it->k,it->l);
        aFamily.eqhkl_normals.push_back(m_geo->getReciVecInRotCry(cry));
      }
    }
  }
  else
  {
    for(HKLList::const_iterator it_hkl=m_cinfo->hklBegin();it_hkl!=it_hklE;++it_hkl)
    {
      double wlthr = calcWavelengthThreshold(CR,it_hkl->dspacing,it_hkl->fsquared);
#if __cplusplus >= 201103L
      m_Reflections->family_list.emplace_back(it_hkl->h,it_hkl->k,it_hkl->l,it_hkl->fsquared,it_hkl->dspacing,wlthr);
#else
      m_Reflections->family_list.push_back(ReflectionFamily(it_hkl->h,it_hkl->k,it_hkl->l,it_hkl->fsquared,it_hkl->dspacing,wlthr));
#endif
      ReflectionFamily& aFamily = m_Reflections->family_list.back();
      aFamily.eqhkl_normals.reserve(it_hkl->demi_normals.size());
      for (std::vector<HKLInfo::Normal>::const_iterator it = it_hkl->demi_normals.begin(); it != it_hkl->demi_normals.end(); ++it)
      {
        aFamily.eqhkl_normals.push_back(m_geo->getReciVecInRotCry(Vector(it->x, it->y, it->z)));
      }
    }
  }
  m_Reflections->shrink_to_fit();
}

void NCrystal::SCBragg::pimpl::updateCache(double wl, const NCrystal::Vector& dir ) const
{
  double dot;
  if ( ncabs(m_cache.wl-wl)<1.0e-9 && (dot=dir.dot(m_cache.dir))>0 ) {
    //wl compatible and dot product is positive.
    double dir_mag2 = dir.mag2();
    if ( dot*dot > dir_mag2*( 1.0 - 1.0e-9 ) ) {//check that cos(angle[m_cache.dir,dir]) ~= 1
      //cache already valid!
      return;
    }
    //not valid!
    m_cache.dir = dir;
    if (dir_mag2!=1.0)
      m_cache.dir *= 1.0/sqrt(dir_mag2);
  } else {
    //not valid!
    m_cache.dir = dir;
    m_cache.dir.normalise();
  }

  //wavelength or direction is new, we must recalculate. Note that m_cache.dir was set to dir.unit() during the check above.

  nc_assert(m_cache.wl!=NCSCBragg_LACKSORIENTATION);

  m_cache.wl = wl;
  m_cache.contribs.clear();
  m_cache.famchan.clear();

  if (wl>m_threshold_wl)
    return;

  const unsigned nfam(m_Reflections->family_list.size());
  double xs = 0.0;
  for(unsigned i=0; i<nfam ;++i)
    {
      if (!channelXSInLab(i, wl, m_cache.dir, xs))
        break;
    }
}

void NCrystal::SCBragg::domain(double& ekin_low, double& ekin_high) const
{
  ekin_low = m_pimpl->m_threshold_ekin;
  ekin_high = infinity;
}

double NCrystal::SCBragg::crossSection(double ekin, const double (&indir)[3] ) const
{
  m_pimpl->updateCache(ekin2wl(ekin), asVect(indir));
  return m_pimpl->m_cache.contribs.empty() ? 0.0 : m_pimpl->m_cache.contribs.back();
}

bool NCrystal::SCBragg::pimpl::channelXSInLab(unsigned index ,
                                              double wavelength,
                                              const NCrystal::Vector& indir,
                                              double &xs_accumulate) const
{
  nc_assert(wavelength>=0);
  const ReflectionFamily* hklInfo = &m_Reflections->family_list[index];

  const double sin_perfect_theta = wavelength * hklInfo->inv2d;

  //TODO for NC2: The threshold here is 1e-5, similar to the fsquared threshold
  //in .nxs/.ncmat factories. We should revisit if we want this threshold as
  //NCMatCfg:
  const double contribution_threshold = 1.0e-5;//barn [NB: hardcoded above as well!]


  //diffraction is impossible
  if(sin_perfect_theta >= 1.)// == 1.0 is backscattering.
    return false;

  if (wavelength<hklInfo->wlthr) {
    //Q(wavelength)<=contribution_threshold, so no way that Q*W can ever exceed
    //the threshold.
    return true;
  }

  const double diff_perfect_distance_sq = 2.0 * ( 1.0 - sin_perfect_theta );


  bool needPT(true);
  double perfectTheta(0.), Q(0.), phi(0.);

  std::vector<Vector>::const_iterator itE = hklInfo->eqhkl_normals.end();
  std::vector<Vector>::const_iterator itB = hklInfo->eqhkl_normals.begin();
  std::vector<Vector>::const_iterator it = itB;
  for(;it!=itE;++it)
    {
      const Vector& normInLab = *it;

      //neutron inverse direction should be on the same side of a reflection plane with
      //the plane normal, so that reflected direction is also on that side.
      //remove 50% of the reflections on average
#ifndef NCRYSTAL_EQREFL_SKIPHALF
      if(normInLab.dot(indir) >= 0.)
        continue;
      const double normsign = 1.0;
#else
      const double dot = normInLab.dot(indir);
      double normsign = 1.0;
      if (dot>=0) {
        if (!dot)
          continue;
        normsign = -1.0;
      }
#endif

      //We expand (normInLab*normsign+indir).mag2() for the benefit of non-opt builds;
      double dr_x = normsign*normInLab.x()+indir.x();
      double dr_y = normsign*normInLab.y()+indir.y();
      double dr_z = normsign*normInLab.z()+indir.z();
      double diff_real_distance_sq = dr_x*dr_x + dr_y*dr_y + dr_z*dr_z;

      //Perform the test
      //  ncabs(sqrt(diff_perfect_distance_sq)-sqrt(diff_real_distance_sq)) > sqrt(m_max_distance_mag2)
      //in a less expensive manner:
      const double x = diff_real_distance_sq + diff_perfect_distance_sq - m_max_distance_mag2;
      if (x >= 0 && x*x > 4*diff_real_distance_sq*diff_perfect_distance_sq)
        continue;

      if (needPT) {
        //only need these calculations at most once per family (depends on wl/2d, not normal).
        //NB: Using sin(2x) = 2*sin(x)*sqrt(1-sin(x)*sin(x)) below to avoid asin+sin calls.
        double tmp1 = wavelength*wavelength*wavelength * hklInfo->Fsqr;
        double tmp2 = 2.0*sin_perfect_theta*m_Reflections->V0 * m_Reflections->numAtom;
        double Q2 = tmp1 * tmp1 / (tmp2*tmp2*(1.0-sin_perfect_theta*sin_perfect_theta));
        //No gain from checking here that Q>=contribution_threshold, since the
        //wavelength threshold above already guarantees that.
        Q = sqrt(Q2);
        perfectTheta = std::asin(sin_perfect_theta);
        phi = M_PI_2-perfectTheta;
        needPT = false;
      }
      double angle_ki_nominal_nor = acos( 1 - 0.5*diff_real_distance_sq );
      double delta_theta=angle_ki_nominal_nor-phi;
      double W = m_gauss_norm * std::exp(m_neginv_2_gauss_sigma_sq*delta_theta*delta_theta);
      double contribution = Q*W;

      if(m_enable_backscattering)
        {
          if(perfectTheta > M_PI_2-2*m_gauss_sigma) //XX. the cutting value "M_PI/2.-2*m_gauss_sigma" was determined by observation
            {
              //TK: Where do we apply the gauss truncation to 5 sigma here??

              //nominal normal is outside the cone
              if(delta_theta > 0)
                {
                  contribution *= ncerf(phi * m_sqrt2_div_sigma);
                }
              else //inside the cone
                {
                  contribution *= ncerf(ncabs(m_sqrt2_div_sigma * (phi+0.5*delta_theta)));
                }
            }
        }

      if (contribution<contribution_threshold)
        continue;

      xs_accumulate += contribution;
      m_cache.contribs.push_back(xs_accumulate);
      m_cache.famchan.push_back(Cache::ContribInfo(index,it-itB,normsign));
    }

  return true;

}

void NCrystal::SCBragg::generateScattering( double ekin, const double (&indir)[3],
                                            double (&outdir)[3], double& de ) const

{
  m_pimpl->genScat(ekin2wl(ekin),asVect(indir),asVect(outdir),de);
}

void NCrystal::SCBragg::pimpl::genScat( double wl, const NCrystal::Vector& indir,
                                        NCrystal::Vector& outdir, double& de ) const
{
  de = 0.0; //always elastic

  updateCache(wl, indir);

  if (m_cache.contribs.empty()) {
    //Scatterings not actually possible at this configuration, so this should
    //only happen if the user did not check for a positive crossSection(..)
    //before calling generateScattering. For lack of a better option, we simply
    //generate an isotropic scattering.
    m_scb->randIsotropicDirection(NC_VECTOR_CAST(outdir));
    return;
  }

  //select a reflection plane at random:
  double rand_contrib = m_scb->rand() * m_cache.contribs.back();

  std::vector<double>::const_iterator it = m_cache.contribs.begin();
  std::vector<double>::const_iterator itE = m_cache.contribs.end();
  //Use binary search to find which plane was selected at random:
  if ( ( it = std::lower_bound( it, itE, rand_contrib ) ) == itE )
    NCRYSTAL_THROW(CalcError,"Inconsistency encountered while selecting reflection plane.");

  const Cache::ContribInfo& ci = m_cache.famchan[it-m_cache.contribs.begin()];

  ReflectionFamily * family = &(m_Reflections->family_list[ci.famidx]);
  const Vector& normalInLab = family->eqhkl_normals[ci.chanidx];

  //Based on m_delta_d and m_gauss_sigma, generate theta_bragg and opening angle
  //of reflection cone (the cone where we must select a random point to get the
  //final reflected direction):

  double d_diff, cone_opening_angle;
  if (m_delta_d) {
    m_scb->randNorm(d_diff,cone_opening_angle);
    d_diff = 1.0 + m_delta_d*d_diff;
  } else {
    //m_delta_d is zero, so no need to sample norm twice:
    m_scb->randNorm(cone_opening_angle);
    d_diff = 1.0;
  }
  cone_opening_angle *= m_gauss_sigma;

  //Note, that since we are generating two random numbers at once, and in
  //principle using just 1 of them if m_delta_d==0, one might be tempted to
  //cache one for the next call. However, such caching is potentially very
  //confusing for users of our library since there would essentially be a hidden
  //random state in addition to the one of their random number generator, thus
  //messing with any reproducibility requirements.

  //Calculate reflected bragg angle, taking into account d_diff:
  double sinthetabragg = wl * family->inv2d/d_diff;

  //Finally, sample a point on the reflection cone and update outdir:
  sampleConeAndReflect(cone_opening_angle,sinthetabragg,ci.normsign,normalInLab,m_cache.dir,outdir);
}

void NCrystal::SCBragg::pimpl::sampleConeAndReflect( double cone_opening_angle, double sinthetabragg, double normsign,
                                                     const NCrystal::Vector& normal, const NCrystal::Vector& indir,
                                                     NCrystal::Vector& outdir ) const
{
  //The cone is parameterized as X(A)=C + cos(th_bragg)* [ cos(A)*U + sin(A)*V ],
  //where V and U are unit vectors, C is the cone centre, A is the cone
  //opening angle.

  const double costhetabragg = std::sqrt(1.0-sinthetabragg*sinthetabragg);//valid for angle in -pi..pi

  Vector V = indir;
  V.cross_inplace(normal);
  V *= (- normsign / V.mag());
  Vector U = V;
  U.cross_inplace(indir);
  Vector C = indir;
  C *= (-1.0 * sinthetabragg / C.mag());
  double cosA,sinA;
  sincos(cone_opening_angle,cosA,sinA);
  U *= -1.0*cosA*costhetabragg;
  V *= sinA*costhetabragg;
  C += U;
  C += V;
  C *= -2.0 * indir.dot(C);
  outdir = indir;
  outdir += C;
}

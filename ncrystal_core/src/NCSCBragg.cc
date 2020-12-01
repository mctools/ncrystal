
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

#include "NCrystal/internal/NCSCBragg.hh"
#include "NCrystal/internal/NCGaussMos.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include "NCrystal/NCSCOrientation.hh"
#include "NCrystal/NCInfo.hh"
#include "NCrystal/internal/NCVector.hh"
#include "NCrystal/internal/NCOrientUtils.hh"
#include "NCrystal/internal/NCPlaneProvider.hh"
#include <functional>//std::greater
namespace NC=NCrystal;

//magic values:
#define NCSCBragg_LACKSORIENTATION (-2.0)
#define NCSCBragg_INVALIDATECACHE (-1.0)

struct NC::SCBragg::pimpl {

  struct ReflectionFamily {

    //A familiy is here taken to be all planes sharing d-spacing and fsquared.

    std::vector<Vector> deminormals;
    double xsfact;// = fsquared / (unit_cell_volume * unit_cell_natoms)
    double inv2d;

    ReflectionFamily(double xsfct, double dspacing)
      : xsfact(xsfct), inv2d(0.5/dspacing) { nc_assert(xsfct>0&&dspacing>0); }
    ~ReflectionFamily() = default;

    //enable move assignment/construction:
    ReflectionFamily & operator= ( ReflectionFamily && ) = default;
    ReflectionFamily( ReflectionFamily && ) = default;
    //disable expensive assignment/copy construction:
    ReflectionFamily & operator= ( const ReflectionFamily & ) = delete;
    ReflectionFamily( const ReflectionFamily & ) = delete;

    bool operator < ( const ReflectionFamily & o ) const
    {
      //sort by d-spacing (secondarily by xsfact for reproducibility):
      if ( o.inv2d!=inv2d ) return o.inv2d > inv2d;
      nc_assert(o.xsfact != xsfact);
      return o.xsfact < xsfact;
    }
  };

  typedef std::map<std::pair<uint64_t,uint64_t>,std::vector<Vector>,
                   std::greater<std::pair<uint64_t,uint64_t> > > SCBraggSortMap;

  pimpl( const NC::Info* cinfo, double mosaicity, double dd,
         const SCOrientation& sco, PlaneProvider * plane_provider,
         double prec, double ntrunc );
  ~pimpl() = default;

  double setupFamilies( const Info * cinfo,
                        const RotMatrix& cry2lab,
                        PlaneProvider * plane_provider,
                        double V0numAtom );
  void genScat( const SCBragg *, Vector& outdir ) const;
  void updateCache(double wl, const Vector& ) const;

  struct Cache {
    Cache() : ekin(NCSCBragg_LACKSORIENTATION) {}
    ~Cache() = default;
    //cache signature:
    double ekin;
    Vector dir;
    //cache contents:
    double wl;
    VectD xs_commul;
    std::vector<GaussMos::ScatCache> scatcache;
  };

  double m_threshold_ekin;
  std::vector<ReflectionFamily> m_reflfamilies;
  GaussMos m_gm;
  //Cache - should be in tread-local storage if calling this in a multi-threaded
  //application (must then also cache unique-id of SCBragg object filling the
  //cache, since thread-local storage implies static storage):
  mutable Cache m_cache;
};

NC::SCBragg::pimpl::pimpl(const NC::Info* cinfo, double mosaicity,
                          double dd, const SCOrientation& sco, PlaneProvider * plane_provider,
                          double prec, double ntrunc)
  : m_threshold_ekin(kInfinity),
    m_gm(mosaicity,true/*mos is FWHM*/,prec,ntrunc)
{
  nc_assert_always(cinfo);
  m_gm.setDSpacingSpread(dd);

  //Start with invalid cache:
  m_cache.ekin = NCSCBragg_LACKSORIENTATION;

  //Always needs structure info:
  if (!cinfo->hasStructureInfo())
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks Structure information.");

  //Setup based on structure info:
  RotMatrix reci_lattice = getReciprocalLatticeRot( *cinfo );
  RotMatrix cry2lab = getCrystal2LabRot( sco, reci_lattice );
  double V0numAtom = cinfo->getStructureInfo().n_atoms * cinfo->getStructureInfo().volume;

  double maxdsp = setupFamilies( cinfo, cry2lab, plane_provider, V0numAtom );

  m_threshold_ekin = wl2ekin(maxdsp * 2.0);

}

NC::SCBragg::SCBragg( const NC::Info* cinfo,
                      const SCOrientation& sco,
                      double mosaicity,
                      double dd,
                      PlaneProvider * plane_provider,
                      double prec, double ntrunc)
  : Scatter("SCBragg"),
    m_pimpl(std::make_unique<pimpl>(cinfo,mosaicity,dd,sco,plane_provider,prec,ntrunc))
{
  validate();
}

NC::SCBragg::~SCBragg() = default;

double NC::SCBragg::pimpl::setupFamilies( const NC::Info * cinfo,
                                          const NC::RotMatrix& cry2lab,
                                          NC::PlaneProvider * plane_provider,
                                          double V0numAtom )
{
  m_cache.ekin = NCSCBragg_INVALIDATECACHE;//Invalidate cache and note that we were initialised

  //expand crystal info
  nc_assert_always(cinfo->hasHKLInfo());
  nc_assert_always(cinfo->hasStructureInfo());
  nc_assert(m_reflfamilies.empty());

  //collect all planes, sorted by (dsp,fsq). To avoid issues connected to
  //floating point number keys, we store dspacing/fsquared as integers, keeping
  //precision down to 1e-10 angstrom and 1e-10 barn respectively.
  SCBraggSortMap planes;
  //but also try to avoid rounding issues when floating point values are not misbehaving:
  std::map<uint64_t,double> origvals_dsp;
  std::map<uint64_t,double> origvals_fsq;

  double dsp, fsq;
  Vector demi_normal;
  const double two30 = 1073741824.0;//2^30 ~= 1.07e9

  std::unique_ptr<PlaneProvider> ppguard;
  if (!plane_provider) {
    //fall back to standard plane provider
    ppguard = createStdPlaneProvider(cinfo);
    plane_provider = ppguard.get();
  } else {
    //use supplied plane provider - we must reset looping since we don't know
    //its current state:
    plane_provider->prepareLoop();
  }

  double maxdspacing(0);

  while (plane_provider->getNextPlane(dsp, fsq, demi_normal)) {

    if (dsp>maxdspacing)
      maxdspacing = dsp;

    nc_assert(dsp>0.0&&fsq>0.0&&dsp<1e7&&fsq<1e7);
    uint64_t ui_dsp = (uint64_t)(dsp*two30+0.5);
    uint64_t ui_fsq = (uint64_t)(fsq*two30+0.5);

    //a bit messy, but nice to preserve values when possible:
    std::map<uint64_t,double>::iterator itOrig = origvals_dsp.find(ui_dsp);
    if (itOrig==origvals_dsp.end()) {
      origvals_dsp[ui_dsp] = dsp;
    } else if (ncabs(dsp-itOrig->second)>1e-12) {
      itOrig->second = -1;//multiple values observed ...!
    }
    itOrig = origvals_fsq.find(ui_fsq);
    if (itOrig==origvals_fsq.end()) {
      origvals_fsq[ui_fsq] = fsq;
    } else if (ncabs(fsq-itOrig->second)>1e-12) {
      itOrig->second = -1;//multiple values observed ...!
    }

    std::pair<uint64_t,uint64_t> key(ui_dsp,ui_fsq);

    SCBraggSortMap::iterator it = planes.find(key);
    if ( it != planes.end() ) {
      it->second.push_back(demi_normal);
    } else {
      std::pair<PairDD,std::vector<Vector> > newentry;
      newentry.first = key;
      newentry.second.push_back(demi_normal);
      planes.insert(it,newentry);
    }
  }

  m_reflfamilies.reserve(planes.size());
  SCBraggSortMap::const_iterator it = planes.begin();
  for (;it!=planes.end();++it) {

    std::map<uint64_t,double>::iterator itOrig = origvals_dsp.find(it->first.first);
    nc_assert(itOrig!=origvals_dsp.end());
    dsp = (itOrig->second > 0 ? itOrig->second : it->first.first / two30);

    itOrig = origvals_fsq.find(it->first.second);
    nc_assert(itOrig!=origvals_fsq.end());
    fsq = (itOrig->second > 0 ? itOrig->second : it->first.second / two30);

    m_reflfamilies.emplace_back(fsq/V0numAtom,dsp);

    //transfer it->second into final vector and put in the lab frame:
    ReflectionFamily& fam = m_reflfamilies.back();
    fam.deminormals.reserve(it->second.size());
    std::vector<Vector>::const_iterator itDN(it->second.begin()), itDNE(it->second.end());
    for (;itDN!=itDNE;++itDN)
      fam.deminormals.push_back( cry2lab * (*itDN) );
    nc_assert(m_reflfamilies.back().deminormals.size()==it->second.size());
  }

  return maxdspacing;
}


namespace NCrystal {
  inline double SCBragg_cacheRound(double x) {
    //Cut off input at 15 decimals, which should be a negligible effect on any
    //realistic value of ekin in eV, but ensures that we don't
    //get call-order irreproducibilities.
    return std::floor(x*1e15+0.5)*1e-15;
  }
}

void NC::SCBragg::pimpl::updateCache(double ekin_raw, const NC::Vector& dir ) const
{
  //We check the cache validity on the rounded ekin value, but for simplicity we
  //keep the direction as it is. We could consider rounding the direction as
  //well...
  //NB: We used to check co-alignment of angles via a dot-product, but that is
  //actually numerically imprecise for small angles, leading to occurances of
  //cache validity where it should have been invalid.
  double ekin = SCBragg_cacheRound(ekin_raw);
  if ( m_cache.ekin==ekin && dir.angle_highres(m_cache.dir)<1.0e-12 ) {
    //cache already valid!
    return;
  }

  //Cache not valid!
  m_cache.dir = dir;
  m_cache.dir.normalise();

  //Energy or direction is new, we must recalculate. Note that m_cache.dir was
  //normalised during the check above.

  nc_assert(m_cache.ekin!=NCSCBragg_LACKSORIENTATION);
  m_cache.ekin = ekin;
  m_cache.wl = ekin2wl(ekin);
  nc_assert(m_cache.wl>=0);
  m_cache.scatcache.clear();
  m_cache.xs_commul.clear();
  if (m_cache.wl==0)
    return;//done, all cross-sections will be zero

  std::vector<ReflectionFamily>::const_iterator it(m_reflfamilies.begin()), itE(m_reflfamilies.end());

  double inv2dcutoff = (1.0-2*std::numeric_limits<double>::epsilon())/m_cache.wl;

  GaussMos::InteractionPars interactionpars;
  for( ; it!=itE; ++it) {
    const ReflectionFamily& fam = *it;
    if( fam.inv2d >= inv2dcutoff )
      break;//stop here, no more families fulfill w<2d requirement.
    interactionpars.set(m_cache.wl, fam.inv2d, fam.xsfact);
    m_gm.calcCrossSections(interactionpars, m_cache.dir, fam.deminormals, m_cache.scatcache,m_cache.xs_commul);
  }

  nc_assert(m_cache.xs_commul.empty()||m_cache.xs_commul.back()>0.0);
}

void NC::SCBragg::pimpl::genScat( const SCBragg* scb, NC::Vector& outdir ) const
{
  nc_assert(!m_cache.xs_commul.empty());
  nc_assert(m_cache.xs_commul.back()>0.0);
  nc_assert(m_cache.xs_commul.size()==m_cache.scatcache.size());

  RandomBase * rng = scb->getRNG();
  std::size_t idx = pickRandIdxByWeight(rng,m_cache.xs_commul);
  nc_assert(idx<m_cache.scatcache.size());
  GaussMos::ScatCache& chosen_scatcache = m_cache.scatcache[idx];

  m_gm.genScat( rng, chosen_scatcache, m_cache.wl, m_cache.dir, outdir );
}

void NC::SCBragg::domain(double& ekin_low, double& ekin_high) const
{
  ekin_low = m_pimpl->m_threshold_ekin;
  ekin_high = kInfinity;
}

double NC::SCBragg::crossSection(double ekin, const double (&indir)[3] ) const
{
  if ( ekin <= m_pimpl->m_threshold_ekin )
    return 0.0;
  m_pimpl->updateCache(ekin, asVect(indir));
  return m_pimpl->m_cache.xs_commul.empty() ? 0.0 : m_pimpl->m_cache.xs_commul.back();
}


void NC::SCBragg::generateScattering( double ekin, const double (&indir)[3],
                                      double (&outdir)[3], double& de ) const

{
  de = 0;

  if ( ekin <= m_pimpl->m_threshold_ekin ) {
    //Scatterings not actually possible at this configuration, so don't change direction:
    asVect(outdir) = asVect(indir);
    return;
  }

  m_pimpl->updateCache(ekin, asVect(indir));

  if (m_pimpl->m_cache.xs_commul.empty()||m_pimpl->m_cache.xs_commul.back()<=0.0) {
    //Again, scatterings are not actually possible here:
    asVect(outdir) = asVect(indir);
    return;
  }

  m_pimpl->genScat(this,asVect(outdir));
}

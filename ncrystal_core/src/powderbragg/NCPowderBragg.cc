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

#include "NCrystal/internal/powderbragg/NCPowderBragg.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <functional>//std::greater

namespace NC = NCrystal;
namespace NCRYSTAL_NAMESPACE {
  namespace {
    constexpr double dspacing_merge_tolerance = 1e-11;
    class CachePowderBragg : public CacheBase {
    public:
      void invalidateCache() override { ekin.dbl() = -1.0; }
      CachePowderBragg() { ekin.dbl() = -1.0; }
      bool cacheOK( NeutronEnergy eee ) const
      {
        nc_assert( eee.dbl() > 0.0 );
        nc_assert( std::isfinite(eee.dbl()) );
        return eee == ekin;
      }
      void updateCache( NeutronEnergy eee,
                        std::size_t lvpi )
      {
        nc_assert( eee.dbl() > 0.0 );
        nc_assert( std::isfinite(eee.dbl()) );
        ekin = eee;
        inv_ekin = 1.0 / eee.dbl();
        lastValidPlaneIdx = lvpi;
      }
      NeutronEnergy ekin;
      double inv_ekin;
      std::size_t lastValidPlaneIdx;
    };

  }
}

void NC::PowderBragg::init( const StructureInfo& si, VectDFM&& data )
{
  nc_assert_always(si.n_atoms>0);
  nc_assert_always(si.volume>0);
  if (!(si.volume>0) || !(si.n_atoms>=1) )
    NCRYSTAL_THROW(BadInput,"Passed structure info object has"
                   " invalid volume or n_atoms fields.");
  init(si.volume * si.n_atoms, std::move(data));
}

void NC::PowderBragg::init( double v0_times_natoms, VectDFM&& origdata )
{
  if (!(v0_times_natoms>0) )
    NCRYSTAL_THROW(BadInput,"v0_times_natoms is not a positive number.");
  double xsectfact = 0.5/v0_times_natoms;
  xsectfact *= wl2ekin(1.0);//Adjust units so we can get cross sections through
                            //multiplication with 1/ekin instead of wl^2.
  VectDFM data = std::move(origdata);
  std::sort(data.begin(),data.end(),std::greater<PairDD>());
  VectD v2dE;
  v2dE.reserve(data.size());
  VectD fdm_commul;
  fdm_commul.reserve(data.size());
  StableSum fdmsum2;
  VectDFM::const_iterator it(data.begin()),itE(data.end());
  double prev_dsp = -kInfinity;
  for (;it!=itE;++it) {
    if (!(it->first>0.0))
      NCRYSTAL_THROW(CalcError,"Inconsistent plane data implies "
                     "non-positive (or NaN) d_spacing.");
    if ( ncabs(prev_dsp-it->first) < dspacing_merge_tolerance ) {
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
    fdm_commul.clear();
    v2dE.clear();
  } else {
    m_threshold = NeutronEnergy{ *(v2dE.begin()) };
  }
  //Transfer while squeezing memory:
  VectD(fdm_commul.begin(),fdm_commul.end()).swap(m_fdm_commul);
  VectD(v2dE.begin(),v2dE.end()).swap(m_2dE);
  nc_assert( m_threshold.get() > 0.0 );
}

NC::PowderBragg::PowderBragg( const StructureInfo& si, VectDFM&&  data)
{
  init(si,std::move(data));
}

NC::PowderBragg::PowderBragg( double v0_times_natoms, VectDFM&&  data)
{
  init(v0_times_natoms,std::move(data));
}

NC::PowderBragg::PowderBragg(const Info&ci)
{
  if (!ci.hasHKLInfo())
    NCRYSTAL_THROW(MissingInfo,"Passed Info object lacks HKL information.");
  if (!ci.hasStructureInfo())
    NCRYSTAL_THROW(MissingInfo,
                   "Passed Info object lacks Structure information.");
  const auto& hklList = ci.hklList();
  VectDFM data;
  data.reserve(hklList.size());

  for ( const auto& hkl : hklList ) {
    double f = hkl.fsquared * hkl.multiplicity;
    if (f<0)
      NCRYSTAL_THROW(CalcError,
                     "Inconsistent data implies negative |F|^2*multiplicity.");
    if (data.empty()||data.back().first!=hkl.dspacing) {
      data.emplace_back(hkl.dspacing,f);
    } else {
      data.back().second += f;
    }
  }
  init(ci.getStructureInfo(),std::move(data));
}

NC::EnergyDomain NC::PowderBragg::domain() const noexcept
{
  return { m_threshold, NeutronEnergy{kInfinity} };
}

std::size_t NC::PowderBragg::findLastValidPlaneIdx( NC::NeutronEnergy ekin) const {
  //Quick binary search to find index of the plane with the smallest d-spacing
  //satisfying wl<=2d, but in energy-space: Finding the index of the plane with
  //the largest value of ekin2wl(2d) satisfying ekin>=ekin2wl(2d).  We already
  //know that ekin>=m_2dE[0], so we search from one past this entry:
  nc_assert( !ncisnan(ekin.dbl()) );
  nc_assert( ekin >= m_threshold );
  nc_assert( std::isfinite( ekin.dbl() ) );
  return (std::upper_bound(m_2dE.begin() + 1,
                           m_2dE.end(),
                           ekin.get()) - m_2dE.begin()) - 1;
}


NC::CrossSect NC::PowderBragg::crossSectionIsotropic( NC::CachePtr& cp,
                                                      NC::NeutronEnergy ekin ) const
{
  if ( ekin < m_threshold || !std::isfinite(ekin.dbl()) )
    return CrossSect{0.0};
  auto& cache = accessCache<CachePowderBragg>(cp);
  if (!cache.cacheOK(ekin))
    cache.updateCache( ekin, findLastValidPlaneIdx(ekin) );
  nc_assert(cache.lastValidPlaneIdx<m_fdm_commul.size());
  return CrossSect{ m_fdm_commul[cache.lastValidPlaneIdx] * cache.inv_ekin };
}

NC::CosineScatAngle NC::PowderBragg::genScatterMu( RNG& rng,
                                                   NeutronEnergy ekin,
                                                   std::size_t last_valid_idx) const
{
  nc_assert( ekin >= m_threshold );
  nc_assert( std::isfinite( ekin.dbl() ) );
  nc_assert(last_valid_idx<m_fdm_commul.size());

  //randomly select one plane by contribution:
  VectD::const_iterator itFCUpper = std::next( m_fdm_commul.begin(), last_valid_idx );
  VectD::const_iterator itFC = std::lower_bound( m_fdm_commul.begin(),
                                                 itFCUpper,
                                                 rng.generate() * (*itFCUpper) );
  std::size_t idx_rand = (std::size_t)( itFC - m_fdm_commul.begin() );
  nc_assert(idx_rand<m_2dE.size());
  double sin_theta_bragg_squared = m_2dE[idx_rand] / ekin.get();

  //scatter angle A=2*theta_bragg, so with x=sin^2(theta_bragg), we have:
  //   x = sin^2(A/2)= (1-cosA)/2 => 1-2x = cosA = mu
  const double mu = 1.0 - 2.0 * sin_theta_bragg_squared;
  nc_assert(ncabs(mu)<=1.0);
  return CosineScatAngle{mu};
}

NC::ScatterOutcomeIsotropic NC::PowderBragg::sampleScatterIsotropic( NC::CachePtr& cp,
                                                                     NC::RNG& rng,
                                                                     NC::NeutronEnergy ekin ) const
{
  //elastic: ekin unchanged
  if ( ekin < m_threshold || !std::isfinite(ekin.dbl()) ) {
    //scatterings not possible here
    return { ekin, CosineScatAngle{1.0} };
  } else {
    auto& cache = accessCache<CachePowderBragg>(cp);
    if (!cache.cacheOK(ekin))
      cache.updateCache( ekin, findLastValidPlaneIdx(ekin) );
    return { ekin, genScatterMu(rng,ekin,cache.lastValidPlaneIdx) };
  }
}

std::shared_ptr<NC::ProcImpl::Process> NC::PowderBragg::createMerged( const Process& oraw,
                                                                      double scale1,
                                                                      double scale2 ) const
{
  nc_assert(scale1>0.0);
  nc_assert(scale2>0.0);
  auto optr = dynamic_cast<const PowderBragg*>(&oraw);
  if (!optr)
    return nullptr;
  auto& o = *optr;

  auto result = std::make_shared<PowderBragg>( no_init );//empty instance
  auto fixThreshold = [&result]()
  {
    result->m_threshold = NeutronEnergy{ result->m_2dE.front() };
  };

  //transfer "a" (2dE) and "b" (fdm_commul) vectors, sorted by a:
  VectD& new_a = result->m_2dE;
  VectD& new_b = result->m_fdm_commul;
  nc_assert( new_a.empty() && new_b.empty() );

  const auto& old1_a = this->m_2dE;
  const auto& old1_b = this->m_fdm_commul;
  const auto& old2_a = o.m_2dE;
  const auto& old2_b = o.m_fdm_commul;
  nc_assert(old1_a.size()==old1_b.size());
  nc_assert(old2_a.size()==old2_b.size());
  new_a.reserve(old1_a.size()+old2_a.size());
  new_b.reserve(old1_b.size()+old2_b.size());

  //Special case empty vectors:
  if ( old1_a.empty() ) {
    new_a = old2_a;
    new_b = vectorTrf(old2_b, [scale2](double x) {return x*scale2;});
    return fixThreshold(), result;
  }
  if ( old2_a.empty() ) {
    new_a = old1_a;
    new_b = vectorTrf(old1_b, [scale1](double x) {return x*scale1;});
    return fixThreshold(), result;
  }

  //Merge lists, sort so m_2dE is ordered by increasing magnitude. And keep in
  //mind that the fdm (_b) vectors are commulative! Try to do it without
  //numerical issues related to subtraction.
  std::size_t i1(0), i1E(old1_a.size());
  std::size_t i2(0), i2E(old2_a.size());

  auto can_merge = []( double a_2dE, double b_2dE ) {
    //NB: could use ekin2wlsq here and save some sqrt's.
    double a_d = 0.5*ekin2wl(a_2dE);
    double b_d = 0.5*ekin2wl(b_2dE);
    return ncabs(a_d-b_d)<dspacing_merge_tolerance;
  };

  auto extractFDM = [](const VectD& commulFDM, std::size_t idx)
  {
    nc_assert(idx<commulFDM.size());
    //The subtraction here is where we could potentially introduce numerical
    //errors, but they seem manageable in practice.
    //TODO: can we do better?  Perhaps cache the 100(?) FDM values at highest
    //d-spacing directly in separate vector?
    return idx ? commulFDM.at(idx) - commulFDM.at(idx-1) : commulFDM.front();
  };

  StableSum new_commulFDMSum;
  auto appendFDMPoint = [&new_commulFDMSum,&extractFDM,&new_b]( const VectD& commulFDM,
                                                                std::size_t idx,
                                                                double scale )
  {
    new_commulFDMSum.add( scale * extractFDM(commulFDM,idx) );
    new_b.push_back(new_commulFDMSum.sum());
  };

  nc_assert(old1_a.back() >= old1_a.front());
  while ( i1 < i1E && i2 < i2E ) {
    if ( can_merge( old1_a.at(i1), old2_a.at(i2) ) ) {
      //Same d-spacing point present in both lists:
      new_a.push_back( 0.5 * ( old1_a.at(i1)+old2_a.at(i2) ) );
      new_commulFDMSum.add(scale1*extractFDM(old1_b,i1));
      new_commulFDMSum.add(scale2*extractFDM(old2_b,i2));
      new_b.push_back(new_commulFDMSum.sum());
      ++i1;
      ++i2;
      continue;
    }
    if ( old1_a.at(i1) < old2_a.at(i2) ) {
      new_a.push_back(old1_a.at(i1));
      appendFDMPoint(old1_b,i1,scale1);
      ++i1;
    } else {
      new_a.push_back(old2_a.at(i2));
      appendFDMPoint(old2_b,i2,scale2);
      ++i2;
    }
  }
  //Transfer any remaining entries:
  while ( i1 < i1E ) {
    new_a.push_back(old1_a.at(i1));
    appendFDMPoint(old1_b,i1,scale1);
    ++i1;
  }
  while ( i2 < i2E ) {
    new_a.push_back(old2_a.at(i2));
    appendFDMPoint(old2_b,i2,scale2);
    ++i2;
  }

  new_a.shrink_to_fit();
  new_b.shrink_to_fit();

  return fixThreshold(), result;
}

NC::Optional<std::string> NC::PowderBragg::specificJSONDescription() const
{
  //Determine max_contrib by looking at the peaks m_2dE:
  double max_contrib(0.0);
  nc_assert(m_2dE.size()==m_fdm_commul.size());
  for ( auto i : ncrange(m_2dE.size()) )
    max_contrib = std::max<double>( max_contrib,m_fdm_commul.at(i) / m_2dE.at(i) );

  std::ostringstream ss;
  {
    std::ostringstream tmp;
    nc_assert(!m_2dE.empty());
    tmp << "nplanes="<<m_2dE.size()
        <<";2dmax="<<m_threshold.wavelength()
        << ";max_contrib="<<CrossSect{max_contrib};
    streamJSONDictEntry( ss, "summarystr", tmp.str(), JSONDictPos::FIRST );
  }
  streamJSONDictEntry( ss, "nhkl", m_2dE.size() );
  streamJSONDictEntry( ss, "max_contrib", max_contrib );
  streamJSONDictEntry( ss, "2dmax", m_threshold.wavelength().dbl(),
                       JSONDictPos::LAST );
  return ss.str();
}

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE

std::pair<NC::CrossSect,NC::ScatterOutcome>
NC::PowderBragg::evalXSAndSampleScatter( CachePtr& cp, RNG& rng,
                                         NeutronEnergy ekin,
                                         const NeutronDirection& dir ) const
{
  if ( ekin < m_threshold || !std::isfinite(ekin.dbl()) ) {
    nc_assert( ekin.dbl()>=0.0 );
    return { CrossSect{0.0},
             { ekin, dir } };
  } else {
    auto& cache = accessCache<CachePowderBragg>(cp);
    if (!cache.cacheOK(ekin))
      cache.updateCache( ekin, findLastValidPlaneIdx(ekin) );
    auto mu = genScatterMu(rng,ekin,cache.lastValidPlaneIdx);
    auto outdir = randNeutronDirectionGivenScatterMu( rng,
                                                      mu.dbl(),
                                                      dir.as<Vector>() );
    return { CrossSect{ m_fdm_commul[cache.lastValidPlaneIdx] * cache.inv_ekin },
             { ekin, outdir } };
  }

}

std::pair<NC::CrossSect,NC::ScatterOutcomeIsotropic>
NC::PowderBragg::evalXSAndSampleScatterIsotropic(CachePtr& cp, RNG& rng,
                                                 NeutronEnergy ekin ) const
{
  if ( ekin < m_threshold || !std::isfinite(ekin.dbl()) ) {
    nc_assert( ekin.dbl()>=0.0 );
    return { CrossSect{0.0},
             { ekin, CosineScatAngle{1.0} } };
  } else {
    auto& cache = accessCache<CachePowderBragg>(cp);
    if (!cache.cacheOK(ekin))
      cache.updateCache( ekin, findLastValidPlaneIdx(ekin) );
    return { CrossSect{ m_fdm_commul[cache.lastValidPlaneIdx] * cache.inv_ekin },
             { ekin, genScatterMu(rng,ekin,cache.lastValidPlaneIdx) } };
  }
}

void NC::PowderBragg::evalManyXSIsotropic( CachePtr&, const double* ekin,
                                           std::size_t N,
                                           double* out_xs ) const
{
  for ( std::size_t i = 0; i < N; ++i ) {
    if ( ekin[i] < m_threshold.dbl() || !std::isfinite(ekin[i]) ) {
      nc_assert( ekin[i] >= 0.0 );
      out_xs[i] = 0.0;
    } else {
      std::size_t idx = findLastValidPlaneIdx(NeutronEnergy{ekin[i]});
      nc_assert(idx<m_fdm_commul.size());
      out_xs[i] = m_fdm_commul[idx] / ekin[i];
    }
  }
}
#endif

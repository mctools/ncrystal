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

#include "NCrystal/internal/extd_utils/NCLCRefModels.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"

namespace NC = NCrystal;

NC::LCBraggRef::LCBraggRef(ProcImpl::ProcPtr scb, LCAxis lcaxis_lab, unsigned nsample)
  : m_sc(std::move(scb)),
    m_lcaxislab(lcaxis_lab.as<Vector>().unit()),
    m_nsample(nsample),
    m_nsampleprime(nsample)
{
  while (!isPrime(m_nsampleprime))
    ++m_nsampleprime;
}

NC::EnergyDomain NC::LCBraggRef::domain() const noexcept
{
  return m_sc->domain();
}

NC::CrossSect NC::LCBraggRef::crossSection(CachePtr& cp, NeutronEnergy ekin, const NeutronDirection& indir_nd ) const
{
  const Vector indir = indir_nd.as<Vector>().unit();
  const Vector lccross = m_lcaxislab.cross(indir);
  const double lcdot = m_lcaxislab.dot(indir);
  StableSum sumxs;
  double dphi = k2Pi / m_nsampleprime;
  for (unsigned i = 0; i<m_nsampleprime; ++i) {
    const PhiRot phirot( i * dphi - kPi );
    auto ndir = phirot.rotateVectorAroundAxis( indir, m_lcaxislab, lccross, lcdot ).as<NeutronDirection>();
    sumxs.add(m_sc->crossSection(cp,ekin,ndir).get());
  }
  return CrossSect{ sumxs.sum()/m_nsampleprime };
}

NC::ScatterOutcome NC::LCBraggRef::sampleScatter(CachePtr& cp, RNG& rng, NeutronEnergy ekin, const NeutronDirection& indir_nd ) const
{
  const Vector indir = indir_nd.as<Vector>().unit();
  const Vector lccross = m_lcaxislab.cross(indir);
  const double lcdot = m_lcaxislab.dot(indir);

  VectD xs;
  std::vector<PhiRot> pr;
  xs.reserve(m_nsample);
  pr.reserve(m_nsample);

  double sumxs = 0.0;

  //Get cross-sections at nsample random phi rotations:
  for (unsigned i = 0; i<m_nsample; ++i) {
    double cosphi,sinphi;
    std::tie(cosphi,sinphi) = randPointOnUnitCircle( rng );
    pr.emplace_back(cosphi,sinphi);
    auto ndir = pr.back().rotateVectorAroundAxis( indir, m_lcaxislab, lccross, lcdot ).as<NeutronDirection>();
    xs.push_back( sumxs += m_sc->crossSection(cp,ekin,ndir).get() );
  }

  if (!sumxs) {
    //no xs, do nothing.
    return { ekin, indir_nd };
  }
  //Select one phi rotation at random:
  const PhiRot& phirot = pr.at(pickRandIdxByWeight(rng,xs));

  //Scatter!
  auto ndir = phirot.rotateVectorAroundAxis( indir, m_lcaxislab, lccross, lcdot ).as<NeutronDirection>();
  auto scatoutcome = m_sc->sampleScatter( cp, rng, ekin, ndir );
  auto outdir = phirot.rotateVectorAroundAxis( scatoutcome.direction.as<Vector>(),
                                               m_lcaxislab, true/*reverse*/).as<NeutronDirection>();
  return { ekin, outdir };
}

NC::LCBraggRndmRot::LCBraggRndmRot(ProcImpl::ProcPtr scb, LCAxis lcaxis_lab, unsigned nsample)
  : m_sc(std::move(scb)),
    m_lcaxislab(lcaxis_lab.as<Vector>().unit()),
    m_nsample(nsample)
{
  nc_assert_always(nsample>0);
}

NC::EnergyDomain NC::LCBraggRndmRot::domain() const noexcept
{
  return m_sc->domain();
}

void NC::LCBraggRndmRot::updateCache(Cache& cache, NeutronEnergy ekin, const Vector& indir) const
{
  cache.neutron_state.first = ekin;
  cache.neutron_state.second = indir;
  cache.rotations.reserve(m_nsample);
  cache.xscommul.reserve(m_nsample);
  cache.rotations.clear();
  cache.xscommul.clear();
  const Vector lccross = m_lcaxislab.cross(indir);
  const double lcdot = m_lcaxislab.dot(indir);
  StableSum sumxs;
  //This reference model is unusual and needs RNG even to calculate cross
  //sections. As a workaround we use a global method to get an RNG stream, which
  //potentially makes this model MT-unsafe (although since getRNG creates an
  //instance unique to the current thread, it will in most cases just behave
  //weirdly and escape user attempts at controlling RNG state, etc.)
  auto rng = getRNG();

  for (unsigned i = 0; i<m_nsample; ++i) {
    double cosphi, sinphi;
    std::tie(cosphi,sinphi) = randPointOnUnitCircle( rng );
    cache.rotations.emplace_back(cosphi, sinphi);
    auto nd = cache.rotations.back().rotateVectorAroundAxis( indir, m_lcaxislab, lccross,
                                                             lcdot ).as<NeutronDirection>();
    sumxs.add(m_sc->crossSection(cache.sc_cacheptr,ekin,nd).get());
    cache.xscommul.push_back(sumxs.sum());
  }
}

NC::CrossSect NC::LCBraggRndmRot::crossSection(CachePtr&cp, NeutronEnergy ekin, const NeutronDirection& indir_nd ) const
{
  const Vector indir = indir_nd.as<Vector>().unit();
  auto& cache = accessCache<Cache>(cp);
  updateCache(cache,ekin,indir);//We always regenerate directions on each cross-section call!
  return CrossSect{ cache.xscommul.back()/m_nsample };
}

NC::ScatterOutcome NC::LCBraggRndmRot::sampleScatter(CachePtr&cp, RNG& rng, NeutronEnergy ekin, const NeutronDirection& indir_nd ) const
{
  const Vector indir = indir_nd.as<Vector>().unit();
  auto& cache = accessCache<Cache>(cp);

  if ( cache.rotations.empty()
       || cache.neutron_state.first != ekin
       || cache.neutron_state.second != indir )
    {
      //trigger generation of random directions and calculate cross sections:
      updateCache(cache,ekin,indir);
    }
  nc_assert(!cache.xscommul.empty());

  if (!cache.xscommul.back()) {
    //no xs, do nothing.
    return {ekin,indir_nd};
  }

  //Select one phi rotation at random:
  const PhiRot& phirot = cache.rotations.at(pickRandIdxByWeight(rng,cache.xscommul));

  //Scatter!
  auto ndir = phirot.rotateVectorAroundAxis( indir, m_lcaxislab).as<NeutronDirection>();
  auto scatoutcome = m_sc->sampleScatter(cache.sc_cacheptr,rng,ekin,ndir);
  auto outdir = phirot.rotateVectorAroundAxis( scatoutcome.direction.as<Vector>(), m_lcaxislab,
                                               true/*reverse*/).as<NeutronDirection>();
  return { ekin, outdir };
}

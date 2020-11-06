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

#include "NCrystal/internal/NCSABScatter.hh"
#include "NCrystal/internal/NCSABFactory.hh"
#include "NCrystal/internal/NCDynInfoUtils.hh"
#include "NCrystal/internal/NCSABFactory.hh"
#include "NCrystal/internal/NCRandUtils.hh"
#include "NCrystal/internal/NCVDOSToScatKnl.hh"
namespace NC = NCrystal;

struct NC::SABScatter::Impl {
  std::shared_ptr<const SAB::SABScatterHelper> m_scathelper_shptr;
};

NC::SABScatter::~SABScatter() = default;

NC::SABScatter::SABScatter( std::shared_ptr<const NC::SAB::SABScatterHelper> sh )
  : ScatterIsotropic("SABScatter"), m_sh(nullptr)
{
  //All other constructors delegate to this one.
  nc_assert_always(!!sh);
  m_impl->m_scathelper_shptr = std::move(sh);
  m_sh = m_impl->m_scathelper_shptr.get();
  nc_assert_always(m_sh);
}

NC::SABScatter::SABScatter( NC::SAB::SABScatterHelper&& sh )
  : SABScatter( std::make_shared<SAB::SABScatterHelper>(std::move(sh)) )
{
}

NC::SABScatter::SABScatter( const DI_ScatKnl& di_sk, unsigned vdoslux, bool useCache )
  : SABScatter( [&di_sk,vdoslux,useCache]()
                {
                  auto sabdata_ptr = extractSABDataFromDynInfo(&di_sk,vdoslux,useCache);
                  nc_assert_always(!!sabdata_ptr);
                  return ( useCache
                           ? SAB::createScatterHelperWithCache( std::move(sabdata_ptr),
                                                                di_sk.energyGrid() )
                           : SAB::createScatterHelper( std::move(sabdata_ptr),
                                                       di_sk.energyGrid() ) );
                }() )
{
}

NC::SABScatter::SABScatter( NC::SABData && sabdata_,
                            const VectD& energyGrid )
  : SABScatter( [&energyGrid](NC::SABData&& sabdata)
                {
                  auto sabdata_shptr = std::make_shared<const SABData>(std::move(sabdata));
                  std::shared_ptr<const VectD> egrid_shptr;
                  if (!energyGrid.empty())
                    egrid_shptr = std::make_shared<const VectD>(energyGrid);
                  return SAB::createScatterHelper( std::move(sabdata_shptr),
                                                   std::move(egrid_shptr) );
                }(std::move(sabdata_)) )
{
}
NC::SABScatter::SABScatter( std::shared_ptr<const SABData> sabdata_shptr,
                            std::shared_ptr<const VectD> egrid_shptr )
  : SABScatter( SAB::createScatterHelper( std::move(sabdata_shptr),
                                          std::move(egrid_shptr) ) )
{
}

double NC::SABScatter::crossSectionNonOriented(double ekin) const
{
  return m_sh->xsprovider.crossSection(ekin);
}

void NC::SABScatter::generateScatteringNonOriented( double ekin, double& angle, double& delta_e ) const
{
  double mu;
  std::tie(delta_e,mu) = m_sh->sampler.sampleDeltaEMu(ekin, *getRNG());
  nc_assert( mu >= -1.0 && mu <= 1.0 );
  angle = std::acos(mu);
}

void NC::SABScatter::generateScattering( double ekin, const double (&indir)[3],
                                         double (&outdir)[3], double& delta_e ) const
{
  double mu;
  RandomBase& rng = *getRNG();
  std::tie(delta_e,mu) = m_sh->sampler.sampleDeltaEMu(ekin, rng);
  nc_assert( mu >= -1.0 && mu <= 1.0 );
  randDirectionGivenScatterMu( &rng, mu, indir, outdir );
}

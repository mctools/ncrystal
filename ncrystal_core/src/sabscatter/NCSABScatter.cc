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

#include "NCrystal/internal/sabscatter/NCSABScatter.hh"
#include "NCrystal/internal/sab/NCSABFactory.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include "NCrystal/internal/sab/NCSABFactory.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/vdos/NCVDOSToScatKnl.hh"
namespace NC = NCrystal;

struct NC::SABScatter::Impl {
  Impl(shared_obj<const SAB::SABScatterHelper> sp) : m_scathelper_shptr(std::move(sp)) {}
  shared_obj<const SAB::SABScatterHelper> m_scathelper_shptr;
};

NC::SABScatter::~SABScatter() = default;

NC::SABScatter::SABScatter( shared_obj<const SAB::SABScatterHelper> sh )
  : m_impl(std::move(sh)), m_sh(m_impl->m_scathelper_shptr.get())
{
  //All other constructors delegate to this one.
}

NC::SABScatter::SABScatter( std::unique_ptr<const SAB::SABScatterHelper> upsh )
  : SABScatter(shared_obj<const SAB::SABScatterHelper>{std::move(upsh)})
{
}

NC::SABScatter::SABScatter( SAB::SABScatterHelper&& sh )
  : SABScatter( makeSO<const SAB::SABScatterHelper>(std::move(sh)) )
{
}

NC::SABScatter::SABScatter( const DI_ScatKnl& di_sk, unsigned vdoslux,
                            bool useCache, uint32_t vdos2sabExcludeFlag )
  : SABScatter( [&di_sk,vdoslux,useCache,vdos2sabExcludeFlag]()
                {
                  auto sabdata_ptr = extractSABDataFromDynInfo(&di_sk,vdoslux,useCache,vdos2sabExcludeFlag);
                  nc_assert_always(!!sabdata_ptr);
                  return ( useCache
                           ? SAB::createScatterHelperWithCache( std::move(sabdata_ptr),
                                                                di_sk.energyGrid() )
                           : SAB::createScatterHelper( std::move(sabdata_ptr),
                                                       di_sk.energyGrid() ) );
                }() )
{
}

NC::SABScatter::SABScatter( SABData && sabdata_, const VectD& energyGrid )
  : SABScatter( [&energyGrid]( SABData&& sabdata)
                {
                  auto sabdata_shptr = makeSO<const SABData>(std::move(sabdata));
                  std::shared_ptr<const VectD> egrid_shptr;
                  if (!energyGrid.empty())
                    egrid_shptr = std::make_shared<const VectD>(energyGrid);
                  return SAB::createScatterHelper( std::move(sabdata_shptr),
                                                   std::move(egrid_shptr) );
                }(std::move(sabdata_)) )
{
}
NC::SABScatter::SABScatter( shared_obj<const SABData> sabdata_shptr,
                            std::shared_ptr<const VectD> egrid_shptr )
  : SABScatter( SAB::createScatterHelper( std::move(sabdata_shptr),
                                          std::move(egrid_shptr) ) )
{
}

NC::CrossSect NC::SABScatter::crossSectionIsotropic( CachePtr&, NeutronEnergy ekin ) const
{
  return CrossSect{ m_sh->xsprovider.crossSection(ekin) };
}

NC::ScatterOutcomeIsotropic NC::SABScatter::sampleScatterIsotropic( CachePtr&, RNG& rng, NeutronEnergy ekin ) const
{
  double delta_e, mu;
  std::tie(delta_e,mu) = m_sh->sampler.sampleDeltaEMu(ekin, rng);
  nc_assert( mu >= -1.0 && mu <= 1.0 );
  return { NeutronEnergy{ncmax(0.0,ekin.get()+delta_e)}, CosineScatAngle{mu} };
}

NC::Optional<std::string> NC::SABScatter::specificJSONDescription() const
{
  return m_sh->specificJSONDescription;
}

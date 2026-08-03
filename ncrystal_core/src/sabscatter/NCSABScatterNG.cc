////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

#include "NCrystal/internal/sabscatter/NCSABScatterNG.hh"
#include "NCrystal/internal/sab/NCSABExtended.hh"
namespace NC = NCrystal;

NC::SABScatterNG::~SABScatterNG() = default;

NC::SABScatterNG::SABScatterNG( shared_obj<const SABExtended> sh, double scale )
  : m_sh(std::move(sh)), m_scale(scale)
{
  nc_assert( m_sh!=nullptr);
  nc_assert( m_scale >= 0.0 );
  nc_assert( std::isfinite(m_scale) );
}

NC::SABScatterNG::SABScatterNG( shared_obj<const SABExtended> sh,
                                SigmaBound sb, double scale )
  : SABScatterNG( std::move(sh), sb.dbl() * scale )
{
}

NC::CrossSect
NC::SABScatterNG::crossSectionIsotropic( CachePtr&, NeutronEnergy ekin ) const
{
  return CrossSect{ m_sh->crossSectionUnitSigmaBound( ekin ).dbl() * m_scale };
}

NC::ScatterOutcomeIsotropic
NC::SABScatterNG::sampleScatterIsotropic( CachePtr&,
                                          RNG& rng, NeutronEnergy ekin ) const
{
  return m_sh->sampleScatter( rng, ekin );
}

std::shared_ptr<NC::ProcImpl::Process>
NC::SABScatterNG::createMerged( const Process& oraw,
                                double scale_self,
                                double scale_other ) const
{
  nc_assert(scale_other>0.0);
  nc_assert(scale_self>0.0);
  auto optr = dynamic_cast<const SABScatterNG*>(&oraw);
  if (!optr)
    return nullptr;
  auto& o = *optr;

  if ( m_sh != o.m_sh )
    return nullptr;

  const double newscale = scale_self*m_scale + scale_other*o.m_scale;
  return std::make_shared<SABScatterNG>( m_sh, newscale );
}

NC::Optional<std::string> NC::SABScatterNG::specificJSONDescription() const
{
  std::string extension_method
    = m_sh->extender().shortName().value_or("custom");
  std::ostringstream ss;
  m_sh->processor().toJSONProcessInfo( ss,
                                       SigmaBound{ m_scale },
                                       extension_method );
  return ss.str();
}

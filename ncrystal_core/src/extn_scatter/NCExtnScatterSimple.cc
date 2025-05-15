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

#include "NCrystal/internal/extn_scatter/NCExtnScatterSimple.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;
namespace NCE = NCrystal::Extinction;

NCE::ExtnScatterSimple::ExtnScatterSimple( PreparedPowderInputData&& data,
                                           Length domainSize )
{
  if ( data.d_fm_list.empty() )
    return;
  m_threshold = NeutronEnergy{ wl2ekin(2.0*data.d_fm_list.front().first) };
  (void)domainSize;//fixme
}

NC::EnergyDomain NCE::ExtnScatterSimple::domain() const noexcept
{
  return { m_threshold, NeutronEnergy{kInfinity} };
}

NC::CrossSect NCE::ExtnScatterSimple::crossSectionIsotropic( CachePtr&,
                                                             NeutronEnergy ekin ) const
{
  return CrossSect{ ekin < m_threshold ? 0.0 : 1.0  };
}

NC::ScatterOutcomeIsotropic
NCE::ExtnScatterSimple::sampleScatterIsotropic( CachePtr&,
                                                RNG&,
                                                NeutronEnergy ekin ) const
{
  return ScatterOutcomeIsotropic::noScat(ekin);
}

NC::Optional<std::string>
NCE::ExtnScatterSimple::specificJSONDescription() const
{
  //fixme: more stuff
  double max_contrib(0.0);//fixme
  std::ostringstream ss;
  {
    std::ostringstream tmp;
    tmp << "nplanes="<<1000//fixmem_2dE.size()
        <<";2dmax="<<m_threshold.wavelength()
        << ";max_contrib="<<CrossSect{max_contrib};
    streamJSONDictEntry( ss, "summarystr", tmp.str(), JSONDictPos::FIRST );
  }
  streamJSONDictEntry( ss, "nhkl", 1000 );//fixme m_2dE.size() );
  streamJSONDictEntry( ss, "max_contrib", max_contrib );
  streamJSONDictEntry( ss, "2dmax", m_threshold.wavelength().dbl(),
                       JSONDictPos::LAST );
  return ss.str();
}

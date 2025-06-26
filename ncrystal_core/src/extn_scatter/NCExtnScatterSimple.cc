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
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/extn_utils/NCExtnUtils.hh"

namespace NC = NCrystal;
namespace NCE = NCrystal::Extn;

NCE::ExtnScatterSimple::ExtnScatterSimple( PowderBraggInput::Data&& data,
                                           Length domainSize )
  : m_domainSizeAa( domainSize.get()/Length::angstrom )
{
  if ( !std::isfinite(m_domainSizeAa) || !( m_domainSizeAa > 0.0 ) )
    NCRYSTAL_THROW2(BadInput,"Invalid domain Size: "<<domainSize);

  if ( data.planes.empty() )
    return;

  m_threshold = NeutronEnergy{ wl2ekin(2.0*data.planes.front().dsp) };
  m_data = std::move(data);
}

NC::EnergyDomain NCE::ExtnScatterSimple::domain() const noexcept
{
  return { m_threshold, NeutronEnergy{kInfinity} };
}

NC::CrossSect NCE::ExtnScatterSimple::crossSectionIsotropic( CachePtr&,
                                                             NeutronEnergy ekin ) const
{
  if ( ekin < m_threshold )
    return CrossSect{ 0.0 };

  //NB: Sabine mu=y=0 in this simple model.

  auto wl = ekin.wavelength();
  const double v0 = m_data.cell.volume;
  const double kkk = ncsquare(wl.get()*m_domainSizeAa/v0);//(unit is Aa^-2)

  double wlhalf = wl.get()*0.5;
  const double factor = ncsquare(wl.get())/(2.0*m_data.cell.volume*m_data.cell.n_atoms);
  StableSum contrib;
  for ( auto& e : m_data.planes ) {
    if ( e.dsp < wlhalf )
      break;
    double sabine_x = kkk * e.fsq * 1e-8;//1e-8 to convert fsq from barn to to Angstrom^2 (fixme absorb on kkk)

    double El = calcSabineEl_y0( sabine_x );
    double Eb = calcSabineEb_y0( sabine_x );

    double sinth_sq = std::min<double>(1.0,ncsquare(0.5 * wl.get() / e.dsp));
    double costh_sq = std::max<double>(1.0 - sinth_sq,0.0);
    double extinction_correction = El * costh_sq + Eb * sinth_sq;
    contrib.add( e.dsp * e.fsq * e.mult * extinction_correction  );
  }
  return CrossSect{ factor * contrib.sum()  };
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

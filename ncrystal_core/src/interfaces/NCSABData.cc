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

#include "NCrystal/interfaces/NCSABData.hh"
namespace NC = NCrystal;


NC::SABData::SABData( VectD&& alphaGrid,
                      VectD&& betaGrid,
                      VectD&& sab,
                      Temperature temperature,
                      SigmaBound boundXS,
                      AtomMass elementMassAMU,
                      double suggestedEmax )
  : m_a(std::move(alphaGrid)),
    m_b(std::move(betaGrid)),
    m_sab(std::move(sab)),
    m_t(DoValidate,temperature),
    m_m(DoValidate,elementMassAMU),
    m_sem(suggestedEmax),
    m_bxs(DoValidate,boundXS)
{
  nc_assert_always( m_a.size() < std::numeric_limits<std::uint16_t>::max() );
  nc_assert_always( m_b.size() < std::numeric_limits<std::uint16_t>::max() );
  nc_assert(std::is_sorted(m_a.begin(),m_a.end()));//nc_is_grid not available here
  nc_assert(std::is_sorted(m_b.begin(),m_b.end()));//nc_is_grid not available here
  nc_assert(m_a.size()*m_b.size()==m_sab.size());
  nc_assert(m_sab.size()>=4);
  nc_assert(m_sem>=0.0);
  nc_assert(*std::min_element(m_sab.begin(),m_sab.end())>=0.0);
  nc_assert(m_a.front()>=0.0);
  nc_assert(m_b.front()<0.0);
}

NC::VDOSData::VDOSData( PairDD egrid,
                        VectD&& density,
                        Temperature temperature,
                        SigmaBound boundXS,
                        AtomMass elementMassAMU )
  : m_e(egrid),
    m_d(std::move(density)),
    m_t(DoValidate,temperature),
    m_m(DoValidate,elementMassAMU),
    m_bxs(DoValidate,boundXS)
{
  nc_assert( m_e.first>=0.0);
  nc_assert( m_e.first < m_e.second );
  nc_assert( m_d.size()>=2 );
  nc_assert( *std::min_element(m_d.begin(),m_d.end())>=0.0);
}

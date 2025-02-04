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

#include "NCrystal/internal/sab/NCSABExtender.hh"
namespace NC = NCrystal;

NC::SAB::SABExtender::~SABExtender() = default;

NC::SAB::SABFGExtender::SABFGExtender( Temperature temp_k, AtomMass mass, NC::SigmaFree sigma )
  : m_xsprovider(temp_k,mass,sigma),
    m_t(DoValidate,temp_k),
    m_m(DoValidate,mass)
{
}

NC::SAB::SABFGExtender::SABFGExtender( Temperature temp_k, AtomMass mass, NC::SigmaBound sigma )
  : m_xsprovider(temp_k,mass,sigma),
    m_t(DoValidate,temp_k),
    m_m(DoValidate,mass)
{
}

NC::SAB::SABFGExtender::~SABFGExtender() = default;

NC::CrossSect NC::SAB::SABFGExtender::crossSection(NeutronEnergy ekin) const
{
  return m_xsprovider.crossSection(ekin);
}

NC::PairDD NC::SAB::SABFGExtender::sampleAlphaBeta( RNG& rng, NeutronEnergy ekin ) const
{
  return FreeGasSampler(ekin, m_t, m_m).sampleAlphaBeta(rng);
}

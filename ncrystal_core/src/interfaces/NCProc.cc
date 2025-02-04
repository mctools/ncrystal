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

#include "NCrystal/interfaces/NCProc.hh"

namespace NC = NCrystal;

NC::Scatter NC::Scatter::clone()
{
  //Check here that move constructors are noexcept. This is important if
  //e.g. using std::vector<Scatter>:
  static_assert(std::is_nothrow_move_constructible<Scatter>::value, "");
  static_assert(std::is_nothrow_move_constructible<Absorption>::value, "");

  return Scatter( m_rngproducer,
                  m_rngproducer->produce(),
                  m_proc );
}

NC::Scatter NC::Scatter::cloneByIdx( RNGStreamIndex idx )
{
  return Scatter( m_rngproducer,
                  m_rngproducer->produceByIdx(idx),
                  m_proc );
}

NC::Scatter NC::Scatter::cloneForCurrentThread()
{
  return Scatter( m_rngproducer,
                  m_rngproducer->produceForCurrentThread(),
                  m_proc );
}

NC::Scatter NC::Scatter::cloneWithIdenticalRNGSettings()
{
  return Scatter( m_rngproducer, m_rng, m_proc );
}

NC::Scatter NC::Scatter::clone( shared_obj<RNGProducer> rp, shared_obj<RNG> r )
{
  return Scatter( std::move(rp), std::move(r), m_proc );
}

void NC::Scatter::replaceRNG( shared_obj<RNG> r, shared_obj<RNGProducer> rp )
{
  m_rngproducer = std::move(rp);
  m_rng = std::move(r);
}

void NC::Scatter::replaceRNGAndUpdateProducer( shared_obj<RNGStream> r )
{
  //Adopt rng here, and reinit producer from the same rng - but make sure it
  //doesn't actually provide the SAME rng to anyone else:
  m_rngproducer->reinit( r, RNGProducer::SkipOriginal::True );
  m_rng = std::move(r);
}

NC::Absorption NC::Absorption::clone() const
{
  return Absorption( m_proc );
}

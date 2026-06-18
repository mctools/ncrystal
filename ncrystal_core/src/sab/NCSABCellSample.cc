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

#include "NCrystal/internal/sab/NCSABCellSample.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"

namespace NC = NCrystal;
namespace NCS = NCrystal::SABUtils;

NCS::LogLinDistSampler::LogLinDistSampler( double a, double fa, double logfa,
                                           double b, double fb, double logfb,
                                           bool force_linlin )
  : m_islinlin(false), m_a(a), m_b(b)
{
  nc_assert( a >= 0.0 );
  nc_assert( b >= a );
  nc_assert( fa >= 0 );
  nc_assert( fb >= 0 );
  nc_assert( std::isfinite(a) );
  nc_assert( std::isfinite(fa) );
  nc_assert( std::isfinite(b) );
  nc_assert( std::isfinite(fb) );
  nc_assert( std::isfinite(logfa) );
  nc_assert( std::isfinite(logfb) );

  const double bma = b-a;

  if ( !( ncmin(fa,fb,bma) > 0.0 ) )
    force_linlin = true;

  if ( !force_linlin ) {
    const double minus_dlnf = logfa-logfb;
    const double df = fb-fa;

    nc_assert(bma>0.0);
    if ( minus_dlnf == 0 || df == 0 ) {
      //fall back to uniform (via linlin code path):
      fa = fb = 1.0;
      force_linlin = true;
    } else {
      m_cache1 = m_b-m_a;
      m_cache2 = fb/fa;
      m_cache3 = logfb-logfa;
    }
  }

  if ( force_linlin ) {
    m_islinlin = true;
    if ( fa!=fb && bma!=0.0 ) {
      //uniform base + triangle top
      m_cache1 = ncmin(fa,fb);//av. height of uniform base
      m_cache2 = 0.5*(fa+fb);//av. height of base + triangle
      m_cache3 = ( fa>fb ? b : a );
      m_cache4 = ( fa>fb ? -bma : bma );
    } else {
      //uniform
      m_cache1 = 1.0;
      m_cache2 = 1.0;
      m_cache3 = a;
      m_cache4 = bma;
    }
  }
}

NCS::FullCellSampler::FullCellSampler( const CellData* cellptr )
  : m_cellptr(cellptr)
{
  //fixme: alternative constructor which simply takes prob1 or (W1,W2), if we
  //already know it?

  nc_assert(m_cellptr!=nullptr);
  const CellData& c = *m_cellptr;
  const double W1 = integrateAlphaInterval_fast(c.a1,c.S[0],c.a2 , c.S[1],
                                                c.logS[0], c.logS[1]);
  const double W2 = integrateAlphaInterval_fast(c.a1,c.S[2], c.a2 , c.S[3],
                                                c.logS[2], c.logS[3]);
  nc_assert_always(W1>0.0||W2>0.0);
  m_prob1 = W1/(W1+W2);
}

NC::PairDD NCS::FullCellSampler::sampleAlphaBeta( RNG& rng )
{
  nc_assert(m_cellptr!=nullptr);
  const CellData& c = *m_cellptr;
  //pick if we sample from the edge @ beta=c.b1 or the edge @ beta=c.b2:
  const int offset = ( rng.generate() <= m_prob1 ? 0 : 2 );
  auto& opt_sampler = m_samplers[offset/2];
  if ( !opt_sampler.has_value() ) {
    opt_sampler.emplace( c.a1, c.S[offset], c.logS[offset],
                         c.a2, c.S[offset+1], c.logS[offset+1] );
  }
  const double alpha = opt_sampler.value().sample(rng);
  //interpolation along beta is purely linear, so use a triangle which is zero
  //at the far side:
  const double rb = ( offset
                      ? ncmax(rng.generate(),rng.generate())
                      : ncmin(rng.generate(),rng.generate()) );//fixme: sqrt or minmax trick?
  const double beta = c.b1*(1.0-rb) + rb*c.b2;
  return { alpha, ncclamp(beta, c.b1, c.b2) };
}

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
#include "NCrystal/internal/utils/NCMsg.hh"

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

  //interpolation along beta is purely linear, so simply sample a triangle which
  //is zero at the far side:

  const double t = std::sqrt(rng.generate());
  const double beta = ( offset
                        ? ( c.b1*(1.0-t) + t*c.b2 )
                        : ( c.b1*t + (1.0-t)*c.b2 ) );
  return { alpha, ncclamp(beta, c.b1, c.b2) };
}

NCS::BoundedCellSampler::Result
NCS::BoundedCellSampler::sampleAlphaBeta( RNG& rng )
{
  //pick side:
  const unsigned offset = ( rng.generate() <= m_prob1 ? 0 : 2 );
  auto& optEdge = m_edge[offset/2];
  if (!optEdge.has_value()) {
    optEdge.emplace();
    initEdge(optEdge.value(),offset);
  }
  const double overlay_val_orig = m_overlay[offset/2];
  const auto& edge = optEdge.value();
  const auto& c = m_cell;

  LogLinDistSampler sampler( m_restrict_a1, edge.s_low, edge.lns_low,
                             m_restrict_a2, edge.s_up, edge.lns_up,
                             edge.islinlin );

  Result res;
  res.ntries = 0;
  double max_overlay_val_seen = -1.0;
  double overlay_val = overlay_val_orig;
  while ( true ) {
    ++res.ntries;
#if 1
    //FIXME: Just testing a way to limit worst performance. Instead of raising
    //the overlay, we could also simply pick the previous best value after N
    //tries? (where "best" includes both the actual value and the R value used
    //for comparison.
    if ( res.ntries == 200 && max_overlay_val_seen == 0.0 ) {
      static Msg::WarnFirstTime wf;
      wf("S(alpha,beta) cell sampling takes unexpectedly long"
         " to find accepted candidate.");
      overlay_val = 0.0;
    }
    //fixme: just do this at init time (for cells with huge k or 1/k), and lower
    //cached overlay values? That way we can also use more values (e.g. 100).
    if ( res.ntries > 20 && max_overlay_val_seen > 0.0 ) {
      const double new_overlay = (res.ntries>20?1.2:2.0)*max_overlay_val_seen;
      overlay_val = ncmin(overlay_val_orig,new_overlay);
    }
    overlay_val = ncmax(overlay_val,1.001*max_overlay_val_seen);
#endif

    nc_assert_always( res.ntries < 10000ull );
    double a = sampler.sample(rng);
    const double dpmb = std::sqrt(m_4e*a);
    double bl = a-dpmb;
    double bu = a+dpmb;
    double bwidth, bmid;
    if ( valueInInterval(c.b1,c.b2,bl)
         && valueInInterval(c.b1,c.b2,bu) ) {
      bwidth = 2.0*dpmb;
      bmid = a;
    } else {
      bl = ncmax(c.b1,bl);
      bu = ncmin(c.b2,bu);
      bwidth = bu-bl;
      bmid = (bl+bu)*0.5;
    }
    if ( !(bwidth>0) )
      continue;
    //Now, the relative beta weight starts at 1.0 at the chosen side,
    //and reduces to 0.0 at the opposite side (due to linear beta
    //interpolation).
    double bmidheight_multdb = ( offset ? (bmid-c.b1) : (c.b2-bmid ) );
    const double actual_val = bwidth*bmidheight_multdb;
    max_overlay_val_seen = ncmax(max_overlay_val_seen,actual_val);
    if ( rng.generate() * overlay_val > actual_val )
      continue;//reject

    //Ok, sample beta.
    const double hbl = ( offset ? (bl-c.b1) : (c.b2-bl ) );
    const double hbu = ( offset ? (bu-c.b1) : (c.b2-bu ) );
    const double hbase = ncmin(hbl,hbu);
    const double havg = (hbl+hbu)*0.5;
    double u;
    if ( rng.generate()*havg <= hbase ) {
      //base
      u = rng.generate();
    } else {
      //triangle (with highest point in appropriate direction)
      u = std::sqrt(rng.generate());
      if ( !offset )
        u = 1.0 - u;
    }
    const double b = ncmin(bu,bl+u*bwidth);
    nc_assert( ncsquare(a-b) <= a*m_4e*(1.0+1e-10) );
    res.alpha = a;
    res.beta = ncclamp( b, c.b1, c.b2 );
    return res;
  }
}

void NCS::BoundedCellSampler::initEdge( EdgeData& edge, unsigned offset ) const
{
  const auto& c = m_cell;
  edge.s_low = c.S[offset];
  edge.s_up = c.S[offset+1];
  edge.islinlin = !(ncmin( edge.s_low, edge.s_up )>0.0);
  if ( edge.islinlin ) {
    edge.lns_low = edge.lns_up = 0.0;
    if ( m_restrict_a1 > m_cell.a1 )
      edge.s_low = interpolate_linlin_NEW(m_cell.a1, c.S[offset],
                                          m_cell.a2, c.S[offset+1],
                                          m_restrict_a1);
    if ( m_restrict_a2 < m_cell.a2 )
      edge.s_up = interpolate_linlin_NEW(m_cell.a1, c.S[offset],
                                         m_cell.a2, c.S[offset+1],
                                         m_restrict_a2);
  } else {
    edge.lns_low = c.logS[offset];
    edge.lns_up = c.logS[offset+1];
    if ( m_restrict_a1 > m_cell.a1 ) {
      auto rs = interpolate_loglin_fast2_NEW(m_cell.a1, c.S[offset],
                                             m_cell.a2, c.S[offset+1],
                                             m_restrict_a1,
                                             c.logS[offset],
                                             c.logS[offset+1]);
      edge.s_low = rs.first;
      edge.lns_low = rs.second;
    }

    if ( m_restrict_a2 < m_cell.a2 ) {
      auto rs = interpolate_loglin_fast2_NEW(m_cell.a1, c.S[offset],
                                             m_cell.a2, c.S[offset+1],
                                             m_restrict_a2,
                                             c.logS[offset],
                                             c.logS[offset+1]);
      edge.s_up = rs.first;
      edge.lns_up = rs.second;
    }
  }
}

NCS::BoundedCellSampler::ProcessedSurveyInfo
NCS::BoundedCellSampler::processSurveyInfo( const SABCellSurvey& survey,
                                            double alpha1, double alpha2,
                                            double beta1, double beta2,
                                            double E_div_kT )
{

  ProcessedSurveyInfo res;
  res.overlay2 = -1.0f;
  res.overlay1 = -1.0f;
  res.alpha_up = -1.0f;
  res.alpha_low = -1.0f;

  auto& regions = survey.regions();
  if ( regions.empty() )
    return res;

  const double b1 = beta1;
  const double b2 = beta2;
  const double e = E_div_kT;

  double o1(-1.0), o2(-1.0);
  auto it = regions.begin();
  auto itE = regions.end();
  SABCellSurvey::Region current = *it++;
  auto updateO1 = [&o1,b2]( double bmid, double bwidth )
  {
    o1 = ncmax(o1, (b2-bmid)*bwidth );

  };
  auto updateO2 = [&o2,b1]( double bmid, double bwidth )
  {
    o2 = ncmax(o2, (bmid-b1)*bwidth );

  };
  auto updateO12 = [&updateO1,&updateO2]( double bmid, double bwidth )
  {
    updateO1(bmid,bwidth);
    updateO2(bmid,bwidth);
  };
  auto processCurrent = [b1,b2,e,&current,
                         &updateO1,&updateO2,&updateO12]()
  {
    const auto& r = current;
    if ( r.is_bounded_by_betaplus ) {
      if ( r.is_bounded_by_betaminus ) {
        //bound by [betaminus(alpha),betaplus(alpha)]
        updateO12( r.alpha_up, 4.0*std::sqrt(e * r.alpha_up) );
        updateO12( r.alpha_low, 4.0*std::sqrt(e * r.alpha_low) );
        //local maximum at a=b2/3 for o1 and at a=b1/3 for o2:
        constexpr double onethird = 1.0/3.0;
        const double amax1 = b2*onethird;
        const double amax2 = b1*onethird;
        if ( valueInInterval( r.alpha_low, r.alpha_up, amax1 ) )
          updateO1( amax1, 4.0*std::sqrt(e * amax1) );
        if ( valueInInterval( r.alpha_low, r.alpha_up, amax2 ) )
          updateO2( amax2, 4.0*std::sqrt(e * amax2) );
      } else {
        //bound by [b1,betaplus(alpha)]
        const double sqrte = std::sqrt(e);
        const double twosqrte = 2.0 * sqrte;
        auto pt = [b1,updateO1,updateO2,twosqrte]( double a,
                                                   bool do1 = true,
                                                   bool do2 = true ) {
          const double bplus = a + twosqrte*std::sqrt(a);
          const double bmid((bplus+b1)*0.5), bwidth(bplus-b1);
          if (do1)
            updateO1(bmid,bwidth);
          if (do2)
            updateO2(bmid,bwidth);
        };
        pt(r.alpha_up);
        pt(r.alpha_low);
        const double amax1 = ncsquare( std::sqrt(e+b2)- sqrte );
        const double amax2 = ncsquare( std::sqrt(e+b1)- sqrte );
        if ( valueInInterval( r.alpha_low, r.alpha_up, amax1 ) )
          pt( amax1, true, false );
        if ( valueInInterval( r.alpha_low, r.alpha_up, amax2 ) )
          pt( amax2, false, true );
      }
    } else {
      if ( r.is_bounded_by_betaminus ) {
        //bound by [betaminus(alpha),b2]
        const double sqrte = std::sqrt(e);
        const double twosqrte = 2.0 * sqrte;
        auto pt = [b2,updateO12,twosqrte]( double a ) {
          const double bminus = a - twosqrte*std::sqrt(a);
          updateO12((bminus+b2)*0.5,b2-bminus);
        };
        pt( r.alpha_up );
        pt( r.alpha_low );
        if ( valueInInterval( r.alpha_low, r.alpha_up, e ) )
          pt( e );
      } else {
        //bound by [b1,b2]
        updateO12( (b1+b2)*0.5, b2-b1 );
      }
    }
  };

  for ( ; it!=itE; ++it ) {
    if ( current.is_bounded_by_betaminus != it->is_bounded_by_betaminus
         || current.is_bounded_by_betaplus != it->is_bounded_by_betaplus ) {
      processCurrent();
      current = *it;
    } else {
      //merge previously split regions:
      assert( current.alpha_low == it->alpha_up );
      current.alpha_low = it->alpha_low;
    }
  }
  processCurrent();

  //finish up, rounding single precision floats the right way and add a bit of
  //safety:
  auto storeFloat = [](float&dest, double val, bool push_is_up=true)
  {
    dest = static_cast<float>( val );
    if ( push_is_up ) {
      while ( dest < val )
        dest = std::nextafter(dest,std::numeric_limits<float>::infinity());
    } else {
      while ( dest > val )
        dest = std::nextafter(dest,-std::numeric_limits<float>::infinity());
    }
  };

  if ( regions.front().alpha_up < alpha2 ) {
    storeFloat(res.alpha_up,regions.front().alpha_up);
    if ( !(static_cast<double>(res.alpha_up)<alpha2) )
      res.alpha_up = -1.0f;
  }

  if ( regions.back().alpha_low > alpha1 ) {
    storeFloat(res.alpha_low,regions.back().alpha_low,false);
    if ( !(static_cast<double>(res.alpha_low)>alpha1) )
      res.alpha_low = -1.0f;
  }

  storeFloat(res.overlay1,o1);
  storeFloat(res.overlay2,o2);
  nc_assert( res.overlay1 >= 0.0f );
  nc_assert( res.overlay2 >= 0.0f );

  constexpr double overlay_safety_factor = 1.0001;
  res.overlay1 *= overlay_safety_factor;
  res.overlay2 *= overlay_safety_factor;
  nc_assert( res.overlay1 > 0.0f );
  nc_assert( res.overlay2 > 0.0f );

  return res;
}

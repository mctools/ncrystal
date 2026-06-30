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
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;
namespace NCS = NCrystal::SABUtils;

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {
    namespace LargeKHandling {
      //For very small k=S2/S1 (or S1/S2, whichever is smaller), alpha sampling
      //will in practice only pick a relative alpha in [0,u] for some small
      //u. To avoid bad acceptance rates later, we will in the
      //BoundedCellSampler use this information to limit the alpha ranges and
      //overlay values to the most significant part of the alpha
      //range. Otherwise the overlay value can end up drastically overshooting
      //the value in the region which is actually sampled.
      //
      // If we define a PDF f(x) = norm*k^t on [0,1], then norm=log(k)/(k-1) (we
      // only deal with small k, so k!=1 is a safe assumption).
      //
      //For a given k<1 we want to determine u so that integral(f,x=0..u)=1-eps.
      //This is equivalent to the condition:
      //
      //       k^u = 1-(1-eps)*(1-k)
      //  <=>  u = ln(1-(1-eps)*(1-k))/ln(k)
      //
      //To avoid the special treatment when not needed, we will only consider it
      //when u<0.5, leading to the threshold value of k given by:
      //     k^0.5 = 1-(1-eps)*(1-k) => k = (eps/(1-eps))^2 ~= eps^2
      //
      //So if for instance eps=1e-4 we only apply the procedure when S2/S1 is
      //outside the interval [1e-8,1e8].
      constexpr double eps = 1e-6;//fixme "luxury"!!
      constexpr double oneminuseps = 1.0 - eps;
      constexpr double threshold_k = eps*eps/(oneminuseps*oneminuseps);
      constexpr double inv_threshold_k = 1.0/threshold_k;
      bool needsLargeKTreatment( double s1, double s2 )
      {
        return ( s1>0.0 && s2>0.0
                 && !valueInInterval( threshold_k*s2,
                                      inv_threshold_k*s2, s1 ) );
      }

      PairDD findRestrictedAlphaRange( double a1, double a2,
                                       double s1, double s2,
                                       double lns1, double lns2 )
      {
        nc_assert( needsLargeKTreatment(s1,s2) );
        nc_assert( floateq( std::log(s1), lns1 ) );
        nc_assert( floateq( std::log(s2), lns2 ) );
        double k,lnk;
        nc_assert( s2>s1 || s1>s2 );
        nc_assert( s1>0.0 && s2>0.0 );
        if ( s2>s1 ) {
          k = s1/s2;
          lnk = lns1-lns2;
        } else {
          k = s2/s1;
          lnk = lns2-lns1;
        }
        nc_assert( k < 1.0 );
        nc_assert( lnk < 0.0 );
        const double t = std::log1p( oneminuseps * ( k - 1.0 ) ) / lnk;
        if ( s2>s1 ) {
          return { intervalPos( a1, a2, 1.0-t ), a2 };
        } else {
          return { a1, intervalPos(a1,a2,t) };
        }
      }

    }
  }
}

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

  LogLinDistSampler sampler( m_restrict_a1[offset/2], edge.s_low, edge.lns_low,
                             m_restrict_a2[offset/2], edge.s_up, edge.lns_up,
                             edge.islinlin );

  Result res;
  res.ntries = 0;
  double max_overlay_val_seen = -1.0;
  double overlay_val = overlay_val_orig;
  while ( true ) {
    ++res.ntries;
    nc_assert_always(res.ntries< 200);//fixme _always + we should be able to
                                      //lower this!
#if 0
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
    if ( res.ntries > 50 && max_overlay_val_seen > 0.0 ) {
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
    if ( m_restrict_a1[offset/2] > m_cell.a1 )
      edge.s_low = interpolate_linlin_NEW(m_cell.a1, c.S[offset],
                                          m_cell.a2, c.S[offset+1],
                                          m_restrict_a1[offset/2]);
    if ( m_restrict_a2[offset/2] < m_cell.a2 )
      edge.s_up = interpolate_linlin_NEW(m_cell.a1, c.S[offset],
                                         m_cell.a2, c.S[offset+1],
                                         m_restrict_a2[offset/2]);
  } else {
    edge.lns_low = c.logS[offset];
    edge.lns_up = c.logS[offset+1];
    if ( m_restrict_a1[offset/2] > m_cell.a1 ) {
      auto rs = interpolate_loglin_fast2_NEW(m_cell.a1, c.S[offset],
                                             m_cell.a2, c.S[offset+1],
                                             m_restrict_a1[offset/2],
                                             c.logS[offset],
                                             c.logS[offset+1]);
      edge.s_low = rs.first;
      edge.lns_low = rs.second;
    }

    if ( m_restrict_a2[offset/2] < m_cell.a2 ) {
      auto rs = interpolate_loglin_fast2_NEW(m_cell.a1, c.S[offset],
                                             m_cell.a2, c.S[offset+1],
                                             m_restrict_a2[offset/2],
                                             c.logS[offset],
                                             c.logS[offset+1]);
      edge.s_up = rs.first;
      edge.lns_up = rs.second;
    }
  }
}

NCS::BoundedCellSampler::BCSData
NCS::BoundedCellSampler::prepareBCSData( double probability_b1_edge,
                                         const SABCellSurvey& survey,
                                         const CellData& cell,
                                         double E_div_kT )
{

  BCSData res;
  res.probability_b1_edge = probability_b1_edge;
  res.overlay[0] = res.overlay[1] = -1.0;
  res.alpha_low[0] = res.alpha_low[1] = -1.0;
  res.alpha_up[0] = res.alpha_up[1] = -1.0;

  auto& regions = survey.regions();
  if ( regions.empty() )
    return res;

  //large k treatment? First do a check across the cell, then a more specific
  //one in the actual alpha-range of the active regions.
  Optional<PairDD> restrict_a[2];
  bool needs_largek[2] = {
    LargeKHandling::needsLargeKTreatment( cell.S[0], cell.S[1] ),
    LargeKHandling::needsLargeKTreatment( cell.S[2], cell.S[3] )
  };

  if ( needs_largek[0] || needs_largek[1] ) {
    const bool region_restricts_a1 = ( regions.back().alpha_low > cell.a1 );
    const bool region_restricts_a2 = ( regions.front().alpha_up < cell.a2 );
    double S[4] = { cell.S[0], cell.S[1], cell.S[2], cell.S[3] };
    double logS[4] = { cell.logS[0], cell.logS[1], cell.logS[2], cell.logS[3] };
    double a1 = cell.a1;
    double a2 = cell.a2;
    //find restricted S/logS values if relevant:
    auto s_interp = [&cell](double a,int offset) {
      return interpolate_loglin_fast2_NEW( cell.a1, cell.S[offset],
                                           cell.a2, cell.S[offset+1],
                                           a,
                                           cell.logS[offset], cell.logS[offset+1] );
    };
    if ( region_restricts_a1 ) {
      auto sls1 = s_interp(regions.back().alpha_low,0);
      auto sls2 = s_interp(regions.back().alpha_low,2);
      S[0] = sls1.first;
      logS[0] = sls1.second;
      S[2] = sls2.first;
      logS[2] = sls2.second;
      a1 = regions.back().alpha_low;
    }
    if ( region_restricts_a2 ) {
      auto sls1 = s_interp(regions.front().alpha_up,0);
      auto sls2 = s_interp(regions.front().alpha_up,2);
      S[1] = sls1.first;
      logS[1] = sls1.second;
      S[3] = sls2.first;
      logS[3] = sls2.second;
      a2 = regions.front().alpha_up;
    }
    //recheck if still needs largek treatment in this reduced range:
    needs_largek[0] = LargeKHandling::needsLargeKTreatment( S[0], S[1] );
    needs_largek[1] = LargeKHandling::needsLargeKTreatment( S[2], S[3] );
    //Find the reduced ranges if needed:
    for ( int i = 0; i < 2; ++i ) {
      auto offset = 2*i;
      if (needs_largek[i]) {
        auto rarange
          = LargeKHandling::findRestrictedAlphaRange( a1, a2,
                                                      S[offset], S[offset+1],
                                                      logS[offset],
                                                      logS[offset+1]);
        restrict_a[i].emplace(rarange);
      }
    }
  }

  const bool restrict_any = ( restrict_a[0].has_value()
                              || restrict_a[1].has_value() );
  const double b1 = cell.b1;
  const double b2 = cell.b2;
  const double e = E_div_kT;

  double o1(-1.0), o2(-1.0);
  auto it = regions.begin();
  auto itE = regions.end();
  SABCellSurvey::Region current = *it++;
  auto updateO1 = [&o1,b2,&restrict_a]( double bmid, double bwidth, double a )
  {
    if ( a>=0.0 && restrict_a[0].has_value() ) {
      if ( a < restrict_a[0].value().first || a > restrict_a[0].value().second )
        return;
    }
    o1 = ncmax(o1, (b2-bmid)*bwidth );
  };
  auto updateO2 = [&o2,b1,&restrict_a]( double bmid, double bwidth, double a )
  {
    if ( a>=0.0 && restrict_a[1].has_value() ) {
      if ( a < restrict_a[1].value().first || a > restrict_a[1].value().second )
        return;
    }
    o2 = ncmax(o2, (bmid-b1)*bwidth );
  };
  auto updateO12 = [&updateO1,&updateO2]( double bmid, double bwidth, double a )
  {
    updateO1(bmid,bwidth,a);
    updateO2(bmid,bwidth,a);
  };

  SmallVector<std::pair<int,double>,4> restriction_pts;//{ sideidx, alpha }
  for ( int i = 0; i < 2; ++i ) {
    if ( restrict_a[i].has_value() ) {
      restriction_pts.emplace_back( i, restrict_a[i].value().first );
      restriction_pts.emplace_back( i, restrict_a[i].value().second );
    }
  }

  auto processCurrent = [b1,b2,e,&current,
                         &restrict_a,restrict_any,&restriction_pts,
                         &updateO1,&updateO2,&updateO12]()
  {
    const auto& r = current;

    SmallVector<std::pair<int,double>,4> restriction_pts_in_region;
    for ( auto& i_a : restriction_pts ) {
      if ( valueInInterval( r.alpha_low, r.alpha_up, i_a.second ) )
        restriction_pts_in_region.emplace_back( i_a );
    }

    if ( r.is_bounded_by_betaplus ) {
      if ( r.is_bounded_by_betaminus ) {
        //bound by [betaminus(alpha),betaplus(alpha)]
        updateO12( r.alpha_up, 4.0*std::sqrt(e * r.alpha_up), r.alpha_up );
        updateO12( r.alpha_low, 4.0*std::sqrt(e * r.alpha_low), r.alpha_low );
        //local maximum at a=b2/3 for o1 and at a=b1/3 for o2:
        constexpr double onethird = 1.0/3.0;
        const double amax1 = b2*onethird;
        const double amax2 = b1*onethird;
        if ( valueInInterval( r.alpha_low, r.alpha_up, amax1 ) )
          updateO1( amax1, 4.0*std::sqrt(e * amax1), amax1 );
        if ( valueInInterval( r.alpha_low, r.alpha_up, amax2 ) )
          updateO2( amax2, 4.0*std::sqrt(e * amax2), amax2 );
        //check for restriction points:
        for ( auto& i_a : restriction_pts_in_region ) {
          const double bw = 4.0*std::sqrt(e * i_a.second);
          ( i_a.first
            ? updateO2( i_a.second, bw, i_a.second )
            :  updateO1( i_a.second, bw, i_a.second ) );
        }
      } else {
        //bound by [b1,betaplus(alpha)]
        const double sqrte = std::sqrt(e);
        const double twosqrte = 2.0 * sqrte;
        auto pt = [b1,updateO1,updateO2,twosqrte]( double a,
                                                   bool do1,
                                                   bool do2 ) {
          const double bplus = a + twosqrte*std::sqrt(a);
          const double bmid((bplus+b1)*0.5), bwidth(bplus-b1);
          if (do1)
            updateO1(bmid,bwidth,a);
          if (do2)
            updateO2(bmid,bwidth,a);
        };
        pt(r.alpha_up,true,true);
        pt(r.alpha_low,true,true);
        const double amax1 = ncsquare( std::sqrt(e+b2)- sqrte );
        const double amax2 = ncsquare( std::sqrt(e+b1)- sqrte );
        if ( valueInInterval( r.alpha_low, r.alpha_up, amax1 ) )
          pt( amax1, true, false );
        if ( valueInInterval( r.alpha_low, r.alpha_up, amax2 ) )
          pt( amax2, false, true );
        //check for restriction points:
        for ( auto& i_a : restriction_pts_in_region ) {
          i_a.first ? pt(i_a.second,true,false) : pt(i_a.second,false,true);
        }
      }
    } else {
      if ( r.is_bounded_by_betaminus ) {
        //bound by [betaminus(alpha),b2]
        const double sqrte = std::sqrt(e);
        const double twosqrte = 2.0 * sqrte;
        auto pt = [b2,updateO12,twosqrte]( double a ) {
          const double bminus = a - twosqrte*std::sqrt(a);
          updateO12((bminus+b2)*0.5,b2-bminus,a);
        };
        pt( r.alpha_up );
        pt( r.alpha_low );
        if ( valueInInterval( r.alpha_low, r.alpha_up, e ) )
          pt( e );
        //check for restriction points:
        for ( auto& i_a : restriction_pts_in_region )
          pt( i_a.second );//potentially calling twice with same point, but not
                           //a problem and extremely rare.
      } else {
        //bound by [b1,b2]
        const double bm = (b1+b2)*0.5;
        const double bw = b2-b1;
        if (!restrict_any) {
          updateO12( bm, bw, -1.0 );
        } else {
          if ( !restrict_a[0].has_value()
               || intervalsOverlap( r.alpha_low, r.alpha_up,
                                    restrict_a[0].value().first,
                                    restrict_a[0].value().second ) ) {
            updateO1( bm, bw, -1.0 );
          }
          if ( !restrict_a[1].has_value()
               || intervalsOverlap( r.alpha_low, r.alpha_up,
                                    restrict_a[1].value().first,
                                    restrict_a[1].value().second ) ) {
            updateO2( bm, bw, -1.0 );
          }
        }
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
  auto storeFloat = [](float&dest, double val, bool push_is_up)
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

#if 0
  if ( regions.front().alpha_up < cell.a2 ) {
    storeFloat(res.alpha_up[0],regions.front().alpha_up,true);
    if ( !(static_cast<double>(res.alpha_up[0])<cell.a2) )
      res.alpha_up[0] = -1.0f;
    res.alpha_up[1] = res.alpha_up[0];//FIXME JUST STORING ONE!
  }

  if ( regions.back().alpha_low > cell.a1 ) {
    storeFloat(res.alpha_low[0],regions.back().alpha_low,false);
    if ( !(static_cast<double>(res.alpha_low[0])>cell.a1) )
      res.alpha_low[0] = -1.0f;
    res.alpha_low[1] = res.alpha_low[0];//FIXME JUST STORING ONE!
  }
#else
  for ( int i = 0; i < 2; ++i ) {
    double aup = ncmin( regions.front().alpha_up, cell.a2 );
    double alow = ncmax( regions.back().alpha_low, cell.a1 );
    if ( restrict_a[i].has_value() ) {
      aup = ncmin( aup, restrict_a[i].value().second );
      alow = ncmax( alow, restrict_a[i].value().first );
    }
    storeFloat(res.alpha_up[i],aup,true);
    if ( !(static_cast<double>(res.alpha_up[i])<cell.a2) )
      res.alpha_up[i] = -1.0f;
    storeFloat(res.alpha_low[i],alow,false);
    if ( !(static_cast<double>(res.alpha_low[i])>cell.a1) )
      res.alpha_low[i] = -1.0f;
  }
#endif

  storeFloat(res.overlay[0],o1,true);
  storeFloat(res.overlay[1],o2,true);
  nc_assert( res.overlay[0] >= 0.0f );
  nc_assert( res.overlay[1] >= 0.0f );

  constexpr double overlay_safety_factor = 1.0001;
  res.overlay[0] *= overlay_safety_factor;
  res.overlay[1] *= overlay_safety_factor;
  nc_assert( res.overlay[0] > 0.0f );
  nc_assert( res.overlay[1] > 0.0f );

  return res;
}

std::size_t NCS::BoundedCellSampler::pack( MixedDataVector& buffer, const BCSData& data )
{
  nc_assert( data.overlay[0]>0.0 );
  nc_assert( data.overlay[1]>0.0 );
  nc_assert( data.probability_b1_edge >= 0.0 );
  nc_assert( data.probability_b1_edge <= 1.0 );
  std::size_t idx = buffer.append<double>( data.probability_b1_edge );

  //Figure out what kind of alpha range limits are available:
  const bool has_alow = ncmax( data.alpha_low[0],
                               data.alpha_low[1]) >= 0.0;
  const bool has_aup = ncmax( data.alpha_up[0],
                              data.alpha_up[1]) >= 0.0;
  auto packalpha = [&buffer]( const float (&a)[2] )
  {
    nc_assert( a[0]>0.0 || a[1] > 0.0 );
    nc_assert( a[0]>0.0 || a[0] == -1.0 );
    nc_assert( a[1]>0.0 || a[1] == -1.0 );
    //encode 2 alpha values. If identical, just put the single
    //value. Otherwise put the two values with negative sign (one might be
    //-1.0 if only one is set).
    if ( a[0]==a[1] ) {
      buffer.append<float>( a[0] );
    } else {
      buffer.append<float>( -a[0] );
      buffer.append<float>( -a[1] );
    }
  };
  if ( has_alow ) {
    if ( has_aup ) {
      //both alow+aup overrides ( signature: -ol[0] -ol[1] ).
      buffer.append<float>( -data.overlay[0] );
      buffer.append<float>( -data.overlay[1] );
      packalpha( data.alpha_low );
      packalpha( data.alpha_up );
    } else {
      //just alow overrides ( signature: -ol[0] +ol[1] ).
      buffer.append<float>( -data.overlay[0] );
      buffer.append<float>( data.overlay[1] );
      packalpha( data.alpha_low );
    }
  } else {
    if ( has_aup ) {
      //just aup overrides ( signature: +ol[0] -ol[1] ).
      buffer.append<float>( data.overlay[0] );
      buffer.append<float>( -data.overlay[1] );
      packalpha( data.alpha_up );
    } else {
      //neither alow or aup overrides ( signature: +ol[0] +ol[1] ).
      buffer.append<float>( data.overlay[0] );
      buffer.append<float>( data.overlay[1] );
    }
  }
  return idx;
}

NCS::BoundedCellSampler::BCSData
NCS::BoundedCellSampler::unpack( const MixedDataVector& buffer,
                                 std::size_t index )
{
  BCSData data;
  data.probability_b1_edge = buffer.extract<double>(index);
  index += sizeof(double);
  auto getfloat = [&buffer,&index]()
  {
    float v = buffer.extract<float>(index);
    index += sizeof(float);
    return v;
  };
  float ol0 = getfloat();
  float ol1 = getfloat();
  data.overlay[0] = ncabs(ol0);
  data.overlay[1] = ncabs(ol1);
  auto unpackalpha = [&getfloat]( float (&a)[2] )
  {
    a[0] = getfloat();
    if ( a[0] > 0.0 ) {
      a[1] = a[0];
    } else {
      a[0] = -a[0];
      a[1] = -getfloat();
    }
  };

  if ( ol0 > 0.0 ) {
    if ( ol1 > 0.0 ) {
      //neither alow or aup overrides ( signature: +ol[0] +ol[1] ).
      data.alpha_low[0] = data.alpha_low[1] = -1.0;
      data.alpha_up[0] = data.alpha_up[1] = -1.0;
    } else {
      //just aup overrides ( signature: +ol[0] -ol[1] ).
      unpackalpha(data.alpha_up);
      data.alpha_low[0] = data.alpha_low[1] = -1.0;
    }
  } else {
    if ( ol1 > 0.0 ) {
      //just alow overrides ( signature: -ol[0] +ol[1] ).
      unpackalpha(data.alpha_low);
      data.alpha_up[0] = data.alpha_up[1] = -1.0;
    } else {
      //both alow+aup overrides ( signature: -ol[0] -ol[1] ).
      unpackalpha(data.alpha_low);
      unpackalpha(data.alpha_up);
    }
  }
  return data;
}

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {
    namespace {
      Optional<PairDD> findActiveAlphaRange( const CellData& c, double E_div_kT )
      {
        double a1 = c.a1;
        double a2 = c.a2;
        const double b1 = c.b1;
        const double b2 = c.b2;
        const double e = E_div_kT;
        if (c.b2 <= -e)
          return NullOpt;//no overlap
        const double twoe = 2.0*e;
        const double tmp = 2.0*std::sqrt(e*(b2+e));
        const double ap2 = twoe+b2+tmp;//alpha^+(b2)
        if ( a1 >= ap2 )
          return NullOpt;//no overlap
        const double am2 = twoe+b2-tmp;//alpha^-(b2)
        double am1(-1.0);
        if ( b1 >= -e ) {
          const double tmp2 = 2.0*std::sqrt(e*(b1+e));
          const double twoe_plus_b1 = twoe+b1;
          am1 = twoe_plus_b1-tmp2;//alpha^-(b1)
        }
        //snap along alpha:
        a2 = ncmin(a2,ap2);
        if ( b2 < 0 ) {
          a1 = ncmax(a1,am2);
        } else if ( b1 > 0 ) {
          a1 = ncmax(a1,am1);
        }
        if (!(a2 > a1))
          return NullOpt;//no overlap
        return PairDD{a1,a2};
      }


      class RCSImpl final : private NoCopyMove {
        CellData m_cell;
        double m_e;//E/kT
        double m_alow = -1.0, m_aup=-1.0;//active alpha range.
        struct AlphaPtGeom {
          double blow, bup, bmid, bwidth;
        };
        VectD m_abinedges;
        VectD m_bincontrib_commul;
        VectD m_bincontrib;

        AlphaPtGeom alphaPtGem( double a ) const
        {
          nc_assert(m_aup>m_alow);
          nc_assert(valueInInterval(m_alow,m_aup,a));
          AlphaPtGeom res;
          double tmp = 2.0*std::sqrt( a * m_e);
          double bl = a - tmp;
          double bu = a + tmp;
          res.blow = ncmax(m_cell.b1,bl);
          res.bup = ncmin(m_cell.b2,bu);
          if ( bl >= m_cell.b1 && bu <= m_cell.b2 ) {
            res.bmid = a;
            res.bwidth = 2.0*tmp;
          } else {
            res.bmid = 0.5*(res.blow+res.bup);
            res.bwidth = res.bup - res.blow;
          }
          return res;
        }

        PairDD findSValsAtA( double a ) const
        {
          auto& S = m_cell.S;
          auto& logS = m_cell.logS;
          double s1 = interpolate_loglin_fallbacklinlin_fast( m_cell.a1, S[0],
                                                              m_cell.a2, S[1],
                                                              a,
                                                              logS[0], logS[1]);
          double s2 = interpolate_loglin_fallbacklinlin_fast( m_cell.a1, S[2],
                                                              m_cell.a2, S[3],
                                                              a,
                                                              logS[2], logS[3]);
          return {s1,s2};
        }

        double contribOfA( double a ) const
        {
          //Contribution is integral over beta at S(alpha=a,beta) multiplied
          //with (b2-b1), since it is a common factor everywhere and we are only
          //looking at sampling anyway.
          const AlphaPtGeom pt = alphaPtGem( a );
          const PairDD svals = findSValsAtA( a );
          StableSum ss, ss2;
          ss.add(m_cell.b2);
          ss.add(-pt.bmid);
          ss.mult(svals.first);
          ss2.add(pt.bmid);
          ss2.add(-m_cell.b1);
          ss2.mult(svals.second);
          ss.add(ss2);
          ss.mult(pt.bwidth);
          return ss.sum();
        }

      public:
        RCSImpl( const CellData& cell, double E_div_kT )
          : m_cell(cell), m_e(E_div_kT)
        {
          {
            auto arange = findActiveAlphaRange( cell, E_div_kT );
            if (!arange.has_value())
              NCRYSTAL_THROW(BadInput,"RefCellSampler can not be used with"
                             " cell with no phasespace overlap");
            m_alow = arange.value().first;
            m_aup = arange.value().second;
          }
          nc_assert_always(m_alow < m_aup );
          nc_assert_always(m_alow >= cell.a1);
          nc_assert_always(m_aup <= cell.a2);

          //The contribution(alpha) might have local maxima at a few special
          //values (found in connection with the BoundedCellSampler), so for
          //added safety, we must also evaluate on these if they fall in a given
          //alpha-bin. Otherwise the overlay value determined purely from the
          //alpha-bin edges could be underestimated.
          SmallVector<double,5> special_avals;
          {
            //candidates (the ncabs(..) inside the sqrt is inserted for safety,
            //it is no harm to add spurious candidates here if that candidate
            //did not actually exist in the given setup.
            const double sqrte = std::sqrt(m_e);
            const double da1 = 2*std::sqrt(ncabs(m_e*(m_cell.b1+m_e)));
            const double da2 = 2*std::sqrt(ncabs(m_e*(m_cell.b2+m_e)));

            SmallVector<double,5> special_avals_candidates
              = { m_e, m_cell.b1/3, m_cell.b2/3,
                  ncsquare( std::sqrt( ncabs(m_e + m_cell.b1) ) - sqrte ),
                  ncsquare( std::sqrt( ncabs(m_e + m_cell.b2) ) - sqrte ),
                  2*m_e+m_cell.b1-da1,
                  2*m_e+m_cell.b1+da1,
                  2*m_e+m_cell.b2-da2,
                  2*m_e+m_cell.b2+da2 };
            for ( auto a : special_avals_candidates ) {
              if ( valueInInterval(m_alow,m_aup,a) )
                special_avals.push_back(a);
            }
          }

          //Now set up a number of alpha bins in which to sample:
          {
            constexpr std::size_t nbins = 1024;
            m_abinedges = linspace( m_alow, m_aup, nbins+1 );
            m_bincontrib_commul.reserve(nbins);
            StableSum contrib_commul_ssum;
            double bin_alow = m_abinedges.front();
            double contrib_low_edge = contribOfA(bin_alow);
            for ( auto ibin : ncrange(nbins) ) {
              double bin_aup = vectAt(m_abinedges,ibin+1);
              double contrib_high_edge = contribOfA(bin_aup);
              double contrib = ncmax( contrib_low_edge, contrib_high_edge );
              for ( auto aspecial : special_avals ) {
                if ( valueInInterval( bin_alow, bin_aup, aspecial ) )
                  contrib = ncmax( contrib, contribOfA(aspecial));
              }
              contrib_commul_ssum.add( contrib );
              m_bincontrib_commul.push_back(contrib_commul_ssum.sum());
              m_bincontrib.push_back(contrib);
              contrib_low_edge = contrib_high_edge;
              bin_alow = bin_aup;
            }
          }
        }

        RefCellSampler::Result sample( RNG& rng )
        {
          RefCellSampler::Result res;
          res.ntries = 0;
          //Sample alpha value first:
          double a(-1.0);
          while (true) {
            ++res.ntries;
            nc_assert_always(res.ntries<1000);
            auto ibin = pickRandIdxByWeight( rng, m_bincontrib_commul );
            const double overlay = vectAt( m_bincontrib, ibin );
            a = randInterval( rng,
                              vectAt( m_abinedges, ibin ),
                              vectAt( m_abinedges, ibin+1 ) );
            const double contrib = contribOfA( a );
            constexpr double safety = 1.05;
            nc_assert_always( contrib <= overlay*1.045 );
            if ( overlay*safety*rng.generate() <= contrib )
              break;
          }
          nc_assert_always( valueInInterval(m_alow,m_aup,a) );
          //Sample beta value. This is done linearly from pt.blow to pt.bup
          //according to s-values:
          const AlphaPtGeom pt = alphaPtGem( a );
          const PairDD svals = findSValsAtA( a );
          const double bl = pt.blow;
          const double bu = pt.bup;
          const double b1 = m_cell.b1;
          const double b2 = m_cell.b2;
          auto contrib_of_b = [&pt,&svals,b1,b2] ( double b )
          {
            //Returns S(b)*(b2_b1) = (b2-b)*s1+(b-b1)*s2;
            StableSum ss;
            ss.add(b2*svals.first);
            ss.add(-b*svals.first);
            ss.add(b*svals.second);
            ss.add(-b1*svals.second);
            return ss.sum();
          };
          const double f_bl = contrib_of_b( bl );
          const double f_bu = contrib_of_b( bu );
          //Sample from [bl,bu] according to linear fct between (bl,f_bl) to
          //(bu,f_bu). To maximise numerical stability, this is done by
          //splitting into a base rectangle or the triangle top.
          const double c_base = ncmin(f_bl,f_bu);//height of rectangle base
          const double c_tot = 0.5*(f_bl+f_bu);//av. height of base+triangle
          double x;
          if ( c_tot*rng.generate()<=c_base ) {
            x = rng.generate();
          } else {
            //just a triangle, slanted the correct way:
            x = std::sqrt(rng.generate());
            if ( f_bl > f_bu )
              x = 1.0 - x;
          }
          res.alpha = a;
          res.beta = intervalPos( bl, bu, x );
          return res;
        }
      };
    }
  }
}

NCS::RefCellSampler::RefCellSampler( const CellData& cell, double E_div_kT )
  : m_impl(static_cast<void*>(new RCSImpl(cell,E_div_kT)))
{
}

NCS::RefCellSampler::~RefCellSampler()
{
  RCSImpl* rcsimpl = static_cast<RCSImpl*>(m_impl);
  delete rcsimpl;
}

NCS::RefCellSampler::Result NCS::RefCellSampler::sampleAlphaBeta( RNG& rng )
{
  return static_cast<RCSImpl*>(m_impl)->sample(rng);
}

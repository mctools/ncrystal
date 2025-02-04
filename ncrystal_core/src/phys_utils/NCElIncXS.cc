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

#include "NCrystal/internal/phys_utils/NCElIncXS.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"

namespace NC = NCrystal;

NC::ElIncXS::ElIncXS( const VectD& elm_msd,
                      const VectD& elm_bixs,
                      const VectD& elm_scale )
{
  set( elm_msd, elm_bixs, elm_scale );
}

double NC::ElIncXS::eval_1mexpmtdivt(double t)
{
  //safe eval of (1-exp(-t))/t for t>=0.0
  nc_assert(t>=0.0);
  if ( t < 0.01 ) {
    //evaluate with Taylor expansion - for numerical stability (gives 10 sign. digits at t=0.01):
    return ( 1 + t * (  -0.5 + t * 0.16666666666666666666666666666666666666666667 * ( 1.-0.25*t ) ) );
  }
  if ( t > 24.0 ) {
    //limiting behaviour at t->inf (~10 significant digits after t>-ln(1e-10)~=23, no need for exp)
    return 1.0 / t;
  }
  //evaluate at intermediate t using actual formula (in principle no need for
  //std::expm1 from C++11 for reasons of precision since t>0.01, but it seems to
  //be significantly faster):
  t = -t;
  return std::expm1(t) / t;
}

NC::CrossSect NC::ElIncXS::evaluateMonoAtomic(NeutronEnergy ekin, double meanSqDisp, SigmaBound bound_incoh_xs)
{
  nc_assert(ekin.get()>=0.0&&meanSqDisp>=0.0&&bound_incoh_xs.get()>=0.0);
  constexpr double kkk = 16.0 * kPiSq * ekin2wlsqinv(1.0);
  return CrossSect{ bound_incoh_xs.dbl() * eval_1mexpmtdivt(kkk * meanSqDisp * ekin.dbl()) };
}

NC::ElIncXS::~ElIncXS() = default;

void NC::ElIncXS::set( const VectD& elm_msd,
                       const VectD& elm_bixs,
                       const VectD& elm_scale )
{
  //sanity check element data:
  nc_assert_always(elm_msd.size()==elm_bixs.size());
  nc_assert_always(elm_msd.size()==elm_scale.size());
  for ( auto i : ncrange( elm_msd.size() ) ) {
    nc_assert_always(elm_msd.at(i)>=0.0&&elm_msd.at(i)<1e6);
    nc_assert_always(elm_bixs.at(i)>=0.0&&elm_bixs.at(i)<1e6);
    nc_assert_always(elm_scale.at(i)>=0.0&&elm_scale.at(i)<=1e6);
  }

  m_elm_data.clear();//releases all memory since it is SmallVector
  m_elm_data.reserve_hint(elm_bixs.size());
  for ( auto i : ncrange( elm_msd.size() ) )
    m_elm_data.emplace_back( elm_msd[i], elm_bixs[i]*elm_scale[i] );
}

NC::CosineScatAngle NC::ElIncXS::sampleMuMonoAtomic( RNG& rng, NeutronEnergy ekin, double meanSqDisp )
{
  nc_assert(ekin.dbl()>=0.0&&meanSqDisp>=0.0);
  constexpr double kkk = 8.0 * kPiSq * ekin2wlsqinv(1.0);
  double twoksq = kkk * ekin.dbl();
  double a = twoksq * meanSqDisp;
  //Must sample mu in [-1,1] according to exp(a*mu). This can either happen with
  //the rejection method which is always numerically stable but slow unless a is
  //tiny, or with the transformation method which is potentially numerically
  //unstable for tiny a.

  if (a<0.01) {
    //Rejection method:

     double maxval = exp_smallarg_approx(a);
     while (true) {
       double mu = rng.generate()*2.0-1.0;
       if (rng.generate()*maxval < exp_smallarg_approx(a*mu))
         return CosineScatAngle{ mu };
     }

  } else {
    //Transformation method:

    // If f(x)=N*exp(a*x) is a normalised distribution on [-1,1], then
    // N=a/(exp(a)-exp(-a)) and the commulative probability function is F(x)=(
    // exp(a*(x+1)) -1 ) / ( exp(2*a) -1 ). With R a uniformly distributed
    // random number in (0,1], solving R=F(x) yields:
    //
    // x(R) = log( 1 + R * ( exp(2*a)-1 ) ) / a - 1
    //
    // Which can preferably be evaluated with expm1/log1p functions.
    return CosineScatAngle{ ncclamp(std::log1p( rng.generate() * std::expm1(2.0*a) ) / a - 1.0,-1.0,1.0) };
  }
}

NC::CrossSect NC::ElIncXS::evaluate(NeutronEnergy ekin) const
{
  //NB: The cross-section code here must be consistent with code in
  //evalXSContribsCommul(), evaluateMany(), and
  //evalXSContribsCommul().
  constexpr double kkk = 16.0 * kPiSq * ekin2wlsqinv(1.0);
  double e = kkk*ekin.dbl();
  double xs = 0.0;
  for ( auto& elmdata : m_elm_data )
    xs += elmdata.second * eval_1mexpmtdivt( elmdata.first * e );
  return CrossSect{ xs };
}

NC::SmallVector<double,32> NC::ElIncXS::evalXSContribsCommul( NeutronEnergy ekin ) const
{
  SmallVector<double,32> res;
  res.resize( m_elm_data.size() );
  constexpr double kkk = 16.0 * kPiSq * ekin2wlsqinv(1.0);
  double e = kkk*ekin.dbl();
  double xs = 0.0;
  unsigned i = 0;
  for ( auto& elmdata : m_elm_data )
    res[i++] = ( xs += elmdata.second * eval_1mexpmtdivt( elmdata.first * e ) );
  return res;
}

void NC::ElIncXS::evaluateMany( Span<const double> ekin, Span<double> tgt) const
{
  constexpr double kkk = 16.0 * kPiSq * ekin2wlsqinv(1.0);
  nc_assert(ekin.size()==tgt.size());
  const std::size_t n = tgt.size();

  constexpr auto nbuf = 2048;
  double buf[nbuf];
  double buf_e[nbuf];

  const double * ekin_ptr = ekin.data();
  double * tgt_ptr = tgt.data();

  std::fill( tgt_ptr, tgt_ptr + n, 0.0 );

  std::size_t nleft = n;
  while ( nleft ) {
    std::size_t nstep = std::min<std::size_t>(nleft,nbuf);
    for( auto i : ncrange(nstep) )
      buf_e[i] = kkk * ekin_ptr[i];
    for ( auto& elmdata : m_elm_data ) {
      for( auto i : ncrange(nstep) )
        buf[i] = elmdata.first * buf_e[i];
      for( auto i : ncrange(nstep) )
        buf[i] = eval_1mexpmtdivt( buf[i] );//NB: One day we want to vectorize
                                            //this!
      for( auto i : ncrange(nstep) )
        tgt_ptr[i] += elmdata.second * buf[i];
    }
    tgt_ptr += nstep;
    ekin_ptr += nstep;
    nleft -= nstep;
  }
}

NC::CosineScatAngle NC::ElIncXS::EPointAnalysis::sampleMu( const ElIncXS& ei,
                                                           RNG& rng ) const
{
  const auto& elm_data = ei.m_elm_data;
  const auto& contribs = this->data;
  nc_assert( contribs.size() == elm_data.size() );
  if ( elm_data.size() == 1 )
    return sampleMuMonoAtomic( rng, this->ekin, elm_data.front().first );
  nc_assert( contribs.size() > 1 );
  auto choiceidx = pickRandIdxByWeight( rng, contribs );
  nc_assert( choiceidx < elm_data.size() );
  return sampleMuMonoAtomic( rng, this->ekin, elm_data[choiceidx].first );
}

NC::CosineScatAngle NC::ElIncXS::sampleMu( RNG& rng, NeutronEnergy ekin )
{
  if ( m_elm_data.size() == 1 )
    return sampleMuMonoAtomic( rng, ekin, m_elm_data.front().first );
  auto contribs = evalXSContribsCommul(ekin);
  nc_assert( contribs.size() > 1 );
  nc_assert( contribs.size() == m_elm_data.size() );
  auto choiceidx = pickRandIdxByWeight( rng, contribs );
  nc_assert( choiceidx < m_elm_data.size() );
  return sampleMuMonoAtomic( rng, ekin, m_elm_data[choiceidx].first );
}

NC::ElIncXS::ElIncXS( const ElIncXS& a, double scale_a, const ElIncXS& b, double scale_b )
{
  //m_elm_data is (msd,boundincohxs*scale)
  nc_assert(scale_a>=0.0);
  using ElmDataVector = decltype(m_elm_data);
  ElmDataVector tmp;
  tmp.reserve_hint( a.m_elm_data.size() + b.m_elm_data.size() );
  auto addData = [&tmp](const ElmDataVector& data, double scale )
  {
    for ( auto& e : data ) {
      double strength = e.second * scale;
      if (strength)
        tmp.emplace_back( e.first, strength );
    }
  };
  addData(a.m_elm_data,scale_a);
  addData(b.m_elm_data,scale_b);
  std::sort(tmp.begin(),tmp.end());
  ElmDataVector result;
  result.reserve_hint(tmp.size());

  for ( const auto& e : tmp ) {
    if ( result.empty() || !floateq(result.back().first,e.first,1e-15) )
      result.emplace_back( e );
    else
      result.back().second += e.second;
  }
  result.shrink_to_fit();
  std::swap(m_elm_data,result);
}

std::size_t NC::ElIncXS::nElements() const
{
  return m_elm_data.size();
}

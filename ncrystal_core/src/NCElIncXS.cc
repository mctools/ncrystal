////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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

#include "NCrystal/internal/NCElIncXS.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCRandUtils.hh"

namespace NC = NCrystal;

NC::ElIncXS::ElIncXS( const VectD& elm_msd,
                      const VectD& elm_bixs,
                      const VectD& elm_scale )
{
  set( elm_msd, elm_bixs, elm_scale );
}

double NC::ElIncXS::evaluate(double ekin) const
{
  //NB: The cross-section code here must be consistent with code in
  //evaluateMonoAtomic() and sampleMu(..)
  constexpr double kkk = 16.0 * kPiSq * ekin2wlsqinv(1.0);
  double e = kkk*ekin;
  double xs = 0.0;
  std::vector<PairDD >::const_iterator it(m_elm_data.begin()), itE(m_elm_data.end());
  for (;it!=itE;++it)
    xs += it->second * eval_1mexpmtdivt(it->first * e);
  return xs;
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

double NC::ElIncXS::evaluateMonoAtomic(double ekin, double meanSqDisp, double bound_incoh_xs)
{
  nc_assert(ekin>=0.0&&meanSqDisp>=0.0&&bound_incoh_xs>=0.0);
  constexpr double kkk = 16.0 * kPiSq * ekin2wlsqinv(1.0);
  return bound_incoh_xs * eval_1mexpmtdivt(kkk * meanSqDisp * ekin);
}

NC::ElIncXS::~ElIncXS() = default;

void NC::ElIncXS::set( const VectD& elm_msd,
                       const VectD& elm_bixs,
                       const VectD& elm_scale )
{
  //sanity check element data:
  nc_assert_always(elm_msd.size()==elm_bixs.size());
  nc_assert_always(elm_msd.size()==elm_scale.size());
  for (std::size_t i = 0; i < elm_msd.size(); ++i) {
    nc_assert(elm_msd.at(i)>=0.0&&elm_msd.at(i)<1e6);
    nc_assert(elm_bixs.at(i)>=0.0&&elm_bixs.at(i)<1e6);
    nc_assert(elm_scale.at(i)>=0.0&&elm_scale.at(i)<=1e6);
  }

  //init:
  {
    //release old memory:
    std::vector<PairDD >().swap(m_elm_data);
  }

  m_elm_data.reserve(elm_bixs.size());
  for (std::size_t i = 0; i < elm_msd.size(); ++i) {
    m_elm_data.push_back(PairDD(elm_msd[i],elm_bixs[i]*elm_scale[i]));
  }
}

double NC::ElIncXS::sampleMuMonoAtomic( RandomBase * rng, double ekin, double meanSqDisp )
{
  nc_assert(ekin>=0.0&&meanSqDisp>=0.0);
  constexpr double kkk = 8.0 * kPiSq * ekin2wlsqinv(1.0);
  double twoksq = kkk * ekin;
  double a = twoksq * meanSqDisp;
  //Must sample mu in [-1,1] according to exp(a*mu). This can either happen with
  //the rejection method which is always numerically stable but slow unless a is
  //tiny, or with the transformation method which is potentially numerically
  //unstable for tiny a.

  if (a<0.01) {
    //Rejection method:

     double maxval = exp_smallarg_approx(a);
     while (true) {
       double mu = rng->generate()*2.0-1.0;
       if (rng->generate()*maxval < exp_smallarg_approx(a*mu))
         return mu;
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
    return ncclamp(std::log1p( rng->generate() * std::expm1(2.0*a) ) / a - 1.0,-1.0,1.0);
  }
}

double NC::ElIncXS::sampleMu( RandomBase * rng, double ekin )
{
  nc_assert(rng!=nullptr);

  const std::size_t nelem = m_elm_data.size();
  nc_assert(nelem!=0);
  if ( nelem == 1 )
    return sampleMuMonoAtomic( rng, ekin, m_elm_data.front().first );

  //Calculate per-element contribution and select accordingly.

  //First a little trick to provide us with an array for caching element-wise
  //cross-sections, without a memory allocation for normal use-cases (but
  //avoiding a hard-coded limit on number of elements).
  constexpr auto nfixed = 8;
  double cache_fixed[nfixed];
  VectD cache_dynamic;
  Span<double> elem_xs;
  if ( nelem > nfixed ) {
    cache_dynamic.resize(nelem);
    elem_xs = cache_dynamic;
  } else {
    elem_xs = Span<double>(cache_fixed).subspan(0,nelem);
  }

  //NB: The cross-section code here must be consistent with code in
  //evaluateMonoAtomic() and evaluate(..):
  constexpr double kkk = 16.0 * kPiSq * ekin2wlsqinv(1.0);
  double e = kkk*ekin;
  double xs = 0.0;
  auto itXS = elem_xs.begin();

  std::vector<PairDD >::const_iterator it(m_elm_data.begin()), itE(m_elm_data.end());
  for (;it!=itE;++it,++itXS)
    *itXS = (xs += it->second * eval_1mexpmtdivt(it->first * e));

  auto choiceidx = pickRandIdxByWeight( rng, elem_xs );//pick index according to weights (values must be commulative)
  nc_assert(choiceidx<nelem);
  return sampleMuMonoAtomic( rng, ekin, m_elm_data[choiceidx].first );
}

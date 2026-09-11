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

#include "NCrystal/internal/vdos/NCVDOSUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC=NCrystal;

NC::PairDD NC::VDOS::rangeXNexpMX(unsigned n, double eps, double accuracy ) {
  //Interval where f(x) = x^n*exp(-x) is above eps*fpeak.

  nc_assert(eps>0.0&&eps<1.0&&eps>1e-200&&n>0&&accuracy>0&&accuracy<=1e-2);

  //The function f(x) = x^n*exp(-x) peaks at x=n and falls off on both
  //sides. Returns the two solutions to f(x)= f(n)*eps, describing the central
  //range around x=n where the function is higher than eps times the peak
  //value.

  //Must solve for x:
  //  x^n*exp(-x) = eps*[n^n*exp(-n)]
  //Raise to 1/n power and get:
  //  x*exp(-x/n) = eps^(1/n)*n*exp(-1)
  //<=> (x/n)*exp(-x/n) =  (1/e)*eps^(1/n) = k
  //
  //Which can be solved numerically for x/n:

  const double fn = static_cast<double>(n);
  const double k = kInvE * std::pow(eps,1.0/fn);
  auto f = [k](double y) { return y*std::exp(-y)-k; };
  return { fn*findRoot2( f, 0.0,   1.0, accuracy ),
           fn*findRoot2( f, 1.0, 700.0, accuracy ) };
}

NC::Optional<NC::PairDD>
NC::VDOS::findExtremeSABPointWithinAlphaPlusCurve( double E_div_kT,
                                                   PairDD alphaRange,
                                                   PairDD betaRange )
{
  //Find the extreme (as in highest alpha, lowest beta) kinematically
  //accessible point in the provided rectangular region in (alpha,beta) space
  //for a neutron with energy/kT= Emax_div_kT. Returns NullOpt in case no point
  //is accessible. Note that we on purpose consider only the kinematic edge
  //given by the alpha+(beta) and beta=-E/kT curves, ignoring the alpha-(beta)
  //curve.
  nc_assert( alphaRange.second > alphaRange.first );
  nc_assert( alphaRange.first >= 0.0 );
  nc_assert( betaRange.second > betaRange.first );
  nc_assert( E_div_kT > 0.0 );

#ifndef NDEBUG
  const bool should_be_accessible
    = sabPointWithinAlphaPlusCurve(E_div_kT,alphaRange.first,betaRange.second);
#endif

  if ( betaRange.second < -E_div_kT ) {
    nc_assert(!should_be_accessible);
    return NullOpt;//no accessible points in region
  }

  auto alphaPlus = [E_div_kT](double beta)
  {
    nc_assert( beta >= -E_div_kT );
    return 2*E_div_kT + beta + 2 * std::sqrt( E_div_kT * ( E_div_kT + beta ) );
  };
  const double apb1 = alphaPlus(betaRange.second);
  if ( apb1 < alphaRange.first ) {
    nc_assert(!should_be_accessible);
    return NullOpt;//no accessible points in region
  }

  nc_assert(should_be_accessible);

  double aup = alphaRange.second;
  double blow = betaRange.first;
  const double bup = betaRange.second;

  //Clip lower beta range at -E/kT:
  blow = ncmax( blow, -E_div_kT );

  const double apb0 = alphaPlus(blow);
  if ( apb0 >= aup )
    return PairDD( aup, blow );//entire rectangle is accessible

  //Cut away excess reach of rectangular region along alpha:
  aup = ncmin(aup,apb1);

  //Cut away excess reach of rectangular region along beta:
  if ( apb0 < aup ) {
    //The next formula follows from inverting the formula
    //alphaPlus(blow) = aup (clamped for numerical imprecisions):
    blow = ncmin( bup, aup - 2.0 * std::sqrt(E_div_kT * aup));
    nc_assert(floateq(alphaPlus(blow), aup));
  }

  //Rectangular region has no excess now, result is given by its extreme
  //corner:
  return PairDD( aup, blow );
}

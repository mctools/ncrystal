////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

#include "NCRandUtils.hh"
#include "NCMath.hh"
namespace NC=NCrystal;

void NC::randIsotropicDirection( NC::RandomBase * rand, double (&dir)[3])
{
  //Very fast method (Marsaglia 1972) for generating points uniformly on the
  //unit sphere, costing approximately ~2.54 calls to rand->generate() and 1
  //call to sqrt().

  //Reference: Ann. Math. Statist. Volume 43, Number 2 (1972), 645-646.
  //           doi:10.1214/aoms/1177692644
  //Available at https://projecteuclid.org/euclid.aoms/1177692644

  double x0,x1,s;
  do {
    x0 = 2.0*rand->generate()-1.0;
    x1 = 2.0*rand->generate()-1.0;
    s = x0*x0 + x1*x1;
  } while (!s||s>=1);
  double t = 2.0*std::sqrt(1-s);
  dir[0] = x0*t;
  dir[1] = x1*t;
  dir[2] = 1.0-2.0*s;
}

void NC::randDirectionGivenScatterMu( NC::RandomBase * rand, double mu, const double(&indir)[3], double(&outdir)[3])
{
  nc_assert(ncabs(mu)<=1.);

  double m2 = indir[0]*indir[0]+indir[1]*indir[1]+indir[2]*indir[2];
  double invm = ( ncabs(m2-1.0)<1e-12 ? 1.0 : 1.0/std::sqrt(m2) );
  double ux = indir[0]*invm;
  double uy = indir[1]*invm;
  double uz = indir[2]*invm;

  //1) Create random unit-vector which is not parallel to indir:
  double tmpdir[3];

  while (true) {
    randIsotropicDirection(rand,tmpdir);
    double dotp = tmpdir[0]*ux+tmpdir[1]*uy+tmpdir[2]*uz;
    double costh2 = dotp*dotp;//tmpdir is normalised vector
    //This cut is symmetric in the parallel plane => does not ruin final
    //phi-angle-flatness:
    if (costh2<0.99)
      break;
  }
  //2) Find ortogonal vector (the randomness thus tracing a circle on the
  //unit-sphere, once normalised)
  double xx = tmpdir[1]*uz - tmpdir[2]*uy;
  double yy = tmpdir[2]*ux - tmpdir[0]*uz;
  double zz = tmpdir[0]*uy - tmpdir[1]*ux;
  double rm2 = xx*xx+yy*yy+zz*zz;

  //3) Use these two vectors to easily find the final direction (the
  //randomness above provides uniformly distributed azimuthal angle):
  double k = std::sqrt((1-mu*mu)/rm2);
  outdir[0] = ux*mu+k*xx;
  outdir[1] = uy*mu+k*yy;
  outdir[2] = uz*mu+k*zz;
}

void NC::randPointOnUnitCircle( NC::RandomBase * rand,  double & x, double& y )
{
  //Sample a random point on the unit circle. This is equivalent to sampling phi
  //randomly in [0,2pi) and letting (x,y)=(cosphi,sinphi).
  double a,b,m2;
  do {
    a = -1.0+rand->generate()*2.0;
    b = -1.0+rand->generate()*2.0;
    m2 = a*a + b*b;
  } while ( !valueInInterval(0.001,1.0,m2) );

  double m = 1.0/std::sqrt(m2);
  x = a * m;
  y = b * m;
}


double NC::randNorm( NC::RandomBase * rand )
{
  //sample a single value from a unit normal distribution via the ratio method.
  //
  //This method is ~30% faster than the randNorm(g1,g2) (using a decent Mersenne
  //Twister RNG), but only provides a single number. TODO for NC2: re-benchmark using our new rng.
  //
  //The loop runs on average ~1.37 times, and log(..) is invoked on average
  //0.167 times per call (hence the speedup).

  double g, g2(0),u,v,invu;
  do {
    u = rand->generate();
    nc_assert(u);
    invu = 1.0/u;
    v = rand->generate();
    g = 1.71552776992141354 * (v-0.5)*invu;
    g2 = g*g;
    if ( g2 <= 5.0 - 5.13610166675096558 * u )
      break;
    if ( g2 >= 1.03696104258356603 * invu )
      continue;
  } while ( g2 >= - 4.0 * std::log(u) );
  return g;
}

void NC::randNorm( NC::RandomBase * rand, double&g1, double&g2)
{
  //sample two independent values from a unit normal distribution via the polar method.
  //
  //The loop runs on average 4/pi ~= 1.27 times.

  double t;
  do {
    g1 = 2.0 * rand->generate() - 1.0;
    g2 = 2.0 * rand->generate() - 1.0;
    t = g1 * g1 + g2 * g2;
  } while ( t >= 1.0 || !t );
  t = std::sqrt( (-2.0 * std::log( t ) ) / t );
  g1 *= t;
  g2 *= t;
}

std::size_t NC::pickRandIdxByWeight( NC::RandomBase * rand, const std::vector<double>& commulvals)
{
  nc_assert(!commulvals.empty());
  std::size_t n = commulvals.size();
  if (n==1)
    return 0;
  std::vector<double>::const_iterator itB = commulvals.begin();
  std::vector<double>::const_iterator it = std::lower_bound( itB,commulvals.end(),
                                                             commulvals.back() * rand->generate() );
  return std::min<std::size_t>((std::size_t)(it-itB),n-1);
}


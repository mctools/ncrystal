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

#include "NCrystal/NCCalcBase.hh"
#include "NCVector.hh"
#include "NCMath.hh"

NCrystal::CalcBase::CalcBase(const char * calculator_type_name)
  : m_name(calculator_type_name)
{
}

NCrystal::CalcBase::~CalcBase()
{
  for (unsigned i=0;i<m_subcalcs.size();++i)
    m_subcalcs[i]->unref();
}

void NCrystal::CalcBase::registerSubCalc(CalcBase*sc)
{
  if (sc) {
    sc->ref();
    m_subcalcs.push_back(sc);
  }
}

void NCrystal::CalcBase::setRandomGenerator(NCrystal::RandomBase* rg)
{
  if (rg)
    rg->ref();
  m_randgen = rg;
  for (unsigned i=0;i<m_subcalcs.size();++i)
    m_subcalcs[i]->setRandomGenerator(rg);
}

double NCrystal::CalcBase::initDefaultRand() const
{
  nc_assert_always(!m_randgen);
  m_randgen = defaultRandomGenerator();
  nc_assert_always(m_randgen.obj());
  return m_randgen.obj()->generate();
}

double NCrystal::CalcBase::randIsotropicScatterAngle() const
{
  //acos() is not fast, but hard to come up with something faster. We could
  //consider downgrading to single-precision results, and using acosf instead of
  //acos (~twice as fast).
  const double mu = -1.0+this->rand()*2.0;
  nc_assert(ncabs(mu)<=1.0);
  return std::acos(mu);
}
double NCrystal::CalcBase::randIsotropicScatterMu() const
{
  //... or, we use randIsotropicScatterMu instead
  const double mu = -1.0+this->rand()*2.0;
  nc_assert(ncabs(mu)<=1.0);
  return mu;
}

void NCrystal::CalcBase::randIsotropicDirection(double (&dir)[3]) const
{
  //Very fast method (Marsaglia 1972) for generating points uniformly on the
  //unit sphere, costing approximately ~2.54 calls to rand() and 1 call to
  //sqrt().
  double x0,x1,s;
  do {
    x0 = 2.0*this->rand()-1.0;
    x1 = 2.0*this->rand()-1.0;
    s = x0*x0 + x1*x1;
  } while (!s||s>=1);
  double t = 2.0*std::sqrt(1-s);
  dir[0] = x0*t;
  dir[1] = x1*t;
  dir[2] = 1.0-2.0*s;
}

void NCrystal::CalcBase::randDirectionGivenScatterAngle(double th, const double(&in)[3], double(&out)[3]) const
{
  double costh,sinth;
  sincos_0pi(th,costh,sinth);
  randDirectionGivenScatterAngle(costh,sinth,in,out);
}

void NCrystal::CalcBase::randDirectionGivenScatterAngle(double costheta,
                                                        double sintheta,
                                                        const double(&indir_raw)[3],
                                                        double(&outdir_raw)[3]) const
{
  nc_assert(sizeof(Vector)==3*sizeof(double));
  const Vector& indir = asVect(indir_raw);
  Vector& outdir = asVect(outdir_raw);

  outdir = indir.unit();

  //1) Create random unit-vector which is not parallel to indir:
  Vector randdir;

  while (true) {
    randIsotropicDirection(reinterpret_cast<double(&)[3]>(randdir));
    double dotp = randdir.dot(outdir);
    double costh2 = dotp*dotp;//randdir is normalised vector
    //This cut is symmetric in the parallel plane => does not ruin final
    //phi-angle-flatness:
    if (costh2<0.99)
      break;
  }
  //2) Find ortogonal unit-vector (the randomness thus tracing a circle on the unit-sphere)
  randdir.cross_inplace(outdir);

  //3) Use these two vectors to easily find the final direction (the
  //randomness above provides the phi-randomness):
  outdir *= costheta;
  randdir *= (sintheta/randdir.mag());
  outdir += randdir;
}

void NCrystal::CalcBase::randNorm(double&g1, double&g2) const
{
  //sample two independent values from a unit normal distribution via the polar method.
  //
  //The loop runs on average 4/pi ~= 1.27 times.

  double t;
  do {
    g1 = 2.0 * this->rand() - 1.0;
    g2 = 2.0 * this->rand() - 1.0;
    t = g1 * g1 + g2 * g2;
  } while ( t >= 1.0 || !t );
  t = std::sqrt( (-2.0 * std::log( t ) ) / t );
  g1 *= t;
  g2 *= t;
}

void NCrystal::CalcBase::randNorm(double&g) const
{
  //sample a single value from a unit normal distribution via the ratio method.
  //
  //This method is ~30% faster than the randNorm(g1,g2) above (using a decent
  //Mersenne Twister RNG), but only provides a single number.
  //
  //The loop runs on average ~1.37 times, and log(..) is invoked on average
  //0.167 times per call (hence the speedup).

  double g2(0),u,v,invu;
  do {
    u = this->rand();
    if (!u)
      continue;
    invu = 1.0/u;
    v = this->rand();
    g = 1.71552776992141354 * (v-0.5)*invu;
    g2 = g*g;
    if ( g2 <= 5.0 - 5.13610166675096558 * u )
      break;
    if ( g2 >= 1.03696104258356603 * invu )
      continue;
  } while ( g2 >= - 4.0 * std::log(u) );

}

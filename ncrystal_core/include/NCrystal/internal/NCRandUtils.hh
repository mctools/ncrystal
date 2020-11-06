#ifndef NCrystal_RandUtils_hh
#define NCrystal_RandUtils_hh

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

#include "NCrystal/NCDefs.hh"

namespace NCrystal {

  //Todo: create Vector interfaces for interfaces below (can be inlined wrappers):
  double NCRYSTAL_API randIsotropicScatterAngle( RandomBase * );//has NCRYSTAL_API attribute since it is used in customphysics example
  double randIsotropicScatterMu( RandomBase * );
  void randIsotropicDirection( RandomBase *, double (&)[3]);//result will be unit vector
  void randDirectionGivenScatterMu( RandomBase *, double mu/*=cos(scatangle)*/, const double(&in)[3], double(&out)[3]);//outdir will be unit vector
  void randPointOnUnitCircle( RandomBase *,  double & x, double& y );//Sample a random point on the unit circle
  double randNorm( RandomBase * );//sample single value from unit Gaussian
  void randNorm( RandomBase *, double&g1, double&g2);//sample two independent values from unit Gaussian.
  double randNormTail(double tail, RandomBase& rng);//sample gaussian tail (tail>=0!), like sampling randNorm until result is >=tail (but more efficient).
  std::size_t pickRandIdxByWeight( RandomBase *, const VectD& commulvals);//pick index according to weights (values must be commulative)

  double randExp( RandomBase& rng );//sample single positive value from exp(-x)
  double randExpInterval( RandomBase& rng, double a, double b, double c );//Samples value in [a,b] from exp(-c*x)

  class RandExpIntervalSampler {
  public:
    //Like randExpInterval fct, but more efficient if need to sample repeatedly with same parameters.
    RandExpIntervalSampler();//invalid instance, must complete with set(..)
    void set(double a, double b, double c);
    RandExpIntervalSampler(double a, double b, double c);
    double sample(RandomBase&rng) const;
    void invalidate();
    bool isValid() const;
  private:
    double m_a, m_c1, m_c2;
  };

  //Sample f(x) = exp(-c*x)/sqrt(x) on [a,b], a>=0 b>a, c>0:
  double randExpDivSqrt( RandomBase&, double c, double a, double b );


}


////////////////////////////
// Inline implementations //
////////////////////////////

inline double NCrystal::randIsotropicScatterAngle( NCrystal::RandomBase * rand )
{
  //acos() is not fast, but hard to come up with something faster. We could
  //consider downgrading to single-precision results, and using acosf instead of
  //acos (~twice as fast).
  const double mu = -1.0+rand->generate()*2.0;
  nc_assert(mu>=-1.0&&mu<=1.0);
  return std::acos(mu);
}

inline double NCrystal::randIsotropicScatterMu( NCrystal::RandomBase * rand )
{
  //... or, just use randIsotropicScatterMu instead
  const double mu = -1.0+rand->generate()*2.0;
  nc_assert(mu>=-1.0&&mu<=1.0);
  return mu;
}

inline double NCrystal::randExp( NCrystal::RandomBase& rng )
{
  return - std::log(rng.generate());
}

inline NCrystal::RandExpIntervalSampler::RandExpIntervalSampler() : m_a(0), m_c1(0), m_c2(0)
{
}

inline void NCrystal::RandExpIntervalSampler::set(double a, double b, double c)
{
  m_a = a;
  m_c1 = -1.0/c;
  m_c2 = std::expm1(-c*(b-a));
  nc_assert(isValid());
}

inline NCrystal::RandExpIntervalSampler::RandExpIntervalSampler(double a, double b, double c)
  : m_a(a), m_c1(-1.0/c), m_c2(std::expm1(-c*(b-a)))
{
  nc_assert(isValid());
}

inline void NCrystal::RandExpIntervalSampler::invalidate()
{
  m_a = m_c1 = m_c2 = 0.0;
}

inline bool NCrystal::RandExpIntervalSampler::isValid() const
{
  return m_c1 < 0.0;
}

inline double NCrystal::RandExpIntervalSampler::sample(RandomBase&rng) const
{
  nc_assert(isValid());
  return m_a + m_c1 * std::log( 1.0 + rng.generate() * m_c2 );
}

inline double NCrystal::randExpInterval( RandomBase& rng, double a, double b, double c )
{
  return RandExpIntervalSampler(a,b,c).sample(rng);
}


#endif

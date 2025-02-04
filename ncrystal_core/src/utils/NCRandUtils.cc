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

#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
namespace NC=NCrystal;

NC::Vector NC::randIsotropicDirection( RNG& rng )
{
  //Very fast method (Marsaglia 1972) for generating points uniformly on the
  //unit sphere, costing approximately ~2.54 calls to rand->generate() and 1
  //call to sqrt().

  //Reference: Ann. Math. Statist. Volume 43, Number 2 (1972), 645-646.
  //           doi:10.1214/aoms/1177692644
  //Available at https://projecteuclid.org/euclid.aoms/1177692644

  double x0,x1,s;
  do {
    x0 = 2.0*rng.generate()-1.0;
    x1 = 2.0*rng.generate()-1.0;
    s = x0*x0 + x1*x1;
  } while (!s||s>=1);
  double t = 2.0*std::sqrt(1-s);
  return { x0*t, x1*t, 1.0-2.0*s };
}

NC::Vector NC::randDirectionGivenScatterMu( RNG& rng, double mu, const Vector& indir )
{
  nc_assert(ncabs(mu)<=1.);

  double m2 = indir.mag2();
  double invm = ( ncabs(m2-1.0)<1e-12 ? 1.0 : 1.0/std::sqrt(m2) );
  Vector u = indir * invm;

  //1) Create random unit-vector which is not parallel to indir:
  Vector tmpdir{ no_init };

  while (true) {
    tmpdir = randIsotropicDirection(rng);
    double dotp = tmpdir.dot(u);
    double costh2 = dotp*dotp;//tmpdir is normalised vector
    //This cut is symmetric in the parallel plane => does not ruin final
    //phi-angle-flatness:
    if (costh2<0.99)
      break;
  }
  //2) Find ortogonal vector (the randomness thus tracing a circle on the
  //unit-sphere, once normalised)
  double xx = tmpdir[1]*u.z() - tmpdir[2]*u.y();
  double yy = tmpdir[2]*u.x() - tmpdir[0]*u.z();
  double zz = tmpdir[0]*u.y() - tmpdir[1]*u.x();
  double rm2 = xx*xx+yy*yy+zz*zz;

  //3) Use these two vectors to easily find the final direction (the
  //randomness above provides uniformly distributed azimuthal angle):
  double k = std::sqrt((1-mu*mu)/rm2);
  u *= mu;
  return { u.x()+k*xx, u.y()+k*yy, u.z()+k*zz };
}

NC::PairDD NC::randPointOnUnitCircle( RNG& rng )
{
  //Sample a random point on the unit circle. This is equivalent to sampling phi
  //randomly in [0,2pi) and letting (x,y)=(cosphi,sinphi).
  double a,b,m2;
  do {
    a = -1.0+rng.generate()*2.0;
    b = -1.0+rng.generate()*2.0;
    m2 = a*a + b*b;
  } while ( !valueInInterval(0.001,1.0,m2) );

  double m = 1.0/std::sqrt(m2);
  return { a * m, b * m };
}

double NC::randNorm( NC::RNG& rng )
{
  //sample a single value from a unit normal distribution via the ratio method.
  //
  //This method is ~30% faster than the randNorm(g1,g2) (using a decent Mersenne
  //Twister RNG), but only provides a single number. TODO: re-benchmark using our new rng.
  //
  //The loop runs on average ~1.37 times, and log(..) is invoked on average
  //0.167 times per call (hence the speedup).

  double g, g2(0),u,v,invu;
  do {
    u = rng.generate();
    nc_assert(u);
    invu = 1.0/u;
    v = rng.generate();
    g = 1.71552776992141354 * (v-0.5)*invu;
    g2 = g*g;
    if ( g2 <= 5.0 - 5.13610166675096558 * u )
      break;
    if ( g2 >= 1.03696104258356603 * invu )
      continue;
  } while ( g2 >= - 4.0 * std::log(u) );
  return g;
}

void NC::randNorm( NC::RNG& rng, double&g1, double&g2)
{
  //sample two independent values from a unit normal distribution via the polar method.
  //
  //The loop runs on average 4/pi ~= 1.27 times.

  double t;
  do {
    g1 = 2.0 * rng.generate() - 1.0;
    g2 = 2.0 * rng.generate() - 1.0;
    t = g1 * g1 + g2 * g2;
  } while ( t >= 1.0 || !t );
  t = std::sqrt( (-2.0 * std::log( t ) ) / t );
  g1 *= t;
  g2 *= t;
}

double NC::randNormTail(double tail, NC::RNG& rng)
{
  nc_assert(tail>=0.0);
  if (tail > 0.8) {
    //"far" out on the tail, it would be inefficient to just generate normally
    //distributed numbers until one would fall above tail. Instead use
    //Marsaglias tail for normal distribution (used when implementing the
    //exponential distribution with ziggurat sampling,
    //cf. https://en.wikipedia.org/wiki/Ziggurat_algorithm):
    //
    //NB: The threshold value 0.8 came out of benchmarks. If e.g. -log(R)
    //would be replaced by something faster, this should be retuned.
    const double minvtail = 1.0/tail;
    while (true) {
      //NB: If we ever implement a faster alternative to -log(R) exponential
      //sampling, we can use it here instead!
      double x = minvtail * randExp(rng);
      double y = randExp(rng);
      if (2*y > x*x)
        return x + tail;
    }
  }
  //not far out on the tail, most efficient is to just use normal alg and
  //reject central values.
  while (true) {
#if 1
    double g1, g2;
    randNorm(rng,g1,g2);//use version providing two numbers at once,
    //i.e. optimize "worst case" performance (tail~=1)
    //rather than "best case" performance (tail~=0).
    if (std::abs(g1)>tail)
      return std::abs(g1);
    if (std::abs(g2)>tail)
      return std::abs(g2);
#else
    double g1 = randNorm(rng);
    if (std::abs(g1)>tail)
      return std::abs(g1);
#endif
  }
}

std::size_t NC::pickRandIdxByWeight( double rand01val, Span<const double> commulvals)
{
  nc_assert(!commulvals.empty());
  const auto n = commulvals.size();
  constexpr decltype(n) binary_search_threshold = 5;//NB: Not tuned, just assumed!
  static_assert( binary_search_threshold > 1,"" );
  if ( n < binary_search_threshold ) {
    if ( n==1 ) {
      //No real choice, no need to consume RNG!
      return 0;
    }
    //Small number of values => linear search:
    const double rand_choice = commulvals.back() * rand01val;
    auto itB = commulvals.begin();
    auto itE = commulvals.end();
    for (auto it = itB; it!=itE; ++it )
      if ( *it > rand_choice )
        return std::distance(itB,it);
    return n-1;
  }
  //Binary search:
  auto itB = commulvals.begin();
  auto it = std::lower_bound( itB,commulvals.end(), commulvals.back() * rand01val );
  return std::min<std::size_t>((std::size_t)(it-itB),n-1);
}

double NC::randExpDivSqrt( RNG& rng, double c, double a, double b )
{
  //Sample f(x) = exp(-c*x)/sqrt(x) on [a,b], a>=0 b>a, c>0:

  nc_assert(a>=0.0);
  nc_assert(b>a);
  nc_assert(c>0.0);

  //Note first that the MC transformation method won't do, as it leads to
  //either std::erf, std::erfc or std::expint usage. The first two have
  //precision issues for larger numbers, and the latter is only available in
  //C++17 (and it might also have precision issues.

  //So we employ the rejection method instead. In short, depending on the
  //parameters, the overlay function will be either an exponential function or
  //1/sqrt(const+x). Note that in all functions below, we liberally add and
  //remove overall normalisation factors, as they are generally unimportant as
  //long as one does not assume functions to be normalised.

  //Define A=c*a, B=c*b, U=C*(b-a) and make the variable change
  //u=c*(x-a)=c*x-A, and sample u with:
  //
  //    f(u) = exp(-u)/sqrt(A+u) on [0,c*(b-a)] = [0,B-A] = [0,U]
  //
  //Now, make another variable change, u=y-A=c*(x-a) and sample instead f(u) \propto exp(-u)/sqrt(A+u)
  //
  const double A = c*a;
  const double large_A_threshold = 0.1;//see discussions below for how to select

  if (A>large_A_threshold) {

    //Case 1, "large A".
    //
    //Scale f(u) to become f(u) = exp(-u)/sqrt(1+u/A), use exp(-u) as overlay
    //function and accept results when R<1/sqrt(1+u/A) <=> (1+u/A)*R^2<1
    //
    //NB: Calculated acceptance rates (integral of exp(-u)*(1/sqrt(1+u/A))
    //    over u=0..inf) is sqrt(pi*A)*exp(A)*erfc(sqrt(A)):
    //
    //     A         Acceptance rate
    //     0.01      15.9%
    //     0.1       40.6% [= our threshold]
    //     0.5       65.6%
    //     1.0       75.8%
    //     5.0       91.1%
    //     20.0      97.7%
    //     100.0     99.5%
    //
    // With a threshold of 0.1, the average amount of samplings needed will be
    // from 1-2.5, depending on A. Each sampling cost 1 call to std::log, in
    // addition to a one-time evaluation of an exponential during
    // initialisation.

    const double U = c*(b-a);
    const double invA = 1.0/A;
    RandExpIntervalSampler expsampler(0,U,1.0);//cost 1 std::expm1

    while (true) {
      double ugen = expsampler.sample(rng);//cost 1 std::log + 1 RNG
      double R = rng.generate();//cost 1 RNG
      if ( (1.0+ugen*invA)*R*R<1.0 )
        return ncclamp( (ugen+A)/c, a, b );//accepted, return corresponding x.
    }

  } else {

    //Case 2, "small A".
    //
    //Sample f(u) = exp(-u)/sqrt(A+u) on [0,U] using 1/sqrt(A+u) as overlay,
    //function and accept results when R<exp(-u). This is the reverse scenario
    //of Case 1 ("large A") above, and since exp(-u) has a much stronger
    //convergence to zero than than 1/sqrt(A+u), we run the risk of having the
    //1/sqrt(A+u) overlay sampling generate a huge amount of relatively high u
    //values, which can in practice never be accepted by exp(-u) (the effect
    //gets worse as A increases). We therefore restrict U to having at most a
    //value of Ulim, where Ulim is chosen as low as possible but high enough
    //that exp(-Ulim) is neglible at the desired precision. The result is that
    //we cut off the tails of the distribution to gain a dramatic increase in
    //sampling efficiency.
    //
    //It can be shown that the acceptance rate in the limit A->0 is
    //sqrt(pi/4)*erf(sqrt(U'))/sqrt(U') with U'=min(U,Ulim). This tends to 1
    //as U'->0 and to 0 as U'->infinity, underlining the need to impose Ulim.
    //
    //Calculated acceptance rates (AR):
    //
    // U' value  :    9.2    13.8    16.1     20.7    34.5
    // exp(-U')  :   1e-4    1e-6    1e-7     1e-9    1e-15
    // AR@A->0   :   29.2%   23.9%   22.1%    19.5%   15.1%
    // AR@A=0.01 :   27.1%   22.0%   20.3%    17.8%   13.8%
    // AR@A=0.1  :   23.5%   18.8%   17.3%    15.1%   11.5% [our threshold]
    // AR@A=0.5  :   19.3%   15.1%   13.8%    11.9%   8.9%
    //
    // Choosing Ulim=16.1 and a threshold of 0.1, we thus get rather uniform
    // acceptance rates of 17.3-22.1%, i.e. we have to try around ~5 times.
    // This would in principle cost us ~5 std::exp evaluations, but we avoid
    // some of them by hardcoding a few checks against known values of exp(u)
    // (see below):

    const double Ulim = 16.1180956509583;//-log(1e-7), i.e. exp(-Ulim)=1e-7.

    const double U = ncmin(c*(b-a),Ulim);
    const double B = U+A;
    if (!(B>A))
      return a;//numerical issues (A and B are both tiny), return most likely value.
    nc_assert( U>0 && B>A && A>=0 );
    const double sqrtA = std::sqrt(A);
    const double sqrtB = std::sqrt(B);
    const double sqrtB_minus_sqrtA = sqrtB - sqrtA;
    const double twosqrtA = 2*sqrtA;
    //Now, how to sample the overlay function 1/sqrt(A+u)? With g(u) \propto
    //1/sqrt(A+u) => G(u) \propto sqrt(A+u). So sampling according to
    //1/sqrt(A+u) on [0,U] can be done simply by the transformation method:
    //
    //  sqrt(A+u)-sqrt(A+0) = R*(sqrt(A+U)-sqrt(A+0))
    //  => u = (R*(sqrt(A+U)-sqrt(A))+sqrt(A))^2-A

    //  => u = (R*(sqrt(B)-sqrt(A))+sqrt(A))^2-A
    //       = (R*(sqrt(B)-sqrt(A)))^2 + 2*sqrt(A)*R*(sqrt(B)-sqrt(A))
    //       = T*(T+2*sqrtA), where T=R*(sqrt(B)-sqrt(A))
    double ugen;
    while (true) {
      const double T = rng.generate()*sqrtB_minus_sqrtA;
      ugen = T*(T+twosqrtA);
      const double Raccept = rng.generate();
      //First see if we can handle the acceptance tests cheaply (this is
      //especially important since we have a rather low acceptance rate):
      if (ugen<2.0) {
        //Out to u=2.0, exp(-u) is bracketed efficiently from above by its 6th order taylor
        //expansion, and from below by the same expansion subtracted 0.020221:
        constexpr double c1 = -1.0;
        constexpr double c2 =  1.0/2.0;
        constexpr double c3 = -1.0/6.0;
        constexpr double c4 =  1.0/24.0;
        constexpr double c5 = -1.0/120.0;
        constexpr double c6 =  1.0/720.0;
        const double taylor6 = 1.0+ugen*(c1+ugen*(c2+ugen*(c3+ugen*(c4+ugen*(c5+ugen*c6)))));
        if ( Raccept > taylor6 ) {
          nc_assert( Raccept >= std::exp(-ugen) );
          continue;//reject cheaply
        }
        if ( Raccept+0.020221 < taylor6 ) {
          nc_assert( Raccept <= std::exp(-ugen) );
          break;//accept cheaply
        }
      } else {
        nc_assert(ugen>=2.0);
        if ( Raccept > 0.135335283236614 ) {//constant is slightly larger than exp(-2.0)
          nc_assert( Raccept >= std::exp(-ugen) );
          continue;//reject cheaply
        }
        if ( ugen > 4.0 && Raccept > 0.0183156388887343 ) {//constant is slightly larger than exp(-4.0)
          nc_assert( Raccept >= std::exp(-ugen) );
          continue;//reject cheaply
        }
      }
      //Expensive check needed:
      if ( Raccept < std::exp(-ugen) )
        break;//accept
    }
    return ncclamp( (ugen+A)/c, a, b );//accepted, return corresponding x.
  }

}

////////////////

//For reference we include here the code with comments which was found on
//2018-03-28 at http://xoroshiro.di.unimi.it/xoroshiro128plus.c (tabs changed to
//2 spaces), for verification and to make it clear the the code is in the public
//domain.  Further down the file we could copy this code directly into the
//member functions of RandXRSRImpl, but we actually instead copy a slightly
//modified version found at https://github.com/skeeto/prng64-shootout (also
//public domain, CC0 1.0), simply because it is slightly more convenient.

//#if 0
///*  Written in 2016 by David Blackman and Sebastiano Vigna (vigna@acm.org)
//
//To the extent possible under law, the author has dedicated all copyright
//and related and neighboring rights to this software to the public domain
//worldwide. This software is distributed without any warranty.
//
//See <http://creativecommons.org/publicdomain/zero/1.0/>. */
//
//#include <stdint.h>
//
///* This is the successor to xorshift128+. It is the fastest full-period
//   generator passing BigCrush without systematic failures, but due to the
//   relatively short period it is acceptable only for applications with a
//   mild amount of parallelism; otherwise, use a xorshift1024* generator.
//
//   Beside passing BigCrush, this generator passes the PractRand test suite
//   up to (and included) 16TB, with the exception of binary rank tests, as
//   the lowest bit of this generator is an LFSR of degree 128. The next bit
//   can be described by an LFSR of degree 8256, but in the long run it will
//   fail linearity tests, too. The other bits needs a much higher degree to
//   be represented as LFSRs.
//
//   We suggest to use a sign test to extract a random Boolean value, and
//   right shifts to extract subsets of bits.
//
//   Note that the generator uses a simulated rotate operation, which most C
//   compilers will turn into a single instruction. In Java, you can use
//   Long.rotateLeft(). In languages that do not make low-level rotation
//   instructions accessible xorshift128+ could be faster.
//
//   The state must be seeded so that it is not everywhere zero. If you have
//   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
//   output to fill s. */
//
//uint64_t s[2];
//
//static inline uint64_t rotl(const uint64_t x, int k) {
//  return (x << k) | (x >> (64 - k));
//}
//
//uint64_t next(void) {
//  const uint64_t s0 = s[0];
//  uint64_t s1 = s[1];
//  const uint64_t result = s0 + s1;
//
//  s1 ^= s0;
//  s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); // a, b
//  s[1] = rotl(s1, 36); // c
//
//  return result;
//}
//
//
///* This is the jump function for the generator. It is equivalent
//   to 2^64 calls to next(); it can be used to generate 2^64
//   non-overlapping subsequences for parallel computations. */
//
//void jump(void) {
//  static const uint64_t JUMP[] = { 0xbeac0467eba5facb, 0xd86b048b86aa9922 };
//
//  uint64_t s0 = 0;
//  uint64_t s1 = 0;
//  for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
//    for(int b = 0; b < 64; b++) {
//      if (JUMP[i] & UINT64_C(1) << b) {
//        s0 ^= s[0];
//        s1 ^= s[1];
//      }
//      next();
//    }
//
//  s[0] = s0;
//  s[1] = s1;
//}
//#endif

//For reference we include here the code which was found on 2018-03-28 at
//http://xoroshiro.di.unimi.it/xorshift/splitmix64.c (tabs changed to 2 spaces), for
//verification and to make it clear the the code is in the public domain.
//Further down the file we use this code in a member function of RandXRSRImpl.

//#if 0
///*  Written in 2015 by Sebastiano Vigna (vigna@acm.org)
//
//To the extent possible under law, the author has dedicated all copyright
//and related and neighboring rights to this software to the public domain
//worldwide. This software is distributed without any warranty.
//
//See <http://creativecommons.org/publicdomain/zero/1.0/>. */
//
//#include <stdint.h>
//
///* This is a fixed-increment version of Java 8's SplittableRandom generator
//   See http://dx.doi.org/10.1145/2714064.2660195 and
//   http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html
//
//   It is a very fast generator passing BigCrush, and it can be useful if
//   for some reason you absolutely want 64 bits of state; otherwise, we
//   rather suggest to use a xoroshiro128+ (for moderately parallel
//   computations) or xorshift1024* (for massively parallel computations)
//   generator. */
//
//uint64_t x; /* The state can be seeded with any value. */
//
//uint64_t next() {
//  uint64_t z = (x += 0x9e3779b97f4a7c15);
//  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
//  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
//  return z ^ (z >> 31);
//}
//#endif

NC::RandXRSRImpl::RandXRSRImpl(uint64_t theseed)
{
  seed(theseed);
}

void NC::RandXRSRImpl::seed(uint64_t theseed)
{
  //Seed the state, using splitmix64 as recommended (note that the call to
  //splitmix64 actually changes the "theseed" variable, so m_s[0] and m_s[1]
  //will not be identical):
  m_s[0] = splitmix64(theseed);
  m_s[1] = splitmix64(theseed);

  //Mix up the state a little bit more, probably not really needed (NB: until
  //NCrystal v2.5 we loop 1000 times here, now only 10):
  for (unsigned i = 0; i<10; i++)
    genUInt64WithBadLowerBits();
}

void NC::RandXRSRImpl::jump()
{
  static const uint64_t JUMP[] = { 0xbeac0467eba5facb, 0xd86b048b86aa9922 };
  constexpr const int njump = static_cast<int>(sizeof(JUMP)/ sizeof(*JUMP));

  uint64_t s0 = 0;
  uint64_t s1 = 0;
  for(int i = 0; i < njump; i++)
    for(int b = 0; b < 64; b++) {
      if (JUMP[i] & UINT64_C(1) << b) {
        s0 ^= m_s[0];
        s1 ^= m_s[1];
      }
      genUInt64WithBadLowerBits();
    }

  m_s[0] = s0;
  m_s[1] = s1;
}

uint64_t NC::RandXRSRImpl::splitmix64(uint64_t& x)
{
  uint64_t z = (x += 0x9e3779b97f4a7c15);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

void NC::RandXRSRImpl::genmanyimpl(int n, double* tgt) ncnoexceptndebug
{
  //At most 1024 at a time. This doesn't actually vectorize much unfortunately.
  constexpr int nvect = 1024;
  uint64_t buf[nvect];
  nc_assert( n <= nvect );
  for ( int i = 0; i < n; ++i )
    buf[i] = this->genUInt64WithBadLowerBits();
  //NB: I tried to vectorize this next one as it has 3 internal statements that
  //can be split out. Only the third statement vectorized, so skipping for now.
  for ( int i = 0; i < n; ++i )
    tgt[i] = randUInt64ToFP01(buf[i]);
}

void NC::RandXRSRImpl::generateMany(std::size_t n, double* tgt) ncnoexceptndebug
{
  nc_assert( n>0 );
  constexpr int nvect = 1024;
  while ( n >= nvect ) {
    this->genmanyimpl(nvect,tgt);
    tgt += nvect;
    n -= nvect;
  }
  //Tail:
  if ( n )
    this->genmanyimpl(n,tgt);
}

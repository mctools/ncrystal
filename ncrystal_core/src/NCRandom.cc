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

#include "NCrystal/NCRandom.hh"
#include "NCrystal/internal/NCVector.hh"
#include "NCrystal/internal/NCMath.hh"
#include <cstdio>
#include <algorithm>

namespace NCrystal {
  static RCHolder<RandomBase> s_default_randgen;
}

void NCrystal::setDefaultRandomGenerator(RandomBase* rg)
{
  s_default_randgen = rg;
}

NCrystal::RandomBase * NCrystal::defaultRandomGenerator(bool trigger_default)
{
  if (!s_default_randgen.obj()) {
    if (!trigger_default)
      return 0;
    s_default_randgen = new RandXRSR;
  }
  return s_default_randgen.obj();
}

//For reference we include here the code with comments which was found on
//2018-03-28 at http://xoroshiro.di.unimi.it/xoroshiro128plus.c (tabs changed to
//2 spaces), for verification and to make it clear the the code is in the public
//domain.  Further down the file we could copy this code directly into the
//member functions of RandXRSR, but we actually instead copy a slightly
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
//Further down the file we use this code in a member function of RandXRSR.

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

NCrystal::RandXRSR::RandXRSR(uint64_t theseed)
{
  seed(theseed);

#define NCrystal_Random_Uint64_to_dbl (5.4210108624275215e-20)
  //NCrystal_Random_Uint64_to_dbl constant is slightly smaller than 1.0/2**64
  //(just enough that when multiplied with 2**64 it is still <1.0 at double
  //precision), which is what we need to generate numbers in [0,1). Use a few
  //asserts to sanity check this:
  uint64_t maxuintgen = std::numeric_limits<uint64_t>::max();
  double maxdblgen = maxuintgen * NCrystal_Random_Uint64_to_dbl ;
  nc_assert_always(maxdblgen < 1.0);
  nc_assert_always(maxdblgen > 0.0);
  nc_assert_always( 1.0 - maxdblgen < 2e-16 );
}

void NCrystal::RandXRSR::seed(uint64_t theseed)
{
  //Seed the state, using splitmix64 as recommended:
  m_s[0] = splitmix64(theseed);
  m_s[1] = splitmix64(theseed);

  //Mix up the state a little bit more, probably not really needed:
  for (unsigned i = 0; i<1000; i++)
    genUInt64();
}


uint64_t NCrystal::RandXRSR::genUInt64()
{
  uint64_t s0 = m_s[0];
  uint64_t s1 = m_s[1];
  uint64_t result = s0 + s1;
  s1 ^= s0;
  m_s[0] = ((s0 << 55) | (s0 >> 9)) ^ s1 ^ (s1 << 14);
  m_s[1] = (s1 << 36) | (s1 >> 28);
  return result;
}

NCrystal::RandXRSR::~RandXRSR()
{
}

double NCrystal::RandXRSR::generate()
{
  //Convert to double prec. floating point uniformly distributed in [0,1).
  return genUInt64() * NCrystal_Random_Uint64_to_dbl;
}

uint64_t NCrystal::RandXRSR::splitmix64(uint64_t& x)
{
  uint64_t z = (x += 0x9e3779b97f4a7c15);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
  z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
  return z ^ (z >> 31);
}

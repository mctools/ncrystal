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

#include "NCrystal/NCRandom.hh"
#include <cstdio>
#include <limits>

namespace NCrystal {
  static RCHolder<RandomBase> s_default_randgen;
}

NCrystal::RandomBase::RandomBase()
{
}

NCrystal::RandomBase::~RandomBase()
{
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
    printf("NCrystal WARNING: No default random generator supplied so will"
           " use the scientifically unsound NCrystal::RandomSimple.\n");
    s_default_randgen = new RandomSimple;
  }
  return s_default_randgen.obj();
}

//RandomSimple implements a very simple multiply-with-carry rand gen
//(http://en.wikipedia.org/wiki/Random_number_generation)

NCrystal::RandomSimple::RandomSimple()
  : m_w(117),/* must not be zero, nor 0x464fffff */
    m_z(11713)/* must not be zero, nor 0x9068ffff */
{
}

NCrystal::RandomSimple::~RandomSimple()
{
}

double NCrystal::RandomSimple::generate()
{
  m_w = 18000 * (m_w & 65535) + (m_w >> 16);
  m_z = 36969 * (m_z & 65535) + (m_z >> 16);
  double r = double((m_z << 16) + m_w)/double((std::numeric_limits<uint32_t>::max)());  /* 32-bit result */
  if(r==1.0)  {
    r = generate();
  }
  return r;
}

#ifndef NCrystal_Random_hh
#define NCrystal_Random_hh

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

#include "NCrystal/NCDefs.hh"

namespace NCrystal {

  //Set the default random generator which all CalcBase classes will use
  //for random number generation (unless overridden explicitly with
  //setRandomGenerator on the instance):

  void setDefaultRandomGenerator(RandomBase*);

  //Returns the global default random generator. If setDefaultRandomGenerator
  //was never called, this will trigger the creation of a RandomKISS generator
  //(see below) as the default unless trigger_default=false:

  RandomBase * defaultRandomGenerator(bool trigger_default = true);

  //Generator implementing the xoroshiro128+ (XOR/rotate/shift/rotate) due to
  //David Blackman and Sebastiano Vigna (released into public domain / CC0
  //1.0). It has a period of 2^128-1, is very fast and passes most statistical
  //tests. The one exception is that the two lowest order bits of the
  //generated integers are not of high quality. Thus, one should not simulate
  //a coin toss with "genUInt64() % 2. We thus keep the genUInt64() method
  //private and only provide callers with double precision floating points
  //uniformly distributed in [0,1):

  class RandXRSR : public RandomBase {
  public:
    RandXRSR(uint64_t seed = 0);//NB: seed = 0 is not a special seed value.
    virtual double generate();
  protected:
    virtual ~RandXRSR();
  private:
    uint64_t genUInt64();
    static uint64_t splitmix64(uint64_t& state);
    uint64_t m_s[2];
  };

}

#endif

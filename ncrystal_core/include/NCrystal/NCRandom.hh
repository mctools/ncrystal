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

#include "NCrystal/NCRCBase.hh"
#include <stdint.h>//cstdint hdr only in C++11

namespace NCrystal {

  class RandomBase : public RCBase {
  public:
    virtual double generate() = 0;//generate numbers uniformly in [0,1[
  protected:
    RandomBase();
    virtual ~RandomBase();
  };

  //Set the default random generator which all CalcBase classes will use
  //for random number generation (unless overridden explicitly with
  //setRandomGenerator on the instance):

  void setDefaultRandomGenerator(RandomBase*);

  //Returns the global default random generator. If setDefaultRandomGenerator
  //was never called, this will trigger the creation of a RandomSimple generator
  //(see below) as the default unless trigger_default=false, and will result in
  //a warning being printed to stdout.

  RandomBase * defaultRandomGenerator(bool trigger_default = true);

  //Very simple implementation which should not be used for important scientific
  //work, but which produces a highly reproducible sequence of pseudo-random
  //numbers which might be useful for validation work:
  class RandomSimple : public RandomBase {
  public:
    RandomSimple();
    virtual double generate();
  protected:
    virtual ~RandomSimple();
  private:
    uint32_t m_w;
    uint32_t m_z;
  };
}

#endif

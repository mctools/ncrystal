#ifndef NCrystal_CalcBase_hh
#define NCrystal_CalcBase_hh

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

/////////////////////////////////////////////////////////////////////
// Base class for all calculators, providing common infrastructure //
// including random sampling of various of distributions.          //
/////////////////////////////////////////////////////////////////////

#include "NCrystal/NCDefs.hh"
#include "NCrystal/NCRandom.hh"

namespace NCrystal {

  //Base class for ref-counted objects, destructors of which should be protected
  //rather than public.

  class NCRYSTAL_API CalcBase : public RCBase {
  public:
    CalcBase(const char * calculator_type_name);
    const char * getCalcName() const;
    void setCalcName(const char * n);

    bool isSubCalc(const CalcBase*) const;

    //Can be used to control and access the random number stream (but note that
    //often the setDefaultRandomGenerator(..) in NCRandom.hh will be an easier
    //way to do this globally). Sub-calcs will have their RNG changed as well:
    void setRandomGenerator(RandomBase* rg);

    //Access current RNG (the first will init and return default RNG if none was
    //set explicitly - this modifies a mutable data member):
    RandomBase* getRNG() const;//always returns valid object
    RandomBase* getRNGNoDefault() const;//returns null ptr if RNG was not set explicitly

  protected:
    //Registering sub-calcs will result in their ref-counts being incremented
    //for the mother calc's lifetime, and also means that future
    //setRandomGenerator calls will propagate to them:
    void registerSubCalc(CalcBase*);
    virtual ~CalcBase();
  private:
    std::vector<CalcBase*> m_subcalcs;
    std::string m_name;
    mutable RCHolder<RandomBase> m_randgen;
    double initDefaultRand() const;
  };
}


////////////////////////////
// Inline implementations //
////////////////////////////

inline const char * NCrystal::CalcBase::getCalcName() const
{
  return m_name.c_str();
}

inline void NCrystal::CalcBase::setCalcName(const char * n)
{
  m_name = n;
}

inline NCrystal::RandomBase* NCrystal::CalcBase::getRNG() const
{
  if (!m_randgen.obj()) {
    m_randgen = defaultRandomGenerator();
    nc_assert(m_randgen.obj());
  }
  return m_randgen.obj();
}

inline NCrystal::RandomBase* NCrystal::CalcBase::getRNGNoDefault() const
{
  return m_randgen.obj();
}


#endif

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

bool NCrystal::CalcBase::isSubCalc(const CalcBase* cb) const
{
  for (unsigned i=0;i<m_subcalcs.size();++i)
    if (cb==m_subcalcs[i])
      return true;
  return false;
}

void NCrystal::CalcBase::setRandomGenerator(NCrystal::RandomBase* rg)
{
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

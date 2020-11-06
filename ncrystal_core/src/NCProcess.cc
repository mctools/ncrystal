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

#include "NCrystal/NCProcess.hh"
#include "NCrystal/NCDefs.hh"
#include "NCrystal/internal/NCVector.hh"
#include "NCrystal/internal/NCMath.hh"

NCrystal::Process::Process(const char * calculator_type_name)
  : CalcBase(calculator_type_name)
{
}

NCrystal::Process::~Process()
{
}

double NCrystal::Process::crossSectionNonOriented(double ekin ) const
{
  if (isOriented())
    NCRYSTAL_THROW(BadInput,"Process::crossSectionNonOriented called for oriented object.");
  double indir[3] = { 0., 0., 1. };
  return crossSection(ekin, indir);
}

void NCrystal::Process::validate()
{
  double test_dir[3] = { 0., 0., 1. };
  double ekin_low;
  double ekin_high;
  double xs_low(0.0),xs_high(0.0);

  domain(ekin_low, ekin_high);
  if ( ekin_low!=ekin_low || ekin_high!=ekin_high || ekin_low<0.0 || !(ekin_low<ekin_high) ) {
    //Invalid domain, unless low=high=inf, which is a special case used to
    //indicate a process with zero cross-section everywhere.
    if ( ! (ncisinf(ekin_low) && ncisinf(ekin_high) ) )
      NCRYSTAL_THROW2(LogicError,getCalcName()<<" returns invalid domain!");
  }
  bool oriented(isOriented());
  double test_ekin((1.0-1.0e-9)*ekin_low);
  if (test_ekin < ekin_low)
    xs_low = oriented ? crossSection(test_ekin, test_dir) : crossSectionNonOriented(test_ekin);
  test_ekin = (1.0+1.0e-9)*ekin_high;
  if (test_ekin > ekin_high)
    xs_high = oriented ? crossSection(test_ekin, test_dir) : crossSectionNonOriented(test_ekin);

  if ( xs_low!=xs_low || xs_high!=xs_high || xs_low!=0.0 || xs_high!=0.0 )
    NCRYSTAL_THROW2(LogicError,getCalcName()<<" returns invalid cross sections outside domain!")
}

bool NCrystal::Process::isNull() const {
  double ekin_low, ekin_high;
  domain(ekin_low, ekin_high);
  return ncisinf(ekin_low) || !(ekin_low<ekin_high);
}

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

#include "NCrystal/NCScatter.hh"
#include "NCrystal/NCException.hh"
#include "NCVector.hh"

NCrystal::Scatter::Scatter(const char * calculator_type_name)
  : Process(calculator_type_name)
{
}

NCrystal::Scatter::~Scatter()
{
}

void NCrystal::Scatter::generateScatteringNonOriented( double ekin, double& angle, double& de ) const
{
  if (isOriented())
    NCRYSTAL_THROW(BadInput,"Scatter::generateScatteringNonOriented called for oriented object.");
  double indir[3] = { 0., 0., 1. };
  double outdir[3];
  generateScattering( ekin, indir, outdir, de );
  angle = asVect(indir).angle(asVect(outdir));
}

NCrystal::NullScatter::NullScatter()
  : Scatter("NullScatter")
{
  validate();
}

NCrystal::NullScatter::~NullScatter()
{
}

double NCrystal::NullScatter::crossSection(double, const double (&)[3] ) const
{
  return 0.0;
}

double NCrystal::NullScatter::crossSectionNonOriented( double ) const
{
  return 0.0;
}

void NCrystal::NullScatter::domain(double& ekin_low, double& ekin_high) const
{
  ekin_low = ekin_high = infinity;
}

void NCrystal::NullScatter::generateScattering( double, const double (&in)[3], double (&out)[3], double& de ) const
{
  out[0]=in[0];out[1]=in[1];out[2]=in[2];de=0;
}

void NCrystal::NullScatter::generateScatteringNonOriented( double, double& angle, double& de ) const
{
  angle = de = 0.0;
}

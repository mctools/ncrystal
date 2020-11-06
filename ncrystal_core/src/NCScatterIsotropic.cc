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

#include "NCrystal/NCScatterIsotropic.hh"
#include "NCrystal/internal/NCRandUtils.hh"

NCrystal::ScatterIsotropic::ScatterIsotropic(const char * calculator_type_name)
  : Scatter(calculator_type_name)
{
}

NCrystal::ScatterIsotropic::~ScatterIsotropic()
{
}

double NCrystal::ScatterIsotropic::crossSection( double ekin, const double (&)[3] ) const
{
  return crossSectionNonOriented(ekin);
}

void NCrystal::ScatterIsotropic::generateScattering( double ekin, const double (&indir)[3],
                                                     double (&outdir)[3], double& de ) const
{
  //Find theta and energy transfer by the non-oriented scatter process:
  double theta;
  generateScatteringNonOriented( ekin, theta, de );

  //Generate random azimuthal angle and pick outdir correspondingly:
  randDirectionGivenScatterMu(getRNG(), std::cos(theta),indir,outdir);
}

double NCrystal::ScatterIsotropic::crossSectionNonOriented( double ) const
{
  NCRYSTAL_THROW(LogicError,"ScatterIsotropic::crossSectionNonOriented must be reimplemented in derived class.");
  return 0.0;
}

void NCrystal::ScatterIsotropic::generateScatteringNonOriented( double, double& a, double& de ) const
{
  NCRYSTAL_THROW(LogicError,"ScatterIsotropic::generateScatteringNonOriented must be reimplemented in derived class.");
  a = de = 0.0;
}

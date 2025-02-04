#ifndef NCrystal_DebyeMSD_hh
#define NCrystal_DebyeMSD_hh

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

#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  //Estimate (isotropic, harmonic) atomic mean-squared-displacement using the
  //Debye Model (eq. 11+12 in R.J. Glauber, Phys. Rev. Vol98 num 6, 1955). Unit
  //of returned MSD value is Aa^2. Input temperatures should be in Kelvin, and
  //input atomic mass should be in amu.

  double debyeIsotropicMSD( DebyeTemperature, Temperature, AtomMass );

  //For debugging purposes, access the two factors in the above:

  //Calculates 1/4+x^2*integral_{0}^(1/x)[x/(exp(x)-1]dx (x=T/Tdebye)
  double calcDebyeMSDShape( double x );

  //Calculates 3*hbar^2/(M*kBoltzmann*debye_temp) (in Aa^2):
  double calcDebyeMSDScale( DebyeTemperature, AtomMass );

  //invert debyeIsotropicMSD to estimate Debye temperature from MSD (NB: Rather slow)
  DebyeTemperature debyeTempFromIsotropicMSD(double msd, Temperature, AtomMass );

}

#endif

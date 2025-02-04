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

#include "NCrystal/internal/utils/NCVector.hh"

void NCrystal::Vector::setMag(double f)
{
  if (f<0)
    NCRYSTAL_THROW(BadInput,"NCVector::setMag(): Can't set negative magnitude.");
  double themag2 = mag2();
  if (!themag2)
    NCRYSTAL_THROW(BadInput,"NCVector::setMag(): Can't scale null-vector.");
  double ff = f / sqrt(themag2);
  *this *= ff;
}

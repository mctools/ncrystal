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

#include "NCrystal/internal/phys_utils/NCDebyeMSD.hh"
#include <iostream>

int main(int,char**) {

  for (double mass: {1.00794,26.981539,238.02891})
    for (double debyetemp: {100.,410.3542,2000.})
      for (double temp: {0.,10.,150.,293.15,1000.})
        std::cout<<"sqrt(MSD[m="<<mass<<"u, Tdebye="<<debyetemp<<"K, T="<<temp<<"K]) = "
                 <<std::sqrt(NCrystal::debyeIsotropicMSD(NCrystal::DebyeTemperature{debyetemp},
                                                         NCrystal::Temperature{temp},
                                                         NCrystal::AtomMass{mass}))<<"Aa"<<std::endl;

}

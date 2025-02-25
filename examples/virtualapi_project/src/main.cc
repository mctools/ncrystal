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

#include "ncrystal_load.hh"
#include <iostream>
#include <vector>

int main()
{
  NCrystalScatProc scat_al("Al_sg225.ncmat;temp=200K");
  std::vector<double> wls = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0,
                              3.5, 4.0, 4.5, 5.0, 5.5, 6.0 };
  for ( auto wl : wls  ) {
    NCrystalScatProc::NeutronState n;
    n.ekin() = 0.081804209605330899 / ( wl*wl );
    n.u()[0] = 1.0;
    n.u()[1] = 0.0;
    n.u()[2] = 0.0;
    std::cout<<" xs("<<wl<<" Aa) = "<<scat_al.crossSection( n )<<std::endl;
  }
  return 0;
}

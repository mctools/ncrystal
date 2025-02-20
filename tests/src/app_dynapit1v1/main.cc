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

#include "NCrystal/internal/dynapi/NCDynAPIFactory.hh"
#include "NCrystal/internal/dynapi/NCDynAPI_Type1_v1.hh"
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iomanip>

int main() {
  auto dynapi = NCrystal::createDynAPI<NCrystalDynamicAPI::DynApi_Type1_v1>();
  nc_assert_always( dynapi != nullptr );

  auto scat_al = dynapi->createScatter(  "stdlib::Al_sg225.ncmat" );
  auto scat_scge = dynapi->createScatter( "stdlib::Ge_sg227.ncmat"
                                          ";dcutoff=0.5;mos=40.0arcsec"
                                          ";dir1=@crys_hkl:5,1,1@lab:0,0,1"
                                          ";dir2=@crys_hkl:0,-1,1@lab:0,1,0" );

  std::vector<double> wls = { 0.5, 1.0, 1.5, 2.0, 2.5, 3.0,
                              3.5, 4.0, 4.5, 5.0, 5.5, 6.0 };
  std::vector<double> refvals_xs_al = {
    1.39245855,
    1.37301271,
    1.37152295,
    1.29456141,
    1.11097728,
    1.05893956,
    1.38597156,
    1.76864975,
    1.40305622,
    0.143286738,
    0.148731135,
    0.154985793
  };

  auto require_flteq = [](double a, double b)
  {
    constexpr double rtol = 1.0e-6;
    constexpr double atol = 1.0e-6;
    if (!( std::fabs(a-b)
           <= 0.5 * rtol * (std::fabs(a) + std::fabs(b)) + atol ) ) {
      std::cout << std::setprecision(7)
                << "ERROR: Expected value: " << a
                << " but got " << b << std::endl;
      throw std::runtime_error("float comparison failed");
    }
  };

  double dir[3] = { 1.0, 0.0, 0.0 };//dummy dir
  auto wl2ekin = [](double wl)
  {
    return 0.081804209605330899 / (wl*wl);
  };
  auto ekin2wl = [](double ekin) {
    return std::sqrt( 0.081804209605330899 / ekin );
  };
  auto it_refval = refvals_xs_al.begin();
  for ( auto wl : wls  ) {
    double xs = dynapi->crossSectionUncached( *scat_al, wl2ekin(wl),
                                              dir[0], dir[1], dir[2] );
    std::cout<<" Al: xs(" << wl << " Aa) = "
             << std::setprecision(7)<<xs << " barn/atom" << std::endl;
    require_flteq( *it_refval++, xs );
  }


  {
    double wl = 1.540;
    double xs = dynapi->crossSectionUncached( *scat_scge, wl2ekin(wl),
                                              0.0, 1.0, 1.0 );
    std::cout<<" GeSC: xs(" << wl << " Aa, dir1) = "
             << std::setprecision(7)<<xs << " barn/atom" << std::endl;
    require_flteq( 591.0263476502018, xs );
  }
  {
    double wl = 1.540;
    double xs = dynapi->crossSectionUncached( *scat_scge, wl2ekin(wl),
                                              1.0, 1.0, 0.0 );
    std::cout<<" GeSC: xs(" << wl << " Aa, dir2) = "
             << std::setprecision(7)<<xs << " barn/atom" << std::endl;
    require_flteq( 1.667600586136298, xs );
  }


  {
    unsigned long rngstate = 1789569706;
    std::function<double()> fakerng = [&rngstate](){
      rngstate = (1103515245 * rngstate + 12345) % 2147483648;
      constexpr double f = 1.0 / 2147483648;
      return rngstate * f;
    };
    fakerng();
    fakerng();
    fakerng();
    const double wl0 = 1.540;
    double ekin = wl2ekin(wl0);
    double dir_x = 0.0;
    double dir_y = 1.0;
    double dir_z = 1.0;
    auto print_state = [&ekin, &dir_x, &dir_y, &dir_z, &ekin2wl]
    {
      std::cout<<"Neutron state: (wl="
               <<std::setprecision(5)
               << ekin2wl(ekin)
               <<" u=("<<dir_x
               <<", "<<dir_y
               <<", "<<dir_z<<")"<<std::endl;
    };

    print_state();
    dynapi->sampleScatterUncached( *scat_scge, fakerng,
                                   ekin, dir_x, dir_y, dir_z );
    print_state();
    dynapi->sampleScatterUncached( *scat_scge, fakerng,
                                   ekin, dir_x, dir_y, dir_z );
    print_state();
    dynapi->sampleScatterUncached( *scat_scge, fakerng,
                                   ekin, dir_x, dir_y, dir_z );
    print_state();
    dynapi->sampleScatterUncached( *scat_scge, fakerng,
                                   ekin, dir_x, dir_y, dir_z );
    print_state();
  }

  dynapi->deallocateScatter(  scat_al );
  dynapi->deallocateScatter(  scat_scge );
  delete dynapi;
  return 0;
}

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

#include "NCrystal/NCrystal.hh"
#include "NCrystal/internal/utils/NCVector.hh"

namespace NC=NCrystal;

void testlc(const std::string& cfg)
{
  printf("----------------- Testing cfg \"%s\"\n",cfg.c_str());
  auto scat = NC::createScatter(cfg);

  const NC::NeutronDirection neutron_dir = {0., 1./std::sqrt(2.), -1./std::sqrt(2.)};//45degree from lcaxis

  NC::VectD wls;
  wls.reserve(45);
  for(unsigned i=0;i<40;i++)
  {
    wls.push_back(0.1 + i*0.2);
  }
  wls.push_back(4.74539); //gives thetabragg=alpha_ideal=45deg for hkl=002
  wls.push_back(2.37269); //gives thetabragg=alpha_ideal=45deg for hkl=004
  wls.push_back(1.58180); //gives thetabragg=alpha_ideal=45deg for hkl=006

  for(auto it = wls.begin(); it != wls.end(); ++it)
  {
    NC::NeutronEnergy ekin = NC::NeutronWavelength{*it};
    auto xs = scat.crossSection(ekin, neutron_dir );
    printf( "xs@%geV (%g Aa): %g\n", ekin.dbl(), *it, xs.dbl() );
    if(xs&& !(std::next(it,3)<wls.end()) )//todo: why this weird thinning?
    {
      for (unsigned j=0;j<5;++j) {
        auto outcome = scat.sampleScatter( ekin , neutron_dir );
        auto v_out = outcome.direction.as<NC::Vector>();
        printf("scatter: (%g,%g,%g)=%gdeg de=%g \n", v_out[0] , v_out[1], v_out[2],
               neutron_dir.as<NC::Vector>().angle(v_out)*NC::kToDeg, (outcome.ekin.dbl()-ekin.dbl()) );
      }
    }
  }

}
int main(int , char**)
{
  //when the angle between neutron_dir and lcaxis{0,0,1} is 45deg, the optimal
  //wavelengths satisfying Bragg scattering of some major reflections are
  //4.74539 , 002
  //2.37269 , 004
  //1.58180 , 006

  std::string cfg = "C_sg194_pyrolytic_graphite.ncmat"
    ";dcutoff=0.5;mos=4deg;dir1=@crys_hkl:0,0,1@lab:0,0,1"
    ";dir2=@crys_hkl:0,1,0@lab:0,1,0;lcaxis=0,0,1;bkgd=0";
  testlc(cfg+";lcmode=-1");
  testlc(cfg+";lcmode=0");
  testlc(cfg+";lcmode=20");

  return 0;
}

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

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/dump/NCDump.hh"
//#include <iostream>

namespace NC = NCrystal;

namespace {
  const char * weirdcfg_1()
  {
    //Test wild cfg-string from https://github.com/mctools/ncrystal/issues/252:
    return "phases<0.98879661*C_sg194_pyrolytic_graphite.ncmat & "
      "0.0109976812*C_sg194_pyrolytic_graphite.ncmat & "
      "5.77146283e-09*freegas::Ag/Ag_is_Ag107/1perAa3 & "
      "5.36206901e-09*freegas::Ag/Ag_is_Ag109/1perAa3 & "
      "2.23266998e-08*freegas::B/B_is_B10/1perAa3 & "
      "8.98697733e-08*freegas::B/B_is_B11/1perAa3 & "
      "2.69571328e-06*freegas::Ca/Ca_is_Ca40/1perAa3 & "
      "1.79921448e-08*freegas::Ca/Ca_is_Ca42/1perAa3 & "
      "3.75408151e-09*freegas::Ca/Ca_is_Ca43/1perAa3 & "
      "5.80070525e-08*freegas::Ca/Ca_is_Ca44/1perAa3 & "
      "1.11231704e-10*freegas::Ca/Ca_is_Ca46/1perAa3 & "
      "5.20008445e-09*freegas::Ca/Ca_is_Ca48/1perAa3 & "
      "6.87760903e-11*freegas::Cd/Cd_is_Cd106/1perAa3 & "
      "4.89683829e-11*freegas::Cd/Cd_is_Cd108/1perAa3 & "
      "6.8720829e-10*freegas::Cd/Cd_is_Cd110/1perAa3 & "
      "7.04258691e-10*freegas::Cd/Cd_is_Cd111/1perAa3 & "
      "1.32765174e-09*freegas::Cd/Cd_is_Cd112/1perAa3 & "
      "6.72356827e-10*freegas::Cd/Cd_is_Cd113/1perAa3 & "
      "1.58070225e-09*freegas::Cd/Cd_is_Cd114/1perAa3 & "
      "4.12099324e-10*freegas::Cd/Cd_is_Cd116/1perAa3 & "
      "3.85055846e-07*freegas::Cl/Cl_is_Cl35/1perAa3 & "
      "1.23071428e-07*freegas::Cl/Cl_is_Cl37/1perAa3 & "
      "1.32454329e-08*freegas::Co/Co_is_Co59/1perAa3 & "
      "1.81648362e-08*freegas::Cr/Cr_is_Cr50/1perAa3 & "
      "3.50287304e-07*freegas::Cr/Cr_is_Cr52/1perAa3 & "
      "3.97190297e-08*freegas::Cr/Cr_is_Cr53/1perAa3 & "
      "9.88704508e-09*freegas::Cr/Cr_is_Cr54/1perAa3 & "
      "2.21712775e-13*freegas::Dy/Dy_is_Dy156/1perAa3 & "
      "3.69525129e-13*freegas::Dy/Dy_is_Dy158/1perAa3 & "
      "8.64689033e-12*freegas::Dy/Dy_is_Dy160/1perAa3 & "
      "6.98767104e-11*freegas::Dy/Dy_is_Dy161/1perAa3 & "
      "9.42653459e-11*freegas::Dy/Dy_is_Dy162/1perAa3 & "
      "9.2011147e-11*freegas::Dy/Dy_is_Dy163/1perAa3 & "
      "1.04131783e-10*freegas::Dy/Dy_is_Dy164/1perAa3 & "
      "1.88912915e-10*freegas::Eu/Eu_is_Eu151/1perAa3 & "
      "2.0622811e-10*freegas::Eu/Eu_is_Eu153/1perAa3 & "
      "3.70803046e-08*freegas::Fe/Fe_is_Fe54/1perAa3 & "
      "5.82085258e-07*freegas::Fe/Fe_is_Fe56/1perAa3 & "
      "1.34434525e-08*freegas::Fe/Fe_is_Fe57/1perAa3 & "
      "1.78896812e-09*freegas::Fe/Fe_is_Fe58/1perAa3 & "
      "7.63722108e-13*freegas::Gd/Gd_is_Gd152/1perAa3 & "
      "8.32453299e-12*freegas::Gd/Gd_is_Gd154/1perAa3 & "
      "5.65149985e-11*freegas::Gd/Gd_is_Gd155/1perAa3 & "
      "7.81670505e-11*freegas::Gd/Gd_is_Gd156/1perAa3 & "
      "5.97615975e-11*freegas::Gd/Gd_is_Gd157/1perAa3 & "
      "9.48547993e-11*freegas::Gd/Gd_is_Gd158/1perAa3 & "
      "8.34744338e-11*freegas::Gd/Gd_is_Gd160/1perAa3 & "
      "6.56618881e-08*freegas::Li/Li_is_Li6/1perAa3 & "
      "7.99457722e-07*freegas::Li/Li_is_Li7/1perAa3 & "
      "9.39982498e-08*freegas::Mn/Mn_is_Mn55/1perAa3 & "
      "6.96476064e-08*freegas::Ni/Ni_is_Ni58/1perAa3 & "
      "2.68281898e-08*freegas::Ni/Ni_is_Ni60/1perAa3 & "
      "1.16624281e-09*freegas::Ni/Ni_is_Ni61/1perAa3 & "
      "3.71839195e-09*freegas::Ni/Ni_is_Ni62/1perAa3 & "
      "9.46959232e-10*freegas::Ni/Ni_is_Ni64/1perAa3 & "
      "1.95555779e-09*freegas::S/S_is_S32/1perAa3 & "
      "1.56562054e-11*freegas::S/S_is_S33/1perAa3 & "
      "8.83719629e-11*freegas::S/S_is_S34/1perAa3 & "
      "4.11995709e-13*freegas::S/S_is_S36/1perAa3 & "
      "1.02872287e-08*freegas::Ti/Ti_is_Ti46/1perAa3 & "
      "9.27721406e-09*freegas::Ti/Ti_is_Ti47/1perAa3 & "
      "9.19248013e-08*freegas::Ti/Ti_is_Ti48/1perAa3 & "
      "6.74601816e-09*freegas::Ti/Ti_is_Ti49/1perAa3 & "
      "6.45912012e-09*freegas::Ti/Ti_is_Ti50/1perAa3 & "
      "1.27595944e-10*freegas::V/V_is_V50/1perAa3 & "
      "5.09128884e-08*freegas::V/V_is_V51/1perAa3 & "
      "0.000133306274*LiquidWaterH2O_T293.6K.ncmat & "
      "1.53326967e-08*LiquidHeavyWaterD2O_T293.6K.ncmat & "
      "6.66381704e-05*freegas::O/O_is_O16/1perAa3 & "
      "2.53315307e-08*freegas::O/O_is_O17/1perAa3>";
  }
}

int main(int,char**)
{
  NC::FactImpl::createScatter(weirdcfg_1());
  auto info = NC::FactImpl::createInfo(weirdcfg_1());
  NC::FactImpl::createAbsorption(weirdcfg_1());
  NC::dump(info);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

#include "NCrystal/internal/NCSCBragg.hh"
#include "NCrystal/NCFactImpl.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCRandUtils.hh"

namespace NC = NCrystal;

int main(int , char**)
{
  NC::SCOrientation orient;
  const NC::CrystalAxis ux = {1.,0.,0.};
  const NC::CrystalAxis uy = {0.,1.,0.};
  orient.setPrimaryDirection( ux, ux.as<NC::LabAxis>() );
  orient.setSecondaryDirection( uy,uy.as<NC::LabAxis>() );
  const NC::MosaicityFWHM mosaicity{ 1 * NCrystal::kDeg };

  //focusing on the  0   2  -2 family
  NCrystal::MatCfg cfg("Ge_sg227.ncmat;dcutoff=2;dcutoffup=2.1;inelas=0;incoh_elas=0");
  auto ci = NC::FactImpl::createInfo(cfg);
  auto sc = NC::makeSO<NC::SCBragg>(ci,orient,mosaicity);
  auto pc = NC::FactImpl::createScatter(cfg);

  auto rng  = NC::getRNG();

  //  NC::NeutronDirection incident{-1,0.1,0};

  nc_assert_always(!ci->hklList().empty());
  double thr = 2 * ci->hklList().front().dspacing;
  auto wlt_vec = NC::linspace(thr*0.9,thr*0.99,10);
  wlt_vec.insert(wlt_vec.begin(),thr*0.6);
  wlt_vec.push_back(thr*0.99999);
  wlt_vec.push_back(thr*1.00001);
  unsigned nLoop = 100000;
  //unit test of xs
  NC::CachePtr cp_sc, cp_pc;

  for(auto wl: wlt_vec)
  {
    double xss=0., xsp=0.;
    NC::NeutronEnergy ekin = NC::NeutronWavelength{wl};
    for(unsigned i=0;i<nLoop;i++)
    {
      xss += sc->crossSection( cp_sc, ekin,
                               NC::randIsotropicNeutronDirection(rng) ).get();
    }
    xss/=nLoop;
    xsp=pc->crossSectionIsotropic(cp_pc, ekin ).get();

    printf("wl %g, pc %g, sgl_cry %g, error %g%% \n", wl , xsp, xss, (xsp?(xss-xsp)/xsp*100.:(xss?9999.:0.)));
    if (wl>thr&&xss>0.0) {
      printf("Forbidden xs beyond threshold seen!\n");
      return 1;
    }
  }

  return 0;
}

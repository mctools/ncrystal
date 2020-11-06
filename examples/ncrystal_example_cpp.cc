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

#include "NCrystal/NCrystal.hh"
#include <iostream>

class MyRandGen : public NCrystal::RandomBase {
public:
  MyRandGen() : RandomBase() {}
  virtual double generate() {
    //using rand from cstdlib for this small example (note,
    //this is not really recommended for scientific work):
    return std::rand() / (RAND_MAX + 1.);
  }
protected:
  virtual ~MyRandGen() {}
};

int main() {

  NCrystal::libClashDetect();//Detect broken installation


  //Setup random generator:
  NCrystal::setDefaultRandomGenerator(new MyRandGen);

  //Create and use polycrystalline aluminium:
  const NCrystal::Scatter * pc = NCrystal::createScatter( "Al_sg225.ncmat;dcutoff=0.5;temp=25C" );
  pc->ref();

  double wl = 2.5;//angstrom
  double ekin = NCrystal::wl2ekin(wl);
  double xsect = pc->crossSectionNonOriented(ekin);
  std::cout << "polycrystal Al x-sect at " << wl << " Aa is " << xsect << " barn" << std::endl;

  double angle,delta_ekin;
  for (unsigned i=0;i<20;++i) {
    pc->generateScatteringNonOriented( ekin, angle, delta_ekin );
    std::cout <<"polycrystal random angle/delta-e at "<<wl<<" Aa is "
              <<angle*57.2957795131<<" degrees and "<<delta_ekin*1e3<<" meV"<<std::endl;
  }

  //Create and use single-crystal germanium
  NCrystal::MatCfg cfg("Ge_sg227.ncmat;dcutoff=0.5");

  cfg.set_mos(40.0*NCrystal::kArcSec);
  double c1[] = {5,1,1};
  double l1[] = {0,0,1};
  double c2[] = {0,-1,1};
  double l2[] = {0,1,0};
  cfg.set_dir1(true,c1,l1);
  cfg.set_dir2(true,c2,l2);

  const NCrystal::Scatter * sc =  NCrystal::createScatter( cfg );
  sc->ref();

  wl = 1.540;//angstrom
  ekin = NCrystal::wl2ekin(wl);
  const double dir1[3] = { 0.0, 1.0, 1.0 };
  xsect = sc->crossSection( ekin, dir1 );
  std::cout << "singlecrystal Ge x-sect at "<<wl<<" Aa is "<<xsect<<" barn (orientation 1)"<<std::endl;

  const double dir2[3] = { 1.0, 1.0, 0.0 };
  xsect = sc->crossSection( ekin, dir2 );
  std::cout << "singlecrystal Ge x-sect at "<<wl<<" Aa is "<<xsect<<" barn (orientation 2)"<<std::endl;

  //Unref to release memory (just to be tidy):
  pc->unref();
  sc->unref();
  return 0;
}

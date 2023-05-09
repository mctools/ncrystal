////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2023 NCrystal developers                                   //
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

//Include all NCrystal public headers:

#include "NCrystal/NCrystal.hh"

//Other includes (<random> is for the std::mt19937_64 example below):
#include <iostream>
#include <random>

namespace NC = NCrystal;

int main() {

  ///////////////////////////////////////////////////////
  // Sanity check the NCrystal installation (optional) //
  ///////////////////////////////////////////////////////

  //This checks that the included NCrystal headers and the linked NCrystal
  //library are from the same release of NCrystal:
  NC::libClashDetect();

  ///////////////////////////////////////
  // Setup random generator (optional) //
  ///////////////////////////////////////

  //NCrystal already ships with a high quality generator, so this is done here
  //merely as an example for users who needs to use their own source of random
  //numbers. The example wraps the C++11 mt19937_64 generator and only
  //implements the actualGenerate function. Additional functions can be
  //overridden to enable more advanced capabilities concerning multi-threading,
  //RNG state manipulation, etc. See the NCRNG.hh file for details. Of course,
  //if we had NOT registered a custom RNG source here, NCrystal would be using
  //the built-in source, which IS multi-thread safe (due to usage of jump-ahead
  //features).

  class CustomRNG : public NC::RNGStream {
    std::mt19937_64 m_gen;
  protected:
    double actualGenerate() override { return NC::randUInt64ToFP01(static_cast<uint64_t>(m_gen())); }
    //For the sake of example, we wrongly claim that this generator is safe and
    //sensible to use multithreaded (see NCRNG.hh for how to correctly deal with
    //MT safety, RNG states, etc.):
    bool useInAllThreads() const override { return true; }
  };

  //The NCrystal makeSO function is similar to std::make_shared
  //and should be used instead of raw calls to new and delete:
  auto rng = NC::makeSO<CustomRNG>();

  //Register:
  NC::setDefaultRNG(rng);

  //////////////////////////////////////
  // Create and use aluminium powder: //
  //////////////////////////////////////

  auto pc = NC::createScatter( "Al_sg225.ncmat;dcutoff=0.5;temp=25C" );

  //Specify neutron energy state using either wavelength (angstrom) or energy
  //(eV), and get the cross section (the material is isotropic so we don't need
  //to specify a direction):
  auto wl = NC::NeutronWavelength{1.8};//Or in eV as e.g. NC::NeutronEnergy{0.0253}
  auto xsect = pc.crossSectionIsotropic( wl );

  //The wl and xsect object are strongly typed, and will print out their values
  //along with units (use wl.get() or xsect.get() to access raw doubles):
  std::cout << "Powder Al x-sect at " << wl << " is " << xsect << std::endl;

  //Carry out some scatterings (again, using the isotropic material interface,
  //to get just mu=cos(theta_scat) instead of actual directions):
  for (unsigned i=0;i<5;++i) {
    auto outcome = pc.sampleScatterIsotropic( wl );
    std::cout <<"Powder Al random angle/delta-e at "<<wl<<" Aa is mu="<<outcome.mu
              <<" ("<<std::acos(outcome.mu.get())*NC::kToDeg<<" degrees) and new energy state is "
              << outcome.ekin << " ("<<outcome.ekin.wavelength()<<")" <<std::endl;
  }


  //It is cheap and easy to clone Scatter objects, which is usually done in
  //order to ensure each thread of a multi-threaded programme will use the same
  //physics model, but with its own independent RNG stream and other caches:
  auto pc_clone = pc.clone();

  (void)pc_clone;//avoid compiler warnings about unused variable

  //////////////////////////////////////////////
  // Create and use single crystal germanium: //
  //////////////////////////////////////////////

  //Create and use single crystal germanium. We orient it carefully to select a
  //scattering on the hkl=511 plane with 1.540 aangstrom:
  NC::MatCfg cfg("Ge_sg227.ncmat;dcutoff=0.5");

  cfg.set_mos( NC::MosaicityFWHM{ 40.0*NC::kArcSec });
  cfg.set_dir1( NC::HKLPoint{5,1,1}, NC::LabAxis{0,0,1} );
  cfg.set_dir2( NC::HKLPoint{0,-1,1}, NC::LabAxis{0,1,0} );

  auto sc =  NC::createScatter( cfg );

  wl = NC::NeutronWavelength{1.540};
  xsect = sc.crossSection( wl, { 0.0, 1.0, 1.0 } );
  std::cout << "single crystal Ge x-sect at "<<wl<<" Aa is "<<xsect<<" barn (orientation 1)"<<std::endl;

  xsect = sc.crossSection( wl, { 1.0, 1.0, 0.0 } );
  std::cout << "singlecrystal Ge x-sect at "<<wl<<" Aa is "<<xsect<<" barn (orientation 2)"<<std::endl;

  auto outcome1 = sc.sampleScatter( wl, { 0.0, 1.0, 1.0 } );
  std::cout << "A random scattering in orientation 1 gives: "
            << outcome1.ekin.wavelength()
            <<" and direction "<<outcome1.direction<<std::endl;

  return 0;
}

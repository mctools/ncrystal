#ifndef NCrystal_MMC_SimEngine_hh
#define NCrystal_MMC_SimEngine_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// A simulation engine which employs just a few "tricks" for variance-        //
// reduction:                                                                 //
//                                                                            //
//   * After an elastic scattering in an isotropic material, the cross        //
//     section is unchanged.                                                  //
//                                                                            //
//   * We always force scatterings and always emit transmitted neutrons       //
//     unconditionally (thus allowing more vectorisation!), but we use        //
//     russian roulette before adding the scattered states to the pending     //
//     stacks for further processing (otherwise the models would blow up and  //
//     spend infinite time on unprobable paths).                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/internal/minimc/NCMMC_RunSim.hh"
#include "NCrystal/internal/minimc/NCMMC_Baskets.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    class SimEngine : NoCopyMove {
    public:
      // A simulation engine is an engine which is able to move the simulation
      // forward one step at a time, processing a pending basket of neutrons and
      // either discard neutrons, pass them on for additional simulation steps,
      // or serve them up for tallies. This is all done via the basket manager
      // and result callback function, and for efficiency everything is done in
      // terms of baskets of neutrons.
      //
      // The class implementing this will most likely accept options including
      // for geometry, material and variance reduction in its constructor.

      virtual ~SimEngine() = default;

      //Produce a clone (most likely to have a different SimEngine object for
      //each thread in a multithreaded simulation):
      virtual shared_obj<SimEngine> clone() const = 0;

      //Advance the simulation one step. This does not need to be multi-thread
      //safe.
      virtual void step( Basket, RNG&, const TallyFct& tallyfct ) = 0;
    };

    //Create a std simulation engine through this factory function:
    shared_obj<SimEngine> createStdSimEngine( GeometryPtr,
                                              MatDef,
                                              shared_obj<BasketMgr>,
                                              const EngineOpts& opts = {} );


  }
}

#endif

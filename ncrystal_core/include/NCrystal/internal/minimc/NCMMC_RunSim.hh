#ifndef NCrystal_MMC_RunSim_hh
#define NCrystal_MMC_RunSim_hh

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

#include "NCrystal/internal/minimc/NCMMC_Tally.hh"
#include "NCrystal/internal/minimc/NCMMC_Geom.hh"
#include "NCrystal/internal/minimc/NCMMC_Source.hh"
#include "NCrystal/factories/NCMatCfg.hh"

// High level interface for a diffraction-pattern MMC application.

namespace NCRYSTAL_NAMESPACE {

  namespace MiniMC {

    struct MatDef {
      ProcPtr scatter;
      ProcPtr absorption;
      NumberDensity numDens;
      MatDef( ProcPtr scatter,
              ProcPtr absorption,
              NumberDensity nd );

      //Initialise from cfg (note, call FactoryThreadPool::enable(..) first if
      //you wish to utilise multithreading to speed up this part):
      MatDef( const MatCfg& cfg );
    };

    struct StdEngineOptions {
      //TODO: The values here are mostly guesses, and assumes initial unit
      //weights of the source particles.
      double roulette_weight_threshold = 1e-2;
      double roulette_survival_probability = 0.1;
      int roulette_nscat_threshold = 2;//particles will only get roulette'd
                                       //after this many scatterings have
                                       //already taken place.
    };

    //Launch simulations:
    void runSim_StdEngine( ThreadCount,
                           GeometryPtr,
                           SourcePtr,
                           TallyPtr,
                           MatDef,
                           StdEngineOptions = {} );

  }
}

#endif

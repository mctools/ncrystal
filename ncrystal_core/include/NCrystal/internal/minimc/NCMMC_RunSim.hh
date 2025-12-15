#ifndef NCrystal_MMC_RunSim_hh
#define NCrystal_MMC_RunSim_hh

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

#include "NCrystal/internal/minimc/NCMMC_Tally.hh"
#include "NCrystal/internal/minimc/NCMMC_Geom.hh"
#include "NCrystal/internal/minimc/NCMMC_Source.hh"
#include "NCrystal/internal/minimc/NCMMC_EngineOpts.hh"
#include "NCrystal/factories/NCMatCfg.hh"

// High level interface for a diffraction-pattern MMC application.

namespace NCRYSTAL_NAMESPACE {

  namespace MiniMC {

    struct MatDef {
      OptionalProcPtr scatter;
      OptionalProcPtr absorption;
      NumberDensity numDens;
      Optional<MatCfg> matcfg;

      MatDef( OptionalProcPtr scatter,
              OptionalProcPtr absorption,
              NumberDensity nd );

      //Initialise from cfg (note, call FactoryThreadPool::enable(..) first if
      //you wish to utilise multithreading to speed up this part):
      MatDef( const MatCfg& cfg );
    };

    struct SimOutputMetadata {
      //Various output of simulation which is not recorded in a particular tally.

      //Particles provided by source:
      ParticleCountSum provided;

      //Particles missing geometry:
      ParticleCountSum miss;
    };
    void simOutMetaDataToJSON(std::ostream&,const SimOutputMetadata&);


    //Launch simulations. Outcomes will mostly be contained in the TallyPtr
    //object afterwards, but additionally the SourcePtr might contain
    //information about total weights generated, etc.
    SimOutputMetadata runSim_StdEngine( GeometryPtr,
                                        SourcePtr,
                                        TallyPtr,
                                        MatDef,
                                        const EngineOpts& = {} );

    //After the simulation has been run, the results can be encoded into JSON by
    //passing the same objects to the following function:
    void resultsToJSON( std::ostream&,
                        GeometryPtr,
                        SourcePtr,
                        TallyPtr,
                        MatDef,
                        const EngineOpts&,
                        const SimOutputMetadata& );


  }
}

#endif

#ifndef NCrystal_MMC_SimMgr_hh
#define NCrystal_MMC_SimMgr_hh

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
#include "NCrystal/internal/minimc/NCMMC_SimEngine.hh"
#include "NCrystal/internal/minimc/NCMMC_EngineOpts.hh"
#include "NCrystal/internal/minimc/NCMMC_CBMgr.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    //Simulation manager class, responsible for actually launching worker
    //threads, cloning/merging managers, etc.

    class SimMgr final : NoCopyMove {
    public:

      SimMgr( const EngineOpts&,
              shared_obj<SimEngine>,
              shared_obj<BasketMgr>,
              shared_obj<TallyMgr>,
              Optional<CB::CBMgrInput> = NullOpt );

      ~SimMgr();

      struct LaunchSimReturnVal {
        ParticleCountSum miss;
        ParticleCountSum tallied;
      };

      LaunchSimReturnVal launchSimulations( ThreadCount nthreads,
                                            std::uint64_t seed );
    private:
      class Impl;
      Impl * m_impl;
    };
  }
}
#endif

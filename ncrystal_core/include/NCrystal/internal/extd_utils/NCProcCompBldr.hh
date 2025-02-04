#ifndef NCrystal_ProcCompBldr_hh
#define NCrystal_ProcCompBldr_hh

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

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/interfaces/NCProcImpl.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace Utils {

    //Helper class for constructing ProcComposition/ComponentList's (usually in
    //Scatter/Absorption factories), with the possibility of creating
    //sub-process components in separate threads (when NCrystal is configured
    //for it), while ensuring a deterministic ordering in the final list.

    class ProcCompBldr final : NoCopyMove {
    public:

      using ProcPtr = ProcImpl::ProcPtr;
      using ComponentList = ProcImpl::ProcComposition::ComponentList;

      ProcCompBldr();
      ~ProcCompBldr();

      //Add info, directly or via a scheduled production function which might be
      //run concurrently:
      void add( ProcPtr, double scale = 1.0 );
      void addfct_cl( std::function<ComponentList()> fct );
      void addfct( std::function<ProcPtr()> fct );
      void add_cl( ComponentList );

      //Wrap up by calling one of the finalise methods (do not call anything
      //else after this):
      ComponentList finalise();
      ProcPtr finalise_scatter();
      ProcPtr finalise_absorption();

    private:
      struct Impl;
      Pimpl<Impl> m_impl;
    };
  }
}

#endif

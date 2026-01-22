#ifndef NCrystal_ExtnFactory_hh
#define NCrystal_ExtnFactory_hh

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

#include "NCrystal/internal/utils/NCExtraTypes.hh"
#include "NCrystal/interfaces/NCProcImpl.hh"

namespace NCRYSTAL_NAMESPACE {

  class Info;

  namespace Extn {

    // Create Bragg diffraction with Extinction effects in an isotropic
    // material. The ExtinctionCfg should indicate that Extinction is enabled,
    // or an exception will be thrown.
    ProcImpl::ProcPtr createIsotropicExtnProc( PowderBraggInput::Data&&,
                                               const Cfg::ExtnCfg& );
  }

}
#endif

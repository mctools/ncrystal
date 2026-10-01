#ifndef NCrystal_DynInfoEx_hh
#define NCrystal_DynInfoEx_hh

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

#include "NCrystal/interfaces/NCInfoTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////////////////////////////////////////////////////
  // Internal intermediate base class allowing the creator of a              //
  // DI_ScatKnlDirect object to attach explicit effective-temperature and    //
  // msd values (the NCMAT v8 keywords), or to forbid automatic msd          //
  // estimation from the kernel (pre-v8 NCMAT files, whose physics must not  //
  // change). It is kept out of the public DI_ScatKnlDirect class for ABI    //
  // stability. Objects of plain DI_ScatKnlDirect provenance (e.g. built by  //
  // plugin code) get modern semantics (estimation allowed). Consumers       //
  // should never probe this class directly, but simply call                 //
  // extractTeffFromDynInfo / extractMSDFromDynInfo (NCDynInfoUtils.hh),     //
  // which encapsulate the resolution.                                       //
  //////////////////////////////////////////////////////////////////////////////

  class DI_ScatKnlDirectEx : public DI_ScatKnlDirect {
  public:
    struct ExData {
      Optional<Temperature> explicit_teff;
      Optional<double> explicit_msd;//[Aa^2], valid at the kernel temperature
      bool allow_msd_estimation = true;
    };
    DI_ScatKnlDirectEx( double fraction, IndexedAtomData atom,
                        Temperature tt, ExData ed )
      : DI_ScatKnlDirect( fraction, std::move(atom), tt ),
        m_exdata( std::move(ed) )
    {
    }
    const ExData& exData() const noexcept { return m_exdata; }
  private:
    ExData m_exdata;
  };

}

#endif

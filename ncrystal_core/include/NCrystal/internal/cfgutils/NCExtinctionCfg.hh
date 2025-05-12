#ifndef NCrystal_ExtinctionCfg_hh
#define NCrystal_ExtinctionCfg_hh

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

#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////////////////////////////////////
  // Infrastructure needed for the variables in NCCfgVars.hh  //
  //////////////////////////////////////////////////////////////

  namespace Cfg {

    class ExtinctionCfg {
    public:

      //Decode cfg data string (throws BadInput in case of syntax issues):
      ExtinctionCfg( const ExtinctionCfgData& );

      //Empty, same as decoding an empty cfg data string:
      ExtinctionCfg();

      //Re-encode into data string:
      ExtinctionCfgData encode() const;

      //Check if extinction is enabled at all:
      bool enabled() const { return m_scale.has_value(); }

      //Actual extinction model parameters (only query if enabled() is true):
      const Length& lengthScale() const { return m_scale.value(); }

    private:
      Optional<Length> m_scale;
    };

  }

}

#endif

#ifndef NCrystal_CfgExtn_hh
#define NCrystal_CfgExtn_hh

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

#include "NCrystal/internal/cfgutils/NCCfgTypes.hh"
#include "NCrystal/internal/utils/NCExtraTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  /////////////////////////////////////////////////////////////////
  // Infrastructure needed for configuration of extinction mdels //
  /////////////////////////////////////////////////////////////////

  namespace Cfg {

    namespace Extn {

      using namespace ExtnTypes;

      //To maximize flexibility, extinction configuration is passed around
      //opaguely in a VarBuf encoded by a CfgKeyValMap. This can be converted to
      //and from either cfg-string value, or model-specific data
      //structures. BadInput exceptions are used everywhere in case of invalid
      //input (like bad cfg-string syntax or inconsistent parameters).

      CfgKeyValMap decode_cfgstr( VarId, StrView );
      void stream_to_cfgstr( std::ostream&, const CfgKeyValMap& );

      ExtnCfg createExtnCfgFromVarBuf( VarBuf&& );
      const VarBuf& accessInternalVarBuf( const ExtnCfg& );

      ////////////////////////////////////////////////
      // General parameters for any model:

      enum class Model { Sabine, BC };

      //Fixme: provide raw data type (from NCExtraTypes?) without Model, and domainSize should be Optional.
      struct ExtnCfg_Base {
        Model model;
        Length domainSize;//fixme: we allow zero, but an Optional<Length> would be better in this case.
        struct Grain {
          Length grainSize;
          MosaicityFWHM angularSpread;//spread of domains inside a grain
        };
        Optional<Grain> grain;
        static ExtnCfg_Base decode( const CfgKeyValMap& );
      };

      ////////////////////////////////////////////////
      // Additional parameters for specific models:

      struct ExtnCfg_Sabine {
        //fixme replace 2 2-state enums with single 3-state enum:
        // enum class Secondary { Correlated,
        //                        UncorrelatedRectangularTilt,
        //                        UncorrelatedTriangularTilt };

        enum class Tilt { Rectangular, Triangular };
        enum class Correlation { Correlated, Uncorrelated };
        Tilt tilt = Tilt::Rectangular;
        Correlation correlation = Correlation::Correlated;
        static ExtnCfg_Sabine decode( const CfgKeyValMap& );
      };

      struct ExtnCfg_BC {
        // enum class RecipeVersion { Std2025 = 0,
        //                            Lux2025 = 1,
        //                            Classic1974 = 2 };
        // enum class SecondaryModel { Gauss, Lorentz, Fresnel };
        BC_RecipeVersion recipeVersion = BC_RecipeVersion::Default;
        BC_ScndComponent secondaryModel = BC_ScndComponent::Default;//fixme parse and allow to change

        //shuqi, type1, type2, tk (FIXME BC_ScndXDef?)

        static ExtnCfg_BC decode( const CfgKeyValMap& );
      };

      ////////////////////////////////////////////////
      // For maximum flexiblity, JSON encoding provided both detailed
      // information available on the decoded ExtnCfg_xxx objects above, as well
      // as information about the corresponding cfg-string value.  Fixme: should
      // we additionally also provide e.g. the Sabine g-value etc.?
      void stream_to_json( std::ostream&, const CfgKeyValMap& );

      ////////////////////////////////////////////////
      // Stream adapters:
      std::ostream& operator<<(std::ostream&, const ExtnCfg_Base& );
      std::ostream& operator<<(std::ostream&, const ExtnCfg_Sabine& );
      std::ostream& operator<<(std::ostream&, const ExtnCfg_BC& );

    }

  }
}

#endif

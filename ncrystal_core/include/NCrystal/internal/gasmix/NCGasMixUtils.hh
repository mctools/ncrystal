#ifndef NCrystal_GasMixUtils_hh
#define NCrystal_GasMixUtils_hh

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
#include "NCrystal/internal/utils/NCAtomUtils.hh"
#include "NCrystal/interfaces/NCAtomData.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace GasMix {

    /////////////////////////////////////////////////////////////////////////////
    // Gas mixture requests (i.e. user specifications). As objects or strings. //
    /////////////////////////////////////////////////////////////////////////////

    using ComponentList = SmallVector_IC<std::pair<double,DecodedChemForm>,6>;

    struct GasMixRequest {

      //all components EXCEPT the water molecules implied by a non-zero relative
      //humidity (it is allowed to add water molecules explicitly here, but then
      //the relativeHumidity field must be zero):
      enum class FracType { MolarFraction, MassFraction };
      FracType fracType = FracType::MolarFraction;
      ComponentList components;

      //Optionally specify a relative humidity:
      double relativeHumidity = 0.0;

      //State parameters, T and either P or rho:
      Temperature temperature;
      Variant<Pressure,Density> targetPresDens;

      //AtomDB (for enriched BF3/He, etc.) in the usual format recognised by
      //NCMAT etc.:
      std::vector<VectS> atomDBLines;
    };

    //Convert to/from a request string (useful for accepting user input, and
    //creating cache keys):
    std::string requestToString( const GasMixRequest& );
    GasMixRequest requestFromString( const std::string& );

    /////////////////////////////
    // Resulting gas mixtures. //
    /////////////////////////////

    //Results (suitable for printing, writing NCMAT data, etc):
    struct GasMixResult {
      Pressure pressure;
      Density density;
      Temperature temperature;
      //Components (always by molar fraction here):
      ComponentList components;
      //Atom data associated with each symbol in the component formulas:
      SmallVector_IC<std::pair<AtomSymbol,AtomDataSP>,8> atomDB;
    };

    //Output GasMixResult (does not print atomDB field):
    std::ostream& operator<<(std::ostream&, const GasMixResult& );

    //Analyse (and possible raise BadInput exception) a request to turn it into
    //an actual mixture (note that requested relative humidity will simply
    //result in H2O components):
    GasMixResult analyseGasMixRequest( const GasMixRequest& );

    //Saturated water pressure (i.e. partial pressure of water at 100%
    //rel. humidity).  Using a formula from
    //https://doi.org/10.1175/JAMC-D-17-0334.1 (seems to be at least precise for
    //-100C to 200C, but possibly also giving somewhat meaningful results at
    //temperatures outside this range.)
    Pressure saturatedVapourPressureOfWater( Temperature );

    //Convert to atomic level molar fractions:
    using AtomicComponentList = SmallVector_IC<std::pair<double,AtomSymbol>,6>;
    AtomicComponentList flattenComponentList( const GasMixResult& );
  }
  std::ostream& operator<<( std::ostream&, const GasMix::AtomicComponentList& );
}

#endif

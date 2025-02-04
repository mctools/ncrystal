#ifndef NCrystal_AtomDB_hh
#define NCrystal_AtomDB_hh

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

#include "NCrystal/interfaces/NCAtomData.hh"

namespace NCRYSTAL_NAMESPACE {

  //Access internal database of element and isotope data. Natural elements can
  //be accessed via Z-value or name ("Al", "H", ...). Isotopes can be accessed
  //via Z+A values or name ("C12", "Al27", "H2"). A few special isotope names
  //are additonally recognised: "D" and "T", which are the same as "H2" and "H3"
  //respectively. If requested data is not available, a null ptr is returned.

  namespace AtomDB {

    //Natural element lookup requires Z value or string like "Al", "H", ...:
    OptionalAtomDataSP getNaturalElement( unsigned Z );
    OptionalAtomDataSP getNaturalElement( const std::string& );

    //Isotope lookup requires Z+A value or string like "Al26", "He3", "D", ...:
    OptionalAtomDataSP getIsotope( unsigned Z, unsigned A );
    OptionalAtomDataSP getIsotope( const std::string& );

    //Lookup either Element or Isotope (A=0 means natural element).
    OptionalAtomDataSP getIsotopeOrNatElem( unsigned Z, unsigned A );
    OptionalAtomDataSP getIsotopeOrNatElem( const std::string& );

    //Note that the getIsotope functions can also provide natural elements if
    //called with A=0 or a chemical symbol without trailing digits.

    //Get all (Z,A) entries in DB (A=0 implies natural elements):
    unsigned getAllEntriesCount();
    std::vector<std::pair<unsigned,unsigned>> getAllEntries();

  }
}

#endif

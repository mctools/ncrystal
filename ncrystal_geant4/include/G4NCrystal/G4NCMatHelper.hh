#ifndef G4NCrystal_MatHelper_hh
#define G4NCrystal_MatHelper_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2023 NCrystal developers                                   //
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

#include "NCrystal/NCMatCfg.hh"
#include "G4Material.hh"

namespace G4NCrystal {

  using NCrystal::MatCfg;

  //Create NCrystal-enabled G4Material directly from an configuration string
  //(see NCMatCfg.hh for format). Note that for oriented crystals (single
  //crystals), any orientations specified will be interpreted in the local frame
  //of the G4LogicalVolume in which the material is installed).:
  NCRYSTAL_API G4Material * createMaterial( const char * cfgstr );
  NCRYSTAL_API G4Material * createMaterial( const G4String& cfgstr );

  //Alternatively create, configure and pass in an NCrystal MatCfg object:
  NCRYSTAL_API G4Material * createMaterial( const MatCfg&  cfg );

  //Set/disable debug output (off by default unless NCRYSTAL_DEBUG_G4MATERIALS
  //was set when the library was loaded):
  NCRYSTAL_API void enableCreateMaterialVerbosity(bool = true);
}

#endif

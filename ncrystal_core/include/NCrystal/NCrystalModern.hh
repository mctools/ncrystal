#ifndef NCrystal_NCrystalModern_hh
#define NCrystal_NCrystalModern_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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

//Special migration header which will ensure that the Modern namespace is
//injected into the NCrystal namespace and the Legacy namespace is NOT injected
//into the NCrystal namespace.

#ifdef NCrystal_NCrystalLegacy_hh
#  error Do not include both NCrystalModern.hh and NCrystalLegacy.hh in the same file.
#endif

#ifdef NCrystal_hh
#  error NCrystalModern.hh must be included BEFORE NCrystal.hh
#endif

#ifndef NCRYSTAL_USE_MODERN_INTERFACES
#  define NCRYSTAL_USE_MODERN_INTERFACES
#endif

#ifdef NCRYSTAL_USE_LEGACY_INTERFACES
#  undef NCRYSTAL_USE_LEGACY_INTERFACES
#endif

#include "NCrystal/NCrystal.hh"

#endif

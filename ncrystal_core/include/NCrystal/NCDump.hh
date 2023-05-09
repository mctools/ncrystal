#ifndef NCrystal_Dump_hh
#define NCrystal_Dump_hh

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

#include "NCrystal/NCDefs.hh"

namespace NCrystal {

  class Info;

  //Dumps info to stdout.
  //
  //  DEFAULT:  Only some information. In particular long atom position and hkl
  //            lists might be shortened. Will not necessarily trigger a full hkl
  //            list construction, so should be fast.
  //  VERBOSE1: All atom positions, both first and last hkl entries (so requires
  //            full hkl list construction), including some equivalent hkls.
  //  VERBOSE2: All hkl entries, all equivalent hkls.
  //
  enum class DumpVerbosity : unsigned { DEFAULT = 0, VERBOSE1 = 1, VERBOSE2 = 2 };
  NCRYSTAL_API void dump(const Info&, DumpVerbosity = DumpVerbosity::DEFAULT );
}

#endif

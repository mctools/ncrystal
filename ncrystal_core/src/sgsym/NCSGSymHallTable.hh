#ifndef NCrystal_SGSymHallTable_hh
#define NCrystal_SGSymHallTable_hh

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

#include "NCrystal/internal/sgsym/NCSpaceGroupHallNumber.hh"

namespace NCRYSTAL_NAMESPACE {

  //Private access to the hardwired table of the 530 space group settings
  //(see NCSGSymHallTable.cc for its origins):

  struct SGSymHallTableEntry {
    std::uint8_t number;//space group number (1..230)
    const char * choice;//setting choice code ("" if only one setting)
    const char * hall;//Hall symbol
  };

  const SGSymHallTableEntry& sgsymHallTableEntry( SpaceGroupHallNumber );

}

#endif

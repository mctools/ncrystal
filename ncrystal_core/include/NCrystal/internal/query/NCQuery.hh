#ifndef NCrystal_Query_hh
#define NCrystal_Query_hh

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

#include "NCrystal/core/NCSmallVector.hh"
#include "NCrystal/internal/utils/NCStrView.hh"

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////////////////////////////////////////////////////
  //
  // Generic query to NCrystal, the result of which will be a JSON encoded
  // string (or an exception in case of issues such as a malformed query). The
  // query takes the form of a string list, neither of which can contain the
  // ASCII BEL char (x07), or start with a dash (or whitespace followed by a
  // dash).
  //
  // This query system is primarily intended as a convenient way to get
  // information from the C++ layer to the Python layer, and the actual
  // arguments and return values are to some degree considered an
  // implementational detail which is bound to change between NCrystal versions.
  //
  using Query = SmallVector<StrView,8>;
  void JSONQuery( std::ostream&, const Query& query );

}

#endif

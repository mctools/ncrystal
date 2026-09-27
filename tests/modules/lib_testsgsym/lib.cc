
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

//Gives Python tests access to the SpaceGroup class and its hardwired table of
//settings (NB: will eventually be replaced by a JSON query).

#include "NCTestUtils/NCTestModUtils.hh"
#include "NCrystal/internal/sgsym/NCSpaceGroup.hh"
#include "NCrystal/internal/sgsym/NCSymOp.hh"

NCTEST_CTYPE_DICTIONARY
{
  return
    "int nctest_sgsym_number( int );"
    "const char * nctest_sgsym_choice( int );"
    "const char * nctest_sgsym_hallsymbol( int );"
    "const char * nctest_sgsym_tostring( int );"
    "int nctest_sgsym_parse( const char * );"
    "const char * nctest_sgsym_symop( const char * );"
    ;
}

namespace {
  NC::SpaceGroup sgFromHallNumber( int hn )
  {
    nc_assert_always( hn >= 1 && hn <= 530 );
    return NC::SpaceGroup(
      NC::SpaceGroupHallNumber{ static_cast<std::uint16_t>( hn ) } );
  }
}

NCTEST_CTYPES int nctest_sgsym_number( int hn )
{
  int res = -1;
  try {
    res = static_cast<int>( sgFromHallNumber( hn ).number() );
  } NCCATCH;
  return res;
}

NCTEST_CTYPES const char * nctest_sgsym_choice( int hn )
{
  const char * res = nullptr;
  try {
    res = sgFromHallNumber( hn ).choice();//static storage
  } NCCATCH;
  return res;
}

NCTEST_CTYPES const char * nctest_sgsym_hallsymbol( int hn )
{
  const char * res = nullptr;
  try {
    res = sgFromHallNumber( hn ).hallSymbol();//static storage
  } NCCATCH;
  return res;
}

NCTEST_CTYPES const char * nctest_sgsym_tostring( int hn )
{
  static std::string s_buf;//NB: Not thread-safe, but fine for this test
  try {
    s_buf = sgFromHallNumber( hn ).toString();
  } NCCATCH;
  return s_buf.c_str();
}

NCTEST_CTYPES int nctest_sgsym_parse( const char * str )
{
  //Returns Hall number, or 0 in case of BadInput:
  int res = -1;
  try {
    nc_assert_always( str );
    try {
      res = NC::SpaceGroup( str ).hallNumber().get();
    } catch ( NC::Error::BadInput& ) {
      res = 0;
    }
  } NCCATCH;
  return res;
}

NCTEST_CTYPES const char * nctest_sgsym_symop( const char * str )
{
  //Parses symmetry operation and returns its string representation (in
  //NCrystal's format) followed by the 9 rotation matrix entries and 3
  //translations (units of 1/24), separated by spaces. Returns "" in case of
  //BadInput.
  static std::string s_buf;//NB: Not thread-safe, but fine for this test
  try {
    nc_assert_always( str );
    try {
      const NC::SymOp op( str );
      std::ostringstream ss;
      ss << op;
      for ( unsigned i = 0; i < 3; ++i )
        for ( unsigned j = 0; j < 3; ++j )
          ss << ' ' << op.rot( i, j );
      for ( unsigned i = 0; i < 3; ++i )
        ss << ' ' << op.trans( i );
      s_buf = ss.str();
    } catch ( NC::Error::BadInput& ) {
      s_buf.clear();
    }
  } NCCATCH;
  return s_buf.c_str();
}

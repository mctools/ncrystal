
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

#include "NCTestUtils/NCTestModUtils.hh"
#include "NCrystal/internal/extn_utils/NCExtnUtils.hh"

namespace NCE = NCrystal::Extn;

NCTEST_CTYPE_DICTIONARY
{
  return
    "double nctest_calcSabineA( double );"
    "double nctest_calcSabineB( double );"
    "double nctest_calcSabineEl( double, double );"
    "double nctest_calcSabineEb( double, double );"
    ;
}

NCTEST_CTYPES double nctest_calcSabineA( double y )
{
  return NCE::calcSabineA(y);
}

NCTEST_CTYPES double nctest_calcSabineB( double y )
{
  return NCE::calcSabineB(y);
}


NCTEST_CTYPES double nctest_calcSabineEl( double x, double y )
{
  return NCE::calcSabineEl(x,y);
}

NCTEST_CTYPES double nctest_calcSabineEb( double x, double y )
{
  return NCE::calcSabineEb(x,y);
}

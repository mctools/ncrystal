
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

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/text/NCTextData.hh"
#include "NCrystal/internal/ncmat/NCParseNCMAT.hh"
#include "NCTestUtils/NCTestModUtils.hh"

NCTEST_CTYPE_DICTIONARY
{
  return
    "const char * nctest_tryParseNCMAT( const char * );"
    "const char * nctest_getTextDataLines( const char * );"
    ;
}

NCTEST_CTYPES const char * nctest_tryParseNCMAT( const char * data )
{
  //Returns ExceptionName::ExceptionMsg in case of exceptions, otherwise empty
  //string.
  static char buf[1048576];//1MB, plenty for tests
  buf[0] = '\0';
  try {
    auto td = NC::TextData( NC::RawStrData ( std::string(data) ),
                            NC::TextData::DataType{"ncmat"} );
    NC::parseNCMATData( td );
  } catch ( std::exception& e ) {
    std::ostringstream ss;
    auto ncerr = dynamic_cast<const NC::Error::Exception *>(&e);
    if ( ncerr ) {
      ss << "NC"<<ncerr->getTypeName();
    } else {
      ss << ( dynamic_cast<const std::runtime_error*>(&e)
              ? "std::runtime_error"
              : "std::exception" );
    }
    ss <<"@@@"<<e.what();
    std::string res_str = ss.str();
    if ( !res_str.empty() )
      std::strncat(buf,res_str.data(),res_str.size());
  }
  return &buf[0];
}

NCTEST_CTYPES const char * nctest_getTextDataLines( const char * data )
{
  static char buf[1048576];//1MB, plenty for tests
  buf[0] = '\0';
  try {
    std::string res_str;
    std::ostringstream ss;
    auto td = NC::TextData( NC::RawStrData ( std::string(data) ),
                            NC::TextData::DataType{"txt"} );
    bool first=true;
    for ( auto& line : td ) {
      if ( !first )
        ss << "<@@@@>";
      first = false;
      ss << line;
    }
    res_str = ss.str();
    if ( !res_str.empty() )
      std::strncat(buf,res_str.data(),res_str.size());
  } NCCATCH;
  return &buf[0];
}

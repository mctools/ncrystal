#ifndef NCrystal_String_hh
#define NCrystal_String_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

#include "NCrystal/NCException.hh"
#include <string>
#include <vector>

//Utilities

namespace NCrystal {

  //All bytes must be in range 32..126 (plus optionally new-lines and tabs).
  bool isSimpleASCII(const std::string&, bool allow_tab=false, bool allow_newline = false);

  //Strip excess whitespace (" \t\r\n") from both ends of string:
  void trim( std::string& );

  //Split input string on separator (default sep=0 means splitting on general
  //whitespace - " \t\r\n"). Results are placed in output vector, which is first
  //cleared of existing contents. Empty parts are only kept when sep!=0 (similar
  //to pythons str.split()). Finally, maxsplit can be used to limit the number
  //of splittings performed:
  void split(std::vector<std::string>& output,
             const std::string& input,
             std::size_t maxsplit = 0,
             char sep = 0 );

  //Get basename and extension from filename:
  std::string basename(const std::string& filename);
  std::string getfileext(const std::string& filename);

  bool startswith(const std::string& str, const std::string& substr);

  //Check if given char or substring (needle) is present in a string (haystack):
  bool contains(const std::string& haystack, char needle );
  bool contains(const std::string& haystack, const std::string& needle);
  //Check if any of the chars in "needles" is present in the string (haystack):
  bool contains_any(const std::string& haystack, const std::string& needles);

  //Convert strings to numbers. In case of problems, a BadInput exception will
  //be thrown (provide err to modify the message in that exception):
  double str2dbl(const std::string&, const char * errmsg = 0);
  int str2int(const std::string&, const char * errmsg = 0);
}

#endif

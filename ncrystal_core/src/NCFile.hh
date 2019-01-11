#ifndef NCrystal_File_hh
#define NCrystal_File_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2019 NCrystal developers                                   //
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

#include <string>

namespace NCrystal {

  //Check if file exists and is readable:
  bool file_exists ( const std::string& filename );

  //Search for the file. If filename does not exist relative to the current
  //working directory, is not an absolute path to the file, first look
  //relatively to the directory NCRYSTAL_DATADIR and secondly relatively to a
  //directory (if any) given at compilation time by the preprocessor define
  //-DNCRYSTAL_DATADIR=/some/dir. Returns empty string if not found:
  std::string find_file( const std::string& filename );

}

#endif

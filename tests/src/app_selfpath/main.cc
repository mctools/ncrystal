////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

#include "NCrystal/internal/NCFileUtils.hh"
#include "NCrystal/internal/NCMath.hh"//isOneOf
#include <iostream>
namespace NC=NCrystal;

int main( int argc, char ** argv )
{
  auto self_path = NC::determine_exe_self_path( argc, argv );
  std::cout<<"determine_exe_self_path: "<<self_path<<std::endl;
  if (!NC::file_exists(self_path)) {
    std::cout<<"Error: did not exist!"<<std::endl;
    return 1;
  }

  if (!NC::path_is_absolute(self_path)) {
    std::cout<<"Error: was not absolute!"<<std::endl;
    return 1;
  }

  nc_assert_always( NC::isOneOf( NC::basename(self_path),
                                 "selfpath",
                                 "selfpath.exe") );
  //  nc_assert_always( NC::basename(NC::dirname(self_path)) == "bin" );
  //fixme figure out what to expect instead of "bin" ^^^^^^


  //tryRealPath
  //std::string normalise(const std::string& path);
  std::cout<<"All looks OK"<<std::endl;
  return 0;
}


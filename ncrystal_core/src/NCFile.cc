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

#include "NCFile.hh"

#include <fstream>
#include <cstdlib>

namespace NCrystal {

  std::string path_join(const std::string p1, const std::string p2)
  {
#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(_WIN32) || defined(__CYGWIN__)
    return p1+'\\'+p2;
#else
    return p1+'/'+p2;
#endif

  }
}
bool NCrystal::file_exists(const std::string& name) {
  //Only portable way in C++98 is to attempt to open the file.
  std::ifstream f(name.c_str());
  return f.good();
}

std::string NCrystal::find_file(const std::string& filename) {

  if (filename.empty())
    return std::string();

  if (file_exists(filename))
    return filename;

  if (filename.at(0)=='/')
    return std::string();//don't look further if absolute path (this test might fail on windows...).

  const char * envpath = std::getenv("NCRYSTAL_DATADIR");
  std::string path = envpath ? envpath : "";
  if (!path.empty()) {
    std::string tmp = path_join(path,filename);
    if (file_exists(tmp))
      return tmp;
  }
#ifdef NCRYSTAL_DATADIR
#  define NCRYSTAL_str(s) #s
#  define NCRYSTAL_xstr(s) NCRYSTAL_str(s)
  path = NCRYSTAL_xstr(NCRYSTAL_DATADIR);
  if (!path.empty()) {
    std::string tmp = path_join(path,filename);
    if (file_exists(tmp))
      return tmp;
  }
#endif


  //not found.
  return std::string();
}

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

#include "NCrystal/interfaces/NCVersion.hh"

int NCrystal::getVersion()
{
  return NCRYSTAL_VERSION;
}

const char* NCrystal::getBuildNameSpace()
{
#ifdef ncrystal_str
#  undef ncrystal_str
#endif
#ifdef ncrystal_xstr
#  undef ncrystal_xstr
#endif
#define ncrystal_str(s) #s
#define ncrystal_xstr(s) ncrystal_str(s)

#ifdef NCRYSTAL_NAMESPACE_PROTECTION
  return ncrystal_xstr(NCRYSTAL_NAMESPACE_PROTECTION);
#else
  return "";
#endif
}

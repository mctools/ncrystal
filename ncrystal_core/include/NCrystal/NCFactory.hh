#ifndef NCrystal_Factory_hh
#define NCrystal_Factory_hh

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

#include "NCrystal/NCMatCfg.hh"

namespace NCrystal {

  class MatCfg;
  class Info;
  class Scatter;
  class Absorption;

  //Generic interface for transforming user configuration (in the form of MatCfg
  //objects) into Info, Scatter or Absorption objects. The interface always
  //returns valid objects with a reference count of zero or one (in cache
  //caching is enabled), and might throw exceptions in case of errors:

  const Info * createInfo( const MatCfg& );
  const Scatter * createScatter( const MatCfg& );
  const Absorption * createAbsorption( const MatCfg& );

  //For convenience, and when the OO interface of the MatCfg class is not
  //needed, factories can work on strings as well:

  inline const Info * createInfo( const char * c ) { return createInfo(MatCfg(c)); }
  inline const Scatter * createScatter( const char * c ) { return createScatter(MatCfg(c)); }
  inline const Absorption * createAbsorption( const char * c ) { return createAbsorption(MatCfg(c)); }

  //To avoid expensive re-generation of Info objects, these are cached behind
  //the scenes The following function can be used to clear the cache and
  //potentially free up some memory:
  void clearInfoCaches();

  //disable and enable caching (default state upon startup is for caching to be
  //enabled, unless the environment variable NCRYSTAL_NOCACHE is set):
  void disableCaching();
  void enableCaching();

  //Note: If trying to debug factory availability and createInfo caching, it
  //might be useful to set the environment variable NCRYSTAL_DEBUGFACTORY=1 in
  //order to get verbose printouts of what goes on behind the scenes.

}


#endif

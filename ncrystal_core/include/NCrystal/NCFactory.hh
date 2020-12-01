#ifndef NCrystal_Factory_hh
#define NCrystal_Factory_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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
  //returns valid objects with a reference count of zero or one (in case
  //caching is enabled), and might throw exceptions in case of errors:

  NCRYSTAL_API const Info * createInfo( const MatCfg& );
  NCRYSTAL_API const Scatter * createScatter( const MatCfg& );
  NCRYSTAL_API const Absorption * createAbsorption( const MatCfg& );

  //For convenience, and when the OO interface of the MatCfg class is not
  //needed, factories can work on strings as well:

  NCRYSTAL_API inline const Info * createInfo( const char * c ) { return createInfo(MatCfg(c)); }
  NCRYSTAL_API inline const Scatter * createScatter( const char * c ) { return createScatter(MatCfg(c)); }
  NCRYSTAL_API inline const Absorption * createAbsorption( const char * c ) { return createAbsorption(MatCfg(c)); }

  //To avoid expensive re-generation of Info objects, these are cached behind
  //the scenes based on the *name* of the input file as well as the values of
  //the MatCfg parameters affecting Info creation. The following function can be
  //used to clear the cache and potentially free up some memory:
  NCRYSTAL_API void clearInfoCaches();

  //Disable and enable caching (default state upon startup is for caching to be
  //enabled, unless the environment variable NCRYSTAL_NOCACHE is set):
  NCRYSTAL_API void disableCaching();
  NCRYSTAL_API void enableCaching();

  //Note: If trying to debug factory availability and createInfo caching, it
  //might be useful to set the environment variable NCRYSTAL_DEBUGFACTORY=1 in
  //order to get verbose printouts of what goes on behind the scenes.

  ////////////////////////////////////////////////////////////////////////////
  //Register in-memory file data. This needs a "filename" and the content of
  //this virtual file. After registering such in-memory "files", they can be
  //used as file names in cfg strings or MatCfg objects (for input types which
  //support it, which certainly includes NCMAT file data, for which even the
  //virtual "filename" should end with ".ncmat"). Registering the same filename
  //more than once, will simply override the content:

  NCRYSTAL_API void registerInMemoryFileData( const std::string& virtual_filename,
                                              const std::string& data );
  NCRYSTAL_API void registerInMemoryFileData( const std::string& virtual_filename,
                                              std::string&& data );

  //Version which just registers the data address but does not copy the data
  //(intended for efficiently hard-coding a large database directly in C/C++
  //code):
  NCRYSTAL_API void registerInMemoryStaticFileData( const std::string& virtual_filename,
                                                    const char* static_data );

  //Call this function ensures that any embedded data (via the standard
  //EMBED_DATA option in NCrystal's CMake) is registered. Calling the usual
  //high-level API with createXXX(..) and/or NCMatCfg objects, will
  //automatically call this function.
  NCRYSTAL_API void ensureEmbeddedDataIsRegistered();

  //WARNING: Calling these functions will the first time result in a global
  //TextInputManager being registered (see NCFile.hh). If your application needs
  //its own custom TextInputManger, you should not call any of these functions.
}


#endif

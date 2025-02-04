#ifndef NCrystal_Fact_hh
#define NCrystal_Fact_hh

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

#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/interfaces/NCProc.hh"
#include "NCrystal/interfaces/NCInfo.hh"

namespace NCRYSTAL_NAMESPACE {

  ///////////////////////////////////////////////////////////////////////////////
  // Factory functions needed by most clients.                                  //
  //                                                                            //
  // These factory functions are thin wrappers around the factory functions     //
  // from the FactImpl namespace (in NCFactImpl.hh). For the case of Scatter    //
  // and Absorption classes, they add caching and (for Scatter) RNG stream      //
  // handling, on top of the process implementations from the ProcImpl          //
  // namespace (in NCProcImpl.hh). For multi-threaded applications, one can use //
  // the cloneXXX() methods of the Absorption and Scatter classes, or simply    //
  // use the ProcImpl classes directly, along with the RNG stream handling from //
  // NCRNG.hh.                                                                  //
  ////////////////////////////////////////////////////////////////////////////////

  NCRYSTAL_API shared_obj<const Info> createInfo( const MatCfg& cfg );
  NCRYSTAL_API Scatter createScatter( const MatCfg& cfg );
  NCRYSTAL_API Absorption createAbsorption( const MatCfg& cfg );

  //////////////////////////////////////////////////////////////////////////
  // For the case of Scatter instance, they can also by created with more //
  // control over the RNG streams (again, this is just thin wrappers over //
  // functionaliy available in NCRNG.hh+NCFactImpl.hh+NCProc.hh):         //
  //////////////////////////////////////////////////////////////////////////

  NCRYSTAL_API Scatter createScatter_RNGByIdx( const MatCfg& cfg, RNGStreamIndex rngidx );
  NCRYSTAL_API Scatter createScatter_RNGForCurrentThread( const MatCfg& cfg );

  //////////////////////////////////////////////////////////////////////////
  // Register in-memory data files which can later be referred to in      //
  // cfg strings. Note that the file NCDataSources.hh provides MANY more  //
  // options and control of input data.                                   //
  //////////////////////////////////////////////////////////////////////////

  NCRYSTAL_API void registerInMemoryFileData( std::string virtualFileName,
                                              std::string&& data );

  //Version which just registers the address of data and does not copy
  //it. Naturally the lifetime of the static_data should be longer than any
  //direct or indirect calls to FactImpl::createTextData (it is intended for
  //efficiently hard-coding file content in C/C++ code):
  NCRYSTAL_API void registerInMemoryStaticFileData( std::string virtualFileName,
                                                    const char* static_data );

}

#endif

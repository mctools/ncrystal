#ifndef NCrystal_MsgCtrl_hh
#define NCrystal_MsgCtrl_hh

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

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////////////////////////////////////////////////////
  //
  // By default, all NCrystal output is written to stdout via usage of the
  // std::cout stream. However, this can be modified by registering a call to
  // the following function:
  //
  // Passing an empty function is the same as restoring the default behaviour,
  // of passing all message strings to std::cout. Specifically, like this:
  //
  // MsgType::Info     :   std::cout << "NCrystal: " << msg << std::endl;
  // MsgType::Warning  :   std::cout << "NCrystal WARNING: " << msg << std::endl;
  // MsgType::RawOutput:   std::cout << "NCrystal: " << msg << std::flush;

  using MsgHandlerFct_t = std::function<void(const char*,MsgType)>;
  NCRYSTAL_API void setMessageHandler( MsgHandlerFct_t = nullptr );

}

 #endif

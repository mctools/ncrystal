#ifndef NCrystal_Msg_hh
#define NCrystal_Msg_hh

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

#include "NCrystal/NCDefs.hh"
#include "NCrystal/NCTypes.hh"
#include <sstream>

namespace NCRYSTAL_NAMESPACE {
  namespace Msg {

    //All NCrystal code should always use the following functions (or the
    //macro's below) in favour of direct std::cout usage. This not only ensures
    //proper capture of all output, but also makes output thread-safe due to
    //projection by a mutex).

    void outputMsg( const char * msg, MsgType = MsgType::Info );

    inline void outputMsg( const std::string& msg, MsgType mt = MsgType::Info )
    {
      outputMsg( msg.c_str(), mt );
    }

    namespace detail {
      inline void outputMsgOSS( std::ostringstream&& msg_ss, MsgType mt )
      {
        outputMsg( msg_ss.str().c_str(), mt );
      }
    }

    //Set msg handler. An empty function indicates the default behaviour.
    using MsgHandlerFct_t = std::function<void(const char*,MsgType)>;
    void setMsgHandler( MsgHandlerFct_t = nullptr );
  }

#ifdef NCRYSTAL_MSG
#  undef NCRYSTAL_MSG
#endif
#ifdef NCRYSTAL_WARN
#  undef NCRYSTAL_WARN
#endif
#ifdef NCRYSTAL_RAWOUT
#  undef NCRYSTAL_RAWOUT
#endif

#define NCRYSTAL_MSG(msg) ::NCRYSTAL_NAMESPACE::Msg::detail:: \
  outputMsgOSS( std::ostringstream() << msg, \
                ::NCRYSTAL_NAMESPACE::MsgType::Info );
#define NCRYSTAL_WARN(msg) ::NCRYSTAL_NAMESPACE::Msg::detail:: \
  outputMsgOSS( std::ostringstream() << msg, \
                ::NCRYSTAL_NAMESPACE::MsgType::Warning );
#define NCRYSTAL_RAWOUT(msg) ::NCRYSTAL_NAMESPACE::Msg::detail:: \
  outputMsgOSS( std::ostringstream() << msg, \
                ::NCRYSTAL_NAMESPACE::MsgType::RawOutput );

}


#endif

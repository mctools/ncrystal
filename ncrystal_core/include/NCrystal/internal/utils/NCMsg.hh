#ifndef NCrystal_Msg_hh
#define NCrystal_Msg_hh

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
#include <sstream>

namespace NCRYSTAL_NAMESPACE {
  namespace Msg {

    //All NCrystal code should always use the following functions (or the
    //macro's below) in favour of direct std::cout usage. This not only ensures
    //proper capture of all output, but also makes output thread-safe due to
    //projection by a mutex).

    void outputMsg( const char * msg, MsgType mt = MsgType::Info );
    void outputMsg( const std::string& msg, MsgType mt = MsgType::Info );

    //Set msg handler. An empty function indicates the default behaviour.
    using MsgHandlerFct_t = std::function<void(const char*,MsgType)>;
    void setMsgHandler( MsgHandlerFct_t = nullptr );

#ifdef NCRYSTAL_MSG
#  undef NCRYSTAL_MSG
#endif
#define NCRYSTAL_MSG(msg) ::NCRYSTAL_NAMESPACE::Msg::detail:: \
  outputMsgMS( ::NCRYSTAL_NAMESPACE::Msg::detail::MsgStream() << msg, \
               ::NCRYSTAL_NAMESPACE::MsgType::Info );

#ifdef NCRYSTAL_WARN
#  undef NCRYSTAL_WARN
#endif
#define NCRYSTAL_WARN(msg) ::NCRYSTAL_NAMESPACE::Msg::detail:: \
  outputMsgMS( ::NCRYSTAL_NAMESPACE::Msg::detail::MsgStream() << msg, \
               ::NCRYSTAL_NAMESPACE::MsgType::Warning );

#ifdef NCRYSTAL_RAWOUT
#  undef NCRYSTAL_RAWOUT
#endif
#define NCRYSTAL_RAWOUT(msg) ::NCRYSTAL_NAMESPACE::Msg::detail:: \
  outputMsgMS( ::NCRYSTAL_NAMESPACE::Msg::detail::MsgStream() << msg, \
               ::NCRYSTAL_NAMESPACE::MsgType::RawOutput );

#ifdef NCPLUGIN_NAME_CSTR
    //The NCPLUGIN_MSG and NCPLUGIN_WARN macros are only available if you are in
    //a plugin using NCPluginBoilerplate.hh, and should be preferred in such
    //plugin code over NCRYSTAL_MSG and NCRYSTAL_WARN:
#  ifdef NCPLUGIN_MSG
#    undef NCPLUGIN_MSG
#  endif
#  ifdef NCPLUGIN_WARN
#    undef NCPLUGIN_WARN
#  endif
#  define NCPLUGIN_MSG(msg) ::NCRYSTAL_NAMESPACE::Msg::detail:: \
    outputMsgMS( ::NCRYSTAL_NAMESPACE::Msg::detail::MsgStream() \
                   <<"plugin::" NCPLUGIN_NAME_CSTR ": "  << msg, \
                 ::NCRYSTAL_NAMESPACE::MsgType::Info );
#  define NCPLUGIN_WARN(msg) ::NCRYSTAL_NAMESPACE::Msg::detail:: \
    outputMsgMS( ::NCRYSTAL_NAMESPACE::Msg::detail::MsgStream() \
                   <<"plugin::" NCPLUGIN_NAME_CSTR ": "  << msg, \
                 ::NCRYSTAL_NAMESPACE::MsgType::Warning );
#endif

  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace Msg {
    //Note that some of the solutions here were inspired by discussions at
    //https://stackoverflow.com/questions/303562.
    namespace detail {
      void outputMsgImpl( const char * msg, MsgType mt );
      struct MsgStream {
        std::ostringstream thestream;
        template<class T>
        MsgStream& operator<<(T const& v) { thestream << v; return *this; }
      };
      inline void outputMsgMS( MsgStream& ms, MsgType mt )
      {
        outputMsgImpl( ms.thestream.str().c_str(), mt );
      }
    }
    inline void outputMsg( const char * msg, MsgType mt )
    {
      detail::outputMsgImpl( msg, mt );
    }
    inline void outputMsg( const std::string& msg, MsgType mt )
    {
      detail::outputMsgImpl( msg.c_str(), mt );
    }
  }
}

#endif

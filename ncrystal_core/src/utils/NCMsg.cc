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

#include "NCrystal/internal/utils/NCMsg.hh"
#include <iostream>

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace Msg {
    namespace detail {
      namespace {
        struct MsgHandler {
          std::mutex mutex;
          MsgHandlerFct_t handler;
        };
        MsgHandler& getMsgHandler()
        {
          static MsgHandler mh;
          return mh;
        }
      }
    }
  }
}

void NC::Msg::setMsgHandler( MsgHandlerFct_t fct )
{
  auto& mh = detail::getMsgHandler();
  NCRYSTAL_LOCK_GUARD(mh.mutex);
  mh.handler = std::move(fct);
}

namespace NCRYSTAL_NAMESPACE {
  namespace Msg {
    namespace {
      void defaultMsgHandler( const char * msg, MsgType mt )
      {
        switch( mt ) {
        case MsgType::RawOutput:
          std::cout << msg << std::flush;//NB: Do not append newline!
          break;
        case MsgType::Info:
          std::cout << "NCrystal: " << msg << std::endl;
          break;
        case MsgType::Warning:
          std::cout << "NCrystal WARNING: " << msg << std::endl;
          break;
        default:
          nc_assert_always(false);
        };
      }
    }
  }
}

void NC::Msg::detail::outputMsgImpl( const char * msg, MsgType mt )
{
  auto& mh = getMsgHandler();
  NCRYSTAL_LOCK_GUARD(mh.mutex);//Protect in case of multithreading.
  if ( mh.handler )
    mh.handler( msg, mt );
  else
    defaultMsgHandler( msg, mt );
}

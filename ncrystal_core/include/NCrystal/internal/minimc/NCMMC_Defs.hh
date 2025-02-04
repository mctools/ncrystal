#ifndef NCrystal_MMC_Defs_hh
#define NCrystal_MMC_Defs_hh

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/interfaces/NCProcImpl.hh"

namespace NCRYSTAL_NAMESPACE {


  namespace MiniMC {

    //Basket size (note, if we migrate to dynamic sizes, we could increase the
    //maximum size to e.g. 65536). It should be a power of two.
    static constexpr size_t basket_N = 4096;
    static constexpr size_t basket_N_almost_Full = basket_N*7/8;

    using ProcImpl::ProcPtr;
    using ProcImpl::Process;

    inline constexpr double macroXS( NumberDensity nd, CrossSect xs ) noexcept
    {
      //returns macroscocopic XS, or attenuation coefficient, in units of
      //1/m. NB, using nc_as_const for pre C++17:
      return 100.0 * nc_as_const(nd).dbl() * nc_as_const(xs).dbl();
    }

  }
}

#include "NCrystal/internal/utils/NCMsg.hh"
#ifdef NCRYSTAL_DEBUGMMCMSG
#  undef NCRYSTAL_DEBUGMMCMSG
#endif

#if 0 //#ifndef NDEBUG (always disable, messes up tests)
#define NCRYSTAL_DEBUGMMCMSG(msg) ::NCRYSTAL_NAMESPACE::Msg::detail:: \
  outputMsgMS( ::NCRYSTAL_NAMESPACE::Msg::detail::MsgStream() << "MMC:" << msg, \
               ::NCRYSTAL_NAMESPACE::MsgType::Info );
#else
#  define NCRYSTAL_DEBUGMMCMSG(msg) {}
#endif

#endif

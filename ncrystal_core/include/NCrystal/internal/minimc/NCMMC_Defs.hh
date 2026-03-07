#ifndef NCrystal_MMC_Defs_hh
#define NCrystal_MMC_Defs_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

    //Basket size (should be power of 2).
    static constexpr std::size_t basket_N = 4096;
    static constexpr std::size_t basket_N_almost_Full = basket_N*7/8;

    //The basket type determines what will be available for tallies and
    //callbacks:
    enum class BasketType : unsigned {
      Invalid = 0,
      Basic = 1,//Neutron parameters + nscat/inelas
      Extended = 2//Basic + initial neutron + (fixme?) neutron id
    };

    //Function signature for acceptance of tally results:
    class Basket;
    using TallyFct = std::function<void(const Basket&)>;

    struct ParticleCountSum {
      //Number and total weight of particles.
      std::size_t count = 0;
      double weight = 0.0;
    };

    struct BasketValBufDbl final {
      using data_t = double[basket_N];
      data_t data;
      constexpr const double& operator[]( std::size_t ) const ncnoexceptndebug;
      //NB: No mutable operator[] since it causes cppcheck false positives
    };

    struct BasketValBufInt final {
      using data_t = int[basket_N];
      data_t data;
      constexpr const int& operator[]( std::size_t ) const ncnoexceptndebug;
      //NB: No mutable operator[] since it causes cppcheck false positives
    };

    //For efficiency, handle larger number of neutrons at once, with each field
    //in a separate array. This has several advantages, mainly to reduce
    //overhead of context switches and branches, as well is making it easier to
    //simd vectorisation, hot caches, etc.

    struct NeutronBasketFields final : NoCopyMove {
      BasketValBufDbl x, y, z, ux, uy, uz, ekin, w;
    };

    struct NeutronBasket final : NoCopyMove {
      NeutronBasketFields fields;
      std::size_t nused = 0;
      constexpr bool full() const noexcept { return nused==basket_N; }
      constexpr bool empty() const noexcept { return nused==0; }
      constexpr std::size_t size() const noexcept { return nused; }
    };

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

////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    inline constexpr const double&
    BasketValBufDbl::operator[]( std::size_t i ) const ncnoexceptndebug
    {
      return nc_assert_rv(i<basket_N), data[i];
    }

    inline constexpr const int&
    BasketValBufInt::operator[]( std::size_t i ) const ncnoexceptndebug
    {
      return nc_assert_rv(i<basket_N), data[i];
    }
  }
}

#endif

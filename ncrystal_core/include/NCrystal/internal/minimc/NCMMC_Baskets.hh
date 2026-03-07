#ifndef NCrystal_MMC_Baskets_hh
#define NCrystal_MMC_Baskets_hh

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

#include "NCrystal/internal/minimc/NCMMC_Geom.hh"
#include "NCrystal/internal/minimc/NCMMC_Source.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    //Non-templated abstraction layer for the
    //MemPool/Basket/BasketHolder/BasketMgr infrastructure. This allows
    //simulation engines to be non-templated and decoupled from the actual
    //basket internals. This abstraction layer should be cheap since we are
    //dealing with many neutrons at once in each call.

    namespace detail { class UBImpl; }

    class UniversalBasket final : MoveOnly {
      //Abstraction of a neutron basket, allowing for multiple implementations
      //of the memory layout etc. behind the scenes, without coupling the basket
      //and client implementations.
      void * internal = nullptr;
      friend class detail::UBImpl;
      void dealloc_warn() noexcept;
    public:
      UniversalBasket() = default;
      ~UniversalBasket();
      UniversalBasket( UniversalBasket&& ) noexcept;
      UniversalBasket& operator=( UniversalBasket&& ) noexcept;
      void swap( UniversalBasket& ) noexcept;
      BasketType basketType() const noexcept;
      std::size_t size() const ncnoexceptndebug;
      std::size_t empty() const ncnoexceptndebug { return size()==0; }
      bool valid() const noexcept { return neutrons != nullptr; }

      //This one is always set for any valid view:
      NeutronBasket * neutrons = nullptr;

      //These might be absent (fixme: guarantee all basic fields?):
      BasketValBufInt* nscat = nullptr;
      BasketValBufInt* nscat_inelas = nullptr;
      BasketValBufDbl* buf1 = nullptr;//generic buffer available to engine (if
                                      //available, values are zero initialised).
      const NeutronBasketFields * neutrons_initial = nullptr;

      //Append single neutron from other compatible basket. Returns idx in this
      //basket of the appended neutron.
      std::size_t append1( const UniversalBasket& o, std::size_t idx_o );
    };

    class UniversalBasketMgr : NoCopyMove {
    public:
      //Universal abstraction of various basket management operations.
      virtual ~UniversalBasketMgr() = default;

      virtual UniversalBasket allocateBasket() = 0;
      virtual void deallocateBasket( UniversalBasket&& b ) = 0;
      virtual void addPendingBasket( UniversalBasket&& b ) = 0;

      //Each MT worker should have its own UniversalBasketMgr, since some memory
      //pools are kept thread-local for efficiency. Before workers are spawned,
      //additional instances should therefore be created with this:
      virtual shared_obj<UniversalBasketMgr> cloneMgrForThread() = 0;
    };

    class InputBasketProvider : NoCopyMove {
    public:
      //Manager providing input baskets for simulation steps, by combining
      //particles from the source with those added back during a simulation step
      //by a call to addPendingBasket.

      virtual ~InputBasketProvider() = default;

      //Get the input baskets. Any source neutrons missing the geometry will be
      //added to the counts in the missCount object, and provided to the
      //TallyFct (if one is provided).
      virtual UniversalBasket getInputBasket( RNG&, const TallyFct&,
                                              ParticleCountSum& missCount ) = 0;

      //In case simulation needs to be halted prematurely for whatever reason, a
      //call to haltSource will cause no more neutrons to be taken from the
      //source:
      virtual void haltSource() = 0;

    };

    //The two closely related managers are always created together by this
    //factory function. Note that the UniversalBasketMgr should be subsequently
    //cloned so each worker thread has its own instance, while the same
    //InputBasketProvider is used by all threads. The implementation details can
    //also depend on the thread count of the simulation:
    using BasketManagementPair = std::pair< shared_obj<UniversalBasketMgr>,
                                            shared_obj<InputBasketProvider> >;
    BasketManagementPair createBasketManagement( GeometryPtr,
                                                 SourcePtr,
                                                 BasketType,
                                                 ThreadCount );

  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    inline UniversalBasket::~UniversalBasket() {
      if ( internal )
        dealloc_warn();
    }

    inline UniversalBasket::UniversalBasket( UniversalBasket&& o ) noexcept
    {
      swap(o);
    }

    inline UniversalBasket& UniversalBasket::operator=( UniversalBasket&& o ) noexcept
    {
      swap(o);//self-assignment is ok, just swapping ptrs with themselves.
      return *this;
    }

    inline void UniversalBasket::swap( UniversalBasket& o ) noexcept
    {
      std::swap(internal,o.internal);
      std::swap(neutrons,o.neutrons);
      std::swap(nscat,o.nscat);
      std::swap(nscat_inelas,o.nscat_inelas);
      std::swap(neutrons_initial,o.neutrons_initial);
      std::swap(buf1,o.buf1);
    }

    inline BasketType UniversalBasket::basketType() const noexcept
    {
      return ( internal ? ( neutrons_initial==nullptr
                            ? BasketType::Basic
                            : BasketType::Extended )
               : BasketType::Invalid );
    }

    inline std::size_t UniversalBasket::size() const ncnoexceptndebug
    {
      nc_assert(valid());
      return neutrons->size();
    }
  }
}

#endif

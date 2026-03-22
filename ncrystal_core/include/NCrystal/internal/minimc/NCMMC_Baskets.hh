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

    namespace detail { class UBImpl; }

    class Basket final : MoveOnly {
      //Abstraction of a neutron basket, allowing for multiple implementations
      //of the memory layout etc. behind the scenes, without coupling the basket
      //and client implementations.
      void * internal = nullptr;
      friend class detail::UBImpl;
      void dealloc_warn() noexcept;
    public:
      Basket() = default;
      ~Basket();
      Basket( Basket&& ) noexcept;
      Basket& operator=( Basket&& ) noexcept;
      void swap( Basket& ) noexcept;
      BasketType basketType() const noexcept;
      std::size_t size() const ncnoexceptndebug;
      std::size_t empty() const ncnoexceptndebug { return size()==0; }
      bool valid() const noexcept { return neutrons != nullptr; }

      //This is always available:
      NeutronBasket * neutrons = nullptr;
      BasketValBufInt* nscat = nullptr;
      BasketValBufInt* nscat_inelas = nullptr;
      BasketValBufDbl* buf1 = nullptr;//generic buffer available to engine
                                      //(values are zero initialised).

      //This is only available for Extended baskets:
      const NeutronBasketFields * neutrons_initial = nullptr;

      //Append single neutron from other compatible basket. Returns idx in this
      //basket of the appended neutron.
      std::size_t append1( const Basket& o, std::size_t idx_o );
    };

    class WorkerToken;

    class BasketMgr : NoCopyMove {
    public:

      //Manager class handling basket allocation's, and shipping baskets of
      //neutrons from the source and between worker threads.

      virtual ~BasketMgr() = default;

      //Get an input basket to be processed by the simulation engine. Will
      //either return a non-empty valid basket, or an invalid one at the end of
      //simulation. In case no baskets are pending and the source has run out,
      //this might mean waiting to see if other worker threads will register new
      //pending baskets.
      //
      //Any neutrons delivered by the source and which are missing the geometry,
      //will be added to the counts in the missCount object, and provided to the
      //TallyFct (if one is provided).
      //
      //For the manager to do its work, it is required that each worker thread
      //should keep a single worker token around, with lifetime ending when the
      //thread is done (after having received an invalid basket).
      virtual Basket getInputBasket( WorkerToken&,
                                     RNG&,
                                     const TallyFct&,
                                     ParticleCountSum& missCount ) = 0;

      //Add a basket needing further processing:
      virtual void addPendingBasket( Basket&& ) = 0;

      //Apart from i/o to the baskets needing processing, basic basket
      //allocation and deallocation is also possible:
      virtual Basket allocateBasket() = 0;
      virtual void deallocateBasket( Basket&& ) = 0;


      //In case simulation needs to be halted prematurely for whatever reason, a
      //call to haltSource will cause no more neutrons to be taken from the
      //source (but already pending neutrons will be continued to be processed).
      virtual void haltSource() = 0;

      //If the simulation is being aborted due to an exception or similar error
      //in one thread, a call to haltError will completely block more baskets
      //coming out of getInputBasket:
      virtual void haltError() = 0;

      //Dispose of a token (done automatically from the token destructor):
      virtual void tokenDispose( WorkerToken& ) = 0;

    };

    //Factory function for the manager:
    shared_obj<BasketMgr> createBasketMgr( GeometryPtr,
                                           SourcePtr,
                                           BasketType );


    class WorkerToken final : MoveOnly {
      std::shared_ptr<BasketMgr> m_bmgr;
      bool m_active = false;
      Basket m_basketbuf;
      friend class detail::UBImpl;
      bool needsDispose() const noexcept
      {
        return m_active || m_basketbuf.valid();
      };
    public:
      WorkerToken( shared_obj<BasketMgr> ibp )
        : m_bmgr( ibp.getsp() )
      {
      }
      void dispose()
      {
        if ( needsDispose() )
          m_bmgr->tokenDispose(*this);
        nc_assert( !needsDispose() );
      }
      ~WorkerToken() {
        if ( needsDispose() ) {
          try {
            dispose();
          } catch (...) {
            std::printf("NCrystal ERROR: Ignoring Exception"
                        " in WorkerToken destructor!\n");
          }
        }
      }
    };
  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    inline Basket::~Basket() {
      if ( internal )
        dealloc_warn();
    }

    inline Basket::Basket( Basket&& o ) noexcept
    {
      swap(o);
    }

    inline Basket& Basket::operator=( Basket&& o ) noexcept
    {
      swap(o);//self-assignment is ok, just swapping ptrs with themselves.
      return *this;
    }

    inline void Basket::swap( Basket& o ) noexcept
    {
      std::swap(internal,o.internal);
      std::swap(neutrons,o.neutrons);
      std::swap(nscat,o.nscat);
      std::swap(nscat_inelas,o.nscat_inelas);
      std::swap(neutrons_initial,o.neutrons_initial);
      std::swap(buf1,o.buf1);
    }

    inline BasketType Basket::basketType() const noexcept
    {
      return ( internal ? ( neutrons_initial==nullptr
                            ? BasketType::Basic
                            : BasketType::Extended )
               : BasketType::Invalid );
    }

    inline std::size_t Basket::size() const ncnoexceptndebug
    {
      nc_assert(valid());
      return neutrons->size();
    }
  }
}

#endif

#ifndef NCrystal_MMC_BasketMgr_hh
#define NCrystal_MMC_BasketMgr_hh

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
#include "NCrystal/core/NCSmallVector.hh"
#include "NCrystal/internal/minimc/NCMMC_Defs.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    struct HeapMem_ensure_alloc_t {};

    template<unsigned N_ALIGN, unsigned N_SIZE>
    class HeapMem : MoveOnly {
    public:
      using ensure_alloc_t = HeapMem_ensure_alloc_t;
      HeapMem() = default;
      HeapMem( HeapMem o, ensure_alloc_t )
        : m_heapptr( std::move(o.m_heapptr) )
      {
        if (!m_heapptr.data)
          allocate();
      }
      HeapMem( ensure_alloc_t )
      {
        if (!m_heapptr.data)
          allocate();
      }
#ifndef NDEBUG
      ~HeapMem()
      {
        assert(!m_heapptr.data);//test against inadvertant deallocation (the
                                //HeapMemPool will deallocate objects before
                                //letting them go out of scope).
      }
#endif

      //Important (for usage in SmallVector) that we are
      //nothrow-move-assignable. The downside is that HeapMem( std::move(o) )
      //does not guarantee to leave o empty.
      HeapMem( HeapMem&& o ) noexcept { m_heapptr.swap( o.m_heapptr ); }
      HeapMem& operator=( HeapMem&& o ) noexcept { m_heapptr.swap( o.m_heapptr ); return *this; }

      void * data() noexcept { return m_heapptr.data; }
      const void * data() const noexcept { return m_heapptr.data; }
      void allocate() { m_heapptr.allocate( N_SIZE ); }
      void deallocate() { m_heapptr.deallocate(); }
      void swap( HeapMem& o ) noexcept { m_heapptr.swap( o.m_heapptr ); }
    private:
      AlignedAlloc::AlignedHeapPtr<N_ALIGN> m_heapptr;
    };

    template<class TBasket>
    class BasketHolder : MoveOnly {
    public:
      static_assert(std::is_trivially_destructible<TBasket>::value,"");
      static_assert(std::is_nothrow_destructible<TBasket>::value,"");
      static_assert(std::is_default_constructible<TBasket>::value,"");
      static_assert(std::is_standard_layout<TBasket>::value,"");
      //      static_assert(std::is_trivially_constructible<TBasket>::value,"");
      static_assert(std::is_nothrow_constructible<TBasket>::value,"");
      using heapmem_t = HeapMem<alignof(TBasket),sizeof(TBasket)>;

      //NB: TBasket is trivially destructible, so OK that we never invoke the
      //destructor.

      BasketHolder( heapmem_t hm )
        : m_heapmem( std::move(hm), typename heapmem_t::ensure_alloc_t() ),
          m_basket(new(m_heapmem.data()) TBasket)
      {
      }

      BasketHolder()
      {
        m_heapmem.allocate();
        //NB: TBasket is trivially destructible, so OK that we never invoke the
        //destructor.
        m_basket = new(m_heapmem.data()) TBasket;
      }

      BasketHolder( no_init_t )
        : m_basket(nullptr)
      {
      }

      // Not used:
      // void deallocate()
      // {
      //   m_basket = nullptr;
      //   m_heapmem.deallocate();
      // }

      void swap( BasketHolder& o ) noexcept
      {
        std::swap(m_heapmem,o.m_heapmem);
        std::swap(m_basket,o.m_basket);
      }

      //Boolean operator can be used to check for invalid object:
      constexpr bool valid() const noexcept { return m_basket!=nullptr; }
      constexpr bool operator()() const noexcept { return m_basket!=nullptr; }
      constexpr bool operator!() const noexcept { return m_basket==nullptr; }

      //Only access basket for valid objects:
      TBasket& basket() noexcept { return *m_basket; }
      const TBasket& basket() const noexcept { return *m_basket; }

      //To use by the memory pool manager (do not access basket() afterwards):
      heapmem_t stealMemory()
      {
        m_basket = nullptr;
        return std::move(m_heapmem);
      }

      //Using swap rather than move to be noexcept:
      BasketHolder( BasketHolder&& o ) noexcept
      {
        this->swap( o );
      }

      BasketHolder& operator=( BasketHolder&& o ) noexcept
      {
        this->swap( o );
        return *this;
      }

    public:
      heapmem_t m_heapmem;
      TBasket * m_basket = nullptr;
    };

    template<class THeapMem, unsigned N_MAX_KEEP = 16 >
    class HeapMemPool : NoCopyMove {
    public:
      using heapmem_t = THeapMem;
      //NB: Not thread-safe!
      bool empty() const { return m_pool.empty(); }
      std::size_t size() const { return m_pool.size(); }
      heapmem_t allocate()
      {
        if ( m_pool.empty() )
          return heapmem_t(typename heapmem_t::ensure_alloc_t());
        auto o = heapmem_t( std::move(m_pool.back()),
                            typename heapmem_t::ensure_alloc_t() );
        m_pool.pop_back();
        return o;
      }
      void deallocate( heapmem_t&& o )
      {
        if ( o.data() && m_pool.size() < N_MAX_KEEP ) {
          // m_pool.reserve(N_MAX_KEEP);
          m_pool.push_back( std::move(o) );
        }
      }

#ifndef NDEBUG
      ~HeapMemPool()
      {
        //Explicitly clear mempool objects, so we can detect if any such objects
        //have data when their destructors are reached (if so, they must have
        //been destructed outside of the pool rather than being first returned
        //to it):
        for ( auto& e : m_pool )
          e.deallocate();
      }
#endif
    private:
      SmallVector<heapmem_t,N_MAX_KEEP> m_pool;
    };

    template<class TBasket>
    class BasketMgr : NoCopyMove {
    public:
      using basket_t = TBasket;
      using basket_holder_t = BasketHolder<TBasket>;
      static_assert(std::is_nothrow_move_constructible<basket_holder_t>::value,"");
      using heapmem_t = typename basket_holder_t::heapmem_t;
      using heapmempool_t = HeapMemPool<heapmem_t,16>;

      basket_holder_t allocateBasket()
      {
        NCRYSTAL_LOCK_GUARD(m_mutex);
        return { m_mempool.allocate() };
      }

      void addPendingBasket( basket_holder_t&& o )
      {
        nc_assert(o.valid());
        NCRYSTAL_LOCK_GUARD(m_mutex);
        NCRYSTAL_DEBUGMMCMSG("Got pending basket with size "<<o.basket().size());
        if ( o.basket().empty() )
          m_mempool.deallocate( o.stealMemory() );
        else
          m_pending_baskets.push_back( std::move(o) );
      }

      //Will return invalid basket_holder_t object if no more pending baskets
      //are available.
      basket_holder_t getPendingBasket()
      {
        basket_holder_t bh{ no_init };
        NCRYSTAL_LOCK_GUARD(m_mutex);
        if ( !m_pending_baskets.empty() ) {
          bh = std::move( m_pending_baskets.back() );
          m_pending_baskets.pop_back();
        }
        return bh;
      }

      basket_holder_t getPendingBasketOrAllocateEmpty( ThreadCount nthreads )
      {
        SmallVector<basket_holder_t,8> bhs_to_merge;
        {
          //Return one basket or move baskets to merge from m_pending_baskets to
          //bhs_to_merge, all while the mutex is locked:
          NCRYSTAL_LOCK_GUARD(m_mutex);
          if ( m_pending_baskets.empty() )
            return { m_mempool.allocate() };

          if ( m_pending_baskets.size()+1 <= static_cast<std::size_t>(nthreads.get()) ) {
            //We have more worker threads than remaining baskets (+1 to avoid
            //spuriously triggering this). So just take one easily accessible
            //basket (better to have N threads with 1000 neutrons each, than N-1
            //idling threads and 1 thread with N*1000 events). The +1 is to
            //avoid triggering this too early or when nthreads=1. However, we do
            //not do this if it would result in ridiculously small amount of
            //particles in the returned basket, since then other per-basket
            //overhead might dominate.
            nc_assert(!m_pending_baskets.empty());
            if ( m_pending_baskets.back().basket().size() >= 32 ) {
              basket_holder_t res = std::move(m_pending_baskets.back());
              m_pending_baskets.pop_back();
              return res;
            }
          }

          //TODO: Now we gave up on the FIFO requirement, we could simplify +
          //optimise the code below somewhat.

          std::size_t ntot = 0;
          auto it = m_pending_baskets.begin();
          auto itE = m_pending_baskets.end();
          while( it!=itE && ntot + it->basket().size() <= basket_N ) {
            ntot += it->basket().size();
            bhs_to_merge.push_back(std::move(*it));
            ++it;
          }
          nc_assert_always(!bhs_to_merge.empty());
          auto nconsumed = bhs_to_merge.size();
          nc_assert_always( nconsumed >= 1 );
          it = m_pending_baskets.begin();
          //itE_consumed = std::next(it,nconsumed);
          for ( it = std::next(m_pending_baskets.begin(),nconsumed);
                it != itE; ++it )
            *std::prev(it,nconsumed) = std::move(*it);
          m_pending_baskets.resize( m_pending_baskets.size() - nconsumed );
          if ( !m_pending_baskets.empty() && ntot < basket_N / 2 ) {
            //TODO: If we have small baskets in front of an almost full basket,
            //we might end up with very few particles to return. So already here
            //(under the lock) we move over some particles from that first
            //basket:
            nc_assert_always(!bhs_to_merge.empty());
            auto& tgt_basket = bhs_to_merge.back().basket();
            auto& src_basket = m_pending_baskets.front().basket();
            std::size_t n_move = basket_N/2 - tgt_basket.size();
            NCRYSTAL_DEBUGMMCMSG("Split-off "<<n_move<<" neutrons from first remaining basket");
            nc_assert_always( src_basket.size() > n_move );//won't become empty
            nc_assert_always( tgt_basket.size()+n_move <= basket_N );
            tgt_basket.appendEntriesFromOther( src_basket, src_basket.size()-n_move, n_move );
            src_basket.neutrons.nused -= n_move;
            nc_assert_always(src_basket.size()>0);
          }
        }//release lock!

        //Merging is performed without a lock:
        nc_assert_always(!bhs_to_merge.empty());

        if (bhs_to_merge.size()==1)
          return std::move(bhs_to_merge.front());
        //First sort, so biggest basket goes first (to copy the least particle
        //states):
        std::stable_sort(bhs_to_merge.begin(),
                         bhs_to_merge.end(),
                         [](const basket_holder_t &a, const basket_holder_t &b)
                         {
                           return a.basket().size() > b.basket().size();
                         });
        //Sanity check that the sort was the correct way:
        nc_assert( bhs_to_merge.front().basket().size()
                   >= bhs_to_merge.back().basket().size() );

        auto itB = bhs_to_merge.begin();
        auto itE = bhs_to_merge.end();
        basket_holder_t bh( std::move(*itB) );
        ++itB;
        for (auto it = itB; it!=itE; ++it ) {
          std::size_t ntransfer = std::min<std::size_t>(basket_N - bh.basket().size(),
                                                        it->basket().size());
          nc_assert(ntransfer>0);//otherwise list was build wrongly.
          const std::size_t o_new_size = it->basket().size() - ntransfer;
          bh.basket().appendEntriesFromOther(  it->basket(), o_new_size, ntransfer );
          it->basket().neutrons.nused = o_new_size;//"pop off" copied entries
        }
        //Ready to return bh, but first properly deal with other baskets in
        //[itB,itE] range:
        {
          NCRYSTAL_LOCK_GUARD(m_mutex);
          unsigned n_returned(0);
          for (auto it = itB; it!=itE; ++it ) {
            nc_assert(it->valid());
            if ( it->basket().empty() ) {
              m_mempool.deallocate( it->stealMemory() );
            } else {
              m_pending_baskets.push_back( std::move(*it) );
              ++n_returned;
            }
          }
          nc_assert_always(n_returned<=1);
        }
        return bh;
      }

      //Please return baskets created with getPendingBasket(..) here once done
      //with them:
      void deallocateBasket( basket_holder_t&& bh )
      {
        NCRYSTAL_LOCK_GUARD(m_mutex);
        m_mempool.deallocate( bh.stealMemory() );
      }

    private:
      heapmempool_t m_mempool;
      std::vector<basket_holder_t> m_pending_baskets;
      std::mutex m_mutex;//for m_mempool and m_pending_baskets
    };
  }
}

#endif

#ifndef NCrystal_MMC_MemPool_hh
#define NCrystal_MMC_MemPool_hh

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

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/core/NCSmallVector.hh"
#include "NCrystal/internal/minimc/NCMMC_Defs.hh"

// Memory pool infrastructure, intended for managing heap allocations for the
// baskets.

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
        assert( o.m_heapptr.data == nullptr );
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
        //We do not deallocate here, since calling code should have done that
        //without letting allocated HeapMem objects go out of scope. But we add
        //an assert to detect implementation flaws.
        assert(!m_heapptr.data);
      }
#endif

      //Important (for usage in SmallVector) that we are
      //nothrow-move-assignable. The downside is that HeapMem( std::move(o) )
      //does not guarantee to leave o empty, so we have the release() method for
      //that instead.
      HeapMem( HeapMem&& o ) noexcept { m_heapptr.swap( o.m_heapptr ); }
      HeapMem& operator=( HeapMem&& o ) noexcept { m_heapptr.swap( o.m_heapptr ); return *this; }

      void * data() noexcept { return m_heapptr.data; }
      const void * data() const noexcept { return m_heapptr.data; }
      void allocate() { m_heapptr.allocate( N_SIZE ); }
      void deallocate() { if (m_heapptr.data) m_heapptr.deallocate(); nc_assert(!m_heapptr.data); }
      void swap( HeapMem& o ) noexcept { m_heapptr.swap( o.m_heapptr ); }
      HeapMem release() noexcept { HeapMem h; this->swap(h); return h; }
      //RAII breaking methods:
      void* release_raw() { return m_heapptr.release_raw(); }
      void set_raw(void*data) ncnoexceptndebug
      { nc_assert(!m_heapptr.data); m_heapptr.set_raw(data); }
    private:
      AlignedAlloc::AlignedHeapPtr<N_ALIGN> m_heapptr;
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
        heapmem_t res;
        nc_assert(res.data()==nullptr);
        if ( m_pool.empty() ) {
          res.allocate();
          nc_assert(res.data());
        } else {
          res = m_pool.back().release();
          nc_assert(m_pool.back().data()==nullptr);
          nc_assert(res.data()!=nullptr);
          m_pool.pop_back();
        }
        nc_assert(res.data());
        return res;
      }
      void deallocate( heapmem_t&& o )
      {
        if ( !o.data() )
          return;
        if ( m_pool.size() < N_MAX_KEEP ) {
          m_pool.emplace_back( o.release() );
          nc_assert( o.data() == nullptr );
          nc_assert( m_pool.back().data() != nullptr );
        } else {
          o.deallocate();
        }
      }

      void clearOutPoolDueToErrors()
      {
        //should normally only be used in case of error aborts.
        for ( auto& e : m_pool )
          e.deallocate();
        m_pool.clear();
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

  }
}

#endif

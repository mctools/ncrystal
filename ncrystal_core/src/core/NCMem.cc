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

#include "NCrystal/core/NCMem.hh"
#include "NCrystal/core/NCDefs.hh"
#include <vector>
#include <mutex>

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    static std::mutex s_cacheCleanerMutex;
    static std::vector<voidfct_t> s_cacheCleanerMutexFcts;
  }
}

void NC::clearCaches()
{
  NCRYSTAL_LOCK_GUARD(s_cacheCleanerMutex);
  for (auto& f : s_cacheCleanerMutexFcts)
    f();
}

void NC::registerCacheCleanupFunction( voidfct_t f )
{
  NCRYSTAL_LOCK_GUARD(s_cacheCleanerMutex);
  s_cacheCleanerMutexFcts.emplace_back(f);
}

void * NC::AlignedAlloc::detail::bigAlignedAlloc( std::size_t alignment, std::size_t nbytes )
{
  nc_assert( alignment > detail::nc_alignof_max_align_t );
  nc_assert( nbytes > 0 );

  //The hard but most portable and standards complient way! Get the
  //alignment through over-allocation, and even overallocate enough that
  //we can have an (aligned) pointer value in front of the returned
  //data. In that pointer value we store the address we must eventually
  //pass to std::free (i.e. the actual result from std::malloc).

  //Alignment big enough for both a pointer and the actual data (NB: use gcd if
  //we stop requiring power of 2 alignments):
  const std::size_t real_align_size = ( alignment < alignof(void*) ? alignof(void*) : alignment );
  //std::malloc with place for nbytes, the ptr, and an overallocation of
  //real_align_size.
  std::size_t nalloc = nbytes + sizeof(void*) + real_align_size;
  void * ptr_allocated = nc_std_malloc( nalloc );
  auto find_first_aligned_address_after_ptr = []( void * data, std::size_t align )
  {
    //Ignore initial unaligned part (I forgot where I found this magic, most
    //likely a small example on stackoverflow ):
    std::uintptr_t pn = reinterpret_cast< std::uintptr_t >( data );
    std::uintptr_t pn_first_aligned_address = ( ( pn + align - 1 ) & - align );
    assert( pn_first_aligned_address % align == 0 );
#ifndef NDEBUG
    if ( reinterpret_cast< std::uintptr_t >( data ) % align == 0 ) {
      //already aligned, do not unintentionally move forward:
      assert( reinterpret_cast<void*>(pn_first_aligned_address) == data );
    }
#endif
    return reinterpret_cast<void*>(pn_first_aligned_address);
  };
  auto check_align_ok = []( void * ptr, std::size_t align )
  {
    (void)ptr;
    (void)align;
    assert( ptr );
    assert( reinterpret_cast< std::uintptr_t >( ptr ) % align == 0 );
  };
  void* ptr1 = find_first_aligned_address_after_ptr( ptr_allocated, real_align_size );
  check_align_ok( ptr1, real_align_size );
  check_align_ok( ptr1, alignment );
  while ( (char*)ptr_allocated+sizeof(void*) > (char*)ptr1 )
    ptr1 = reinterpret_cast<void*>( reinterpret_cast<std::uintptr_t>(ptr1) + real_align_size );
  void * ptr0 = (void*)((char*)ptr1 - sizeof(void*));
  check_align_ok( ptr0, alignof(void*) );
  check_align_ok( ptr1, alignment );
  void * data = ptr1;
  assert( ptr0 >= ptr_allocated );
  assert( ptr1 >= ptr_allocated );
  *reinterpret_cast<void**>(ptr0) = ptr_allocated;
  assert ( (char*)data + nbytes <= (char*)ptr_allocated + nalloc );
  assert ( (char*)ptr0 + sizeof(void*) == (char*)ptr1 );

  return data;
}

void NC::AlignedAlloc::detail::freeBigAlignedAlloc( void* addr )
{
  std::uintptr_t pn_addr = reinterpret_cast< std::uintptr_t >( addr );
  std::uintptr_t pn_first_aligned_address = pn_addr - sizeof(void*);
  void * tofree = *reinterpret_cast<void**>( pn_first_aligned_address );
  nc_std_free(tofree);
}


#ifdef NCRYSTAL_TRACKALIGNEDALLOC

namespace NCRYSTAL_NAMESPACE {
  namespace AlignedAlloc {
    namespace detail {
      namespace {
        static std::set<void*> trackalloc_active;
        static bool do_trackalignedalloc = false;
      }
      void trackalloc_enable( bool status )
      {
        if ( do_trackalignedalloc != status ) {
          trackalloc_active.clear();
          do_trackalignedalloc = status;
        }
      }
      void trackalloc_malloc(void* addr)
      {
        if (!do_trackalignedalloc)
          return;
        nc_assert_always(!trackalloc_active.count(addr));
        trackalloc_active.insert(addr);
      }

      void trackalloc_free(void* addr)
      {
        if (!do_trackalignedalloc)
          return;
        nc_assert_always(trackalloc_active.count(addr));
        trackalloc_active.erase(addr);
      }
      std::size_t trackalloc_count()
      {
        return trackalloc_active.size();
      }
    }
  }
}

#endif

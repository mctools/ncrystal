#ifndef NCrystal_SmallVector_hh
#define NCrystal_SmallVector_hh

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

#ifndef NCrystal_Defs_hh
#  include "NCrystal/core/NCDefs.hh"
#endif
#include <initializer_list>
#include <type_traits>

namespace NCRYSTAL_NAMESPACE {

  //This file provides two container classes, SmallVector and SmallVector_IC.
  //
  //SmallVector is a generic vector class similar in many ways to std::vector,
  //but using small buffer optimisation to keep contained objects on the
  //SmallVector itself when the number of such objects is <= NSMALL, thus
  //avoiding an allocation and a level of indirection at the expense of a larger
  //memory footprint. If more elements are added to the vector, it will fall
  //back to dynamically allocated heap storage. In any case, the elements are
  //kept in contiguous memory.
  //
  //For obvious reasons, the memory footprint of the class will be
  //NSMALL*sizeof(TValue) in addition to a small constant overhead which is
  //typically 8 bytes or 16 bytes depending amongst other things on the MODE and
  //the alignment requirements of TValue.
  //
  //If a SmallVector is determined to have fast access, a data member is added
  //which caches access to the beginning of the storage (local or heap), meaning
  //that access will be as fast as it is for std::vector (or better, due to the
  //better cache-locality when size()<=NSMALL). If MODE is LOWFOOTPRINT, fast
  //access is only enabled if the additional data member does not actually
  //increase the memory footprint (which it will unless TValue has very high
  //alignment requirements). If mode is FASTACCESS, fast access is always
  //enabled as the name implies. If fast access is not added, the consequence is
  //that each call to data(), begin(), operator[], etc. requires a test and
  //branch. This is probably not a big issues for most applications, but can be
  //worked around by only minimizing the number of such calls. E.g. instead of
  //"for (size_t i = 0; i<v.size();++i) v[i] =...;" instead do "auto
  //it=v.begin(); auto itE = it+v.size(); for (;it!=itE;++it) *it = ...;
  //
  //To additionally ensure good efficiency, contained TValue objects are
  //required to be noexcept-move-constructible and noexcept-destructible.
  //
  //The SmallVector_IC class is for all practical purposes identical to the
  //SmallVector class, with the only difference being that it is also implicitly
  //copyable (hence the suffix _IC).

  enum class SVMode : std::uint32_t { FASTACCESS=0, LOWFOOTPRINT=1 };

  namespace detail {
    class SV_CacheBegin;
    class SV_Empty;
    struct SV_Dummy {};
    template<class TValue, SVMode MODE>
    inline constexpr bool SVUseFast();
  }

  struct SVAllowCopy_t {};
  constexpr SVAllowCopy_t SVAllowCopy = SVAllowCopy_t{};
  struct SVCountConstruct_t {};
  constexpr SVCountConstruct_t SVCountConstruct = SVCountConstruct_t{};


  template<class TValue, std::size_t NSMALL, SVMode MODE = SVMode::FASTACCESS>
  class SmallVector : private MoveOnly,
                      protected std::conditional<detail::SVUseFast<TValue,MODE>(),
                                                 detail::SV_CacheBegin,
                                                 detail::SV_Empty>::type
  {
  public:
    typedef TValue element_type;
    typedef typename std::remove_cv< TValue >::type value_type;
    typedef TValue &       reference;
    typedef TValue *       pointer;
    typedef TValue const * const_pointer;
    typedef TValue const & const_reference;
    typedef pointer        iterator;
    typedef const_pointer  const_iterator;
    typedef std::ptrdiff_t difference_type;

    using size_type = decltype(NSMALL);
    static constexpr size_type nsmall = NSMALL;
    static constexpr SVMode mode = MODE;

    //TODO: Add custom iterators classes.

    ///////////////////////////////////////////////////////////////////////////
    //Default construction (as empty vector) or destruction are noexcept:
    constexpr SmallVector() noexcept;
    ~SmallVector() noexcept;

    ///////////////////////////////////////////////////////////////////////////
    //Move semantics and swap. This is a very cheap pointer swap when
    //size()>NSMALL (i.e. objects are on heap). If not, the contained objects
    //must be moved (but this will be at most NSMALL objects, and no memory
    //allocations will be required):
    SmallVector( SmallVector&& ) noexcept;
    SmallVector& operator=( SmallVector&& ) noexcept;
    void swap( SmallVector& ) noexcept;//also used by std::swap

    ///////////////////////////////////////////////////////////////////////////
    //No implict copies (if you need those, use the SmallVector_IC class
    //instead):
    SmallVector& operator=( const SmallVector& ) = delete;
    SmallVector( const SmallVector& ) = delete;

    ///////////////////////////////////////////////////////////////////////////
    //Comparison (first on size, then element-wise - this is different than
    //std::vector):
    bool operator==( const SmallVector& ) const noexcept;
    bool operator<( const SmallVector& ) const noexcept;

    ///////////////////////////////////////////////////////////////////////////
    //Braced initializer list initialisation (requires TValue to be copy
    //constructible):
    SmallVector( std::initializer_list<TValue> );

    ///////////////////////////////////////////////////////////////////////////
    //Can also move or copy (if TValue is copy constructible) from existing
    //containers. The SVAllowCopy trait is used to protect against accidental
    //copies:
    template <class TIter>
    SmallVector( SVAllowCopy_t, TIter it_begin, TIter it_end );
    template <class TIter>
    void setByCopy( TIter it_begin, TIter it_end );
    template <class TIter>
    void setByMove( TIter it_begin, TIter it_end );
    template<class TOther>
    SmallVector( SVAllowCopy_t, const TOther& o );

    //Construct from other container (the argument must be an r-value which will
    //be moved-from):
    template<class TOther, typename = typename
             std::enable_if<std::is_rvalue_reference<TOther&&>::value>::type >
    SmallVector( TOther&& o ) : SmallVector() { setByMove(o.begin(),o.end()); }

    ////////////////////////////////////////////////////////////////////////////
    //Access contents. Note that if isFastAccess() is false (which happens when
    //MODE is LOWFOOTPRINT* and TValue does not have unusually large alignment
    //requirements), each call has a small overhead since it needs to check if
    //heap or local storage is used. To ensure this overhead is only paid once
    //when looping over all contents, one can use code a'la:
    //
    // auto it = v.begin();
    // auto itE = it+v.size();
    // for ( ; it!=itE; ++it ) { ... }
    //
    //Or just create a Span from the vector and loop over that, which is
    //equivalent.

    ncnodiscard17 ncconstexpr17 TValue* data() noexcept;
    ncnodiscard17 ncconstexpr17 const TValue* data() const noexcept;
    ncnodiscard17 ncconstexpr17 TValue* begin() noexcept;
    ncnodiscard17 ncconstexpr17 TValue* end() noexcept;
    ncnodiscard17 constexpr const TValue* begin() const noexcept;
    ncnodiscard17 constexpr const TValue* end() const noexcept;
    ncnodiscard17 constexpr const TValue* cbegin() const noexcept { return begin(); }
    ncnodiscard17 constexpr const TValue* cend() const noexcept { return end(); }
    ncnodiscard17 ncconstexpr17 const value_type& operator[]( size_type ) const noexcept;
    ncnodiscard17 ncconstexpr17 value_type& operator[]( size_type ) noexcept;
    ncnodiscard17 ncconstexpr17 const value_type& at( size_type ) const;//checks and throws std::out_of_range
    ncnodiscard17 ncconstexpr17 value_type& at( size_type );//checks and throws std::out_of_range
    ncnodiscard17 ncconstexpr17 TValue& front() noexcept { return *begin(); }
    ncnodiscard17 constexpr const TValue& front() const noexcept { return *begin(); }
    ncnodiscard17 ncconstexpr17 TValue& back() noexcept { return *std::prev(end()); }
    ncnodiscard17 constexpr const TValue& back() const noexcept { return *std::prev(end()); }

    ///////////////////////////////////////////////////////////////////////////
    //Various standard methods:
    ncnodiscard17 constexpr size_type size() const noexcept;
    ncnodiscard17 constexpr bool empty() const noexcept;
    ncnodiscard17 ncconstexpr17 size_type capacity() const noexcept;

    ///////////////////////////////////////////////////////////////////////////
    //Insert or remove elements:
    void push_back( const TValue& value );
    void push_back( TValue&& value );
    template<typename ...Args>
    TValue& emplace_back( Args&& ... );
    void pop_back() noexcept;
    void clear() noexcept;

    //Calls to resize(n) requires TValue to have noexcept default construction
    //or noexcept copy construction, while calls to resize(n,val) needs noexcept
    //copy construction:
    template<class U = TValue, typename = typename
             std::enable_if<std::is_nothrow_default_constructible<U>::value>::type >
    void resize( size_type );

    template<class U = TValue, typename = typename
             std::enable_if<(!std::is_nothrow_default_constructible<U>::value
                             &&std::is_nothrow_copy_constructible<U>::value)>::type >
    void resize( size_type n, detail::SV_Dummy * = nullptr );

    template<class U = TValue, typename = typename
             std::enable_if<std::is_nothrow_copy_constructible<U>::value>::type >
    void resize( size_type, const TValue& value );

    ///////////////////////////////////////////////////////////////////////////
    //Capacity modification methods. By design, the capacity is always NSMALL
    //when size()<=NSMALL, so the following methods only have any effect when
    //size()>NSMALL. As above, state is unchanged in case of memory allocation
    //error:
    void shrink_to_fit();
    void reserve_hint( size_type n );

    //SmallVector specifics:
    constexpr bool isLocalStorage() const noexcept;
    static constexpr bool isFastAccess() noexcept;

    ///////////////////////////////////////////////////////////////////////////
    //SmallVector does not have constructors like std::vector(count) or
    //std::vector(count,value). Deleting these explicitly gives better error
    //diagnostics for the user (specifying various bit widths for robustness):
    SmallVector( uint8_t ) = delete;
    SmallVector( uint16_t ) = delete;
    SmallVector( uint32_t ) = delete;
    SmallVector( uint64_t ) = delete;
    SmallVector( int8_t ) = delete;
    SmallVector( int16_t ) = delete;
    SmallVector( int32_t ) = delete;
    SmallVector( int64_t ) = delete;

    ///////////////////////////////////////////////////////////////////////////
    //The std::vector(size_type) and std::vector(size_type,const TValue&)
    //constructors are known to cause confusion for integer-convertible types
    //(i.e. does vector{ 10 } give 10 default constructed elements or a single
    //element equal 10?). Thus for SmallVector we protect against accidental
    //invocation by requiring the SVCountConstruct trait. To be used like:
    //  SmallVector<MyType,4> myvect( SVCountConstruct, 20 );
    //  SmallVector<MyType,4> myvect( SVCountConstruct, 20, MyType("bla") )
    //SmallVector( SVCountConstruct_t, size_type count );
    SmallVector( SVCountConstruct_t, size_type count, const TValue& );

    //To make std::swap work (also has additional global functions inlined
    //below):
    friend void swap( SmallVector& a, SmallVector& b ) noexcept { a.swap(b); }

  private:
    size_type m_count = 0;
    union {
      struct {
        TValue * data;
        size_type capacity;
      } large;
      //NB: We do not use std::array for local storage, since we rely on address
      //of m_data being the same as the beginning of the local storage, which is
      //certainly guaranteed when we use a C-style array:
      alignas(TValue) uint8_t small_data[ NSMALL * sizeof(TValue) ];
    } m_data;
    struct Impl;
  };

  /////////////////////////////////////////////////////////////////////////////
  //We implement the SmallVector_IC class by inheritance from the SmallVector
  //class, and then implementing copy/move as required:

  template<class TValue, std::size_t NSMALL, SVMode MODE = SVMode::FASTACCESS>
  class SmallVector_IC final : public SmallVector<TValue,NSMALL,MODE>
  {
  public:
    using SmallVector<TValue,NSMALL,MODE>::SmallVector;
    constexpr SmallVector_IC() noexcept {}
    ~SmallVector_IC() noexcept {}
    SmallVector_IC( SmallVector_IC&& ) noexcept;
    SmallVector_IC& operator=( SmallVector_IC&& ) noexcept;
    SmallVector_IC( const SmallVector_IC& );
    SmallVector_IC& operator=( const SmallVector_IC& );
    //It can be access as a reference to the corresponding SmallVector:
    using base_t = SmallVector<TValue,NSMALL,MODE>;
    base_t& asBase() noexcept;
    const base_t& asBase() const noexcept;
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace detail {
    using SV_size_type = std::size_t;
    class NCRYSTAL_API SV_CacheBegin {
      void * m_begin;
    protected:
      SV_CacheBegin( void * b ) : m_begin{b} {}
      ncconstexpr17 void setBeginPtr( void* p ) noexcept { m_begin = p; }
      ncconstexpr17 void* beginPtr() const noexcept { return m_begin; }
    };
    class NCRYSTAL_API SV_Empty {
    protected:
      SV_Empty( void * ) {}
      ncconstexpr17 void setBeginPtr( void* ) noexcept {}
      ncconstexpr17 void* beginPtr() const noexcept { return nullptr; }
    };

    template<class TValue>
    class SVTestCB : SV_CacheBegin { SV_size_type a; struct alignas(TValue) { char b; } c; };
    template<class TValue>
    class SVTestE  : SV_Empty { SV_size_type a; struct alignas(TValue) { char b; } c; };

    template<class TValue, SVMode MODE>
    inline constexpr bool SVUseFast() {
      return ( MODE==SVMode::FASTACCESS || sizeof(SVTestE<TValue>)==sizeof(SVTestCB<TValue>) );
    }
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  struct SmallVector<TValue,NSMALL,MODE>::Impl {

    using SVBase = typename std::conditional<detail::SVUseFast<TValue,MODE>(),
                                             detail::SV_CacheBegin,
                                             detail::SV_Empty>::type;

    static_assert( std::is_same<detail::SV_size_type,size_type>::value, "");
    static_assert(NSMALL>0,"SmallVector can not have NSMALL=0.");
    static_assert( (NSMALL*sizeof(TValue)) < (size_type)100000000ull,
                  "SmallVector local storage configured to use more than 100MB.");
    static_assert(std::is_nothrow_destructible<TValue>::value,
                  "SmallVector can only keep objects with noexcept destructors");
    static_assert(std::is_move_constructible<TValue>::value,
                  "SmallVector can only keeps objects with move constructors.");
    static_assert(std::is_nothrow_move_constructible<TValue>::value,
                  "SmallVector can only keep objects with noexcept move-constructors.");

    static constexpr bool large(const SmallVector* THIS) noexcept { return THIS->m_count > NSMALL; }
    static constexpr bool small(const SmallVector* THIS) noexcept { return THIS->m_count <= NSMALL; }

    class DetachedHeap {
      TValue * m_begin;
      TValue * m_end;
      size_type m_capacity;
    public:
      DetachedHeap(TValue*b,TValue*e,size_type cc) : m_begin(b), m_end(e), m_capacity(cc) {}
      ncconstexpr17 TValue * begin() noexcept
      {
        assert( m_begin != nullptr );
        return m_begin;
      }
      ncconstexpr17 TValue * end() noexcept
      {
        assert( m_end != nullptr );
        return m_end;
      }
      constexpr size_type capacity() const noexcept
      {
#if nc_cplusplus >= 201402L // assert in constexpr requires C++14
        assert( m_begin != nullptr);
#endif
        return m_capacity;
      }
      TValue * release() { TValue * d = m_begin; m_begin=m_end=nullptr; return d; }
      template<typename ...Args>
      void emplace_back( Args&& ...args )
      {
        assert( m_begin != nullptr );
        //NB: calling code is responsible for ensuring adequate capacity.
        new (m_end) TValue(std::forward<Args>(args)...);
        ++m_end;//on line after TValue constructor (in case it throws)
        assert( m_end <= m_begin + m_capacity );
      }

      ~DetachedHeap()
      {
        if ( m_begin ) {
          auto it(m_begin), itE(m_end);
          for ( ; it!=itE; ++it )
            it->~TValue();
          AlignedAlloc::freeAlignedAlloc<TValue>( m_begin );
        }
      }
    };

    static void adoptHeap(SmallVector* THIS,DetachedHeap& heap) noexcept
    {
      size_type count = std::distance( heap.begin(), heap.end() );
      THIS->clear();
      assert(count>NSMALL);
      THIS->m_data.large.capacity = heap.capacity();
      THIS->setBeginPtr( THIS->m_data.large.data = heap.release() );
      THIS->m_count = count;
    }

    static DetachedHeap createNewDetachedHeap( size_type capacity )
    {
      //Can throw std::bad_alloc
      TValue * b = AlignedAlloc::alignedAlloc<TValue>(capacity);
      return DetachedHeap(b,b,capacity);
    }

    static DetachedHeap detachHeapDataAndClear(SmallVector * THIS) noexcept
    {
      assert(large(THIS));
      TValue * b = THIS->m_data.large.data;
      size_type capacity = THIS->m_data.large.capacity;
      size_type count = THIS->m_count;
      THIS->m_data.large.data = nullptr;
      THIS->m_count = 0;
      setBeginPtrSmallData(THIS);
      return DetachedHeap(b,b+count,capacity);
    }

    static void resizeLargeCapacity( SmallVector* THIS, size_type n )
    {
      //leaves THIS unchanged if createNewDetachedHeap throws bad_alloc
      assert( large(THIS) );
      assert( n >= THIS->m_count );
      auto heap = createNewDetachedHeap(n);
      for ( auto&& e : *THIS )
        heap.emplace_back(std::move(e));
      adoptHeap(THIS,heap);
    }

    static constexpr void * smallDataBegin(SmallVector * THIS) noexcept
    {
      //We rely on the fact that the address of m_data.small_data is guaranteed
      //to be the same as the address of m_data, irrespective of which data
      //member is active.
      return static_cast<void*>(&(THIS->m_data));
    }
    static ncconstexpr17 void setBeginPtrSmallData(SmallVector * THIS) noexcept
    {
      THIS->setBeginPtr(smallDataBegin(THIS));
    }

    static void clear(SmallVector * THIS ) noexcept
    {
      //Call destructors, release heap alloction (if any) and set count to
      //0. It is noexcept since destructors should not throw.
      if ( THIS->m_count == 0 )
        return;

      if ( large(THIS) ) ncunlikely {
        //Since v3.7.0 (summer 2023) we release the large area like this, to
        //avoid what looks like a spurious error with gcc 12 (see
        //https://github.com/mctools/ncrystal/issues/125). Doing it like this
        //seems to avoid this issue for some obscure reason, perhaps because it
        //lets gcc focus on one RAII class at a time. Also, it is important
        //(apparently!) that the "THIS->m_count == 0" check above comes first.
        Impl::detachHeapDataAndClear(THIS);
        return;
      }

      auto itE = THIS->end();
      for ( auto it = THIS->begin(); it!=itE ; ++it )
        it->~TValue();
      THIS->m_count = 0;
      setBeginPtrSmallData(THIS);
    }

    static void resizeDown( SmallVector * THIS, size_type n ) noexcept
    {
      assert ( THIS->m_count >= n );
      if ( THIS->m_count == n )
        return;
      if  ( n > NSMALL || THIS->m_count <= NSMALL ) {
        //Can peel entries off without changing storage mode:
        auto it = THIS->begin() + n;
        auto itE = THIS->end();
        assert( itE > it );
        for ( ; it != itE; ++it )
          it->~TValue();
        THIS->m_count = n;
      } else {
        assert( n <= NSMALL );
        resizeDown( THIS, NSMALL+1 );//peel off heap storage entries
        THIS->pop_back();//peel of one entry which switches storages from heap to local
        resizeDown( THIS, n );//peel off local storage entries
      }
    }


    template<typename ...Args>
    static TValue& grow_and_emplace_back( SmallVector * THIS, Args&& ...args )
    {
      assert( THIS->m_count == THIS->capacity() );
      //In this edge-case we actually construct the object here on the stack
      //first, and then afterwards increase capacity and move it over. This is
      //to make sure the SmallVector state is unchanged in case the constructor
      //throws:
      TValue newvalue(std::forward<Args>(args)...);
      if ( THIS->m_count == NSMALL ) {
        //Special case, moving from small to large buffer!
        auto heap = createNewDetachedHeap( NSMALL*2 );//might throw bad_alloc
        //Ok, done with everything that might throw, it is now safe to start
        //modifying our state:
        for ( auto&& e : *THIS )
          heap.emplace_back(std::move(e));//noexcept move constructor
        heap.emplace_back(std::move(newvalue));
        TValue * last = std::prev(heap.end());
        adoptHeap( THIS, heap );
        return *last;
      } else {
        resizeLargeCapacity( THIS, THIS->m_count*2 );
        return THIS->emplace_back(std::move(newvalue));
      }
    }

    template<typename ...Args>
    static TValue& emplace_back( SmallVector * THIS, Args&& ...args )
    {
      if ( THIS->m_count < THIS->capacity() ) nclikely {
        //If TValue constructor throws we must not change state. This is easily
        //handled by only incrementing m_count on the line AFTER trying the
        //construction (a partially constructed object should not have its
        //destructor run so there won't be any cleanup):
        TValue * newobjaddr = THIS->end();
        new( (void*)(newobjaddr) ) TValue(std::forward<Args>(args)...);
        ++(THIS->m_count);
        return *newobjaddr;
      } else {
        return grow_and_emplace_back(THIS, std::forward<Args>(args)...);
      }
    }
  };

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline SmallVector<TValue,NSMALL,MODE>::SmallVector( SmallVector&& o ) noexcept
    : SmallVector()
  {
    *this = std::move(o);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline SmallVector<TValue,NSMALL,MODE>& SmallVector<TValue,NSMALL,MODE>::operator=( SmallVector&& o ) noexcept
  {
    if ( this == &o )
      return *this;
    if ( m_count )
      clear();
    if ( Impl::small(&o) ) {
      //Move values:
      auto it = begin();
      for ( auto&& e : o )
        new(it++) TValue(std::move(e));
      m_count = o.m_count;
      o.clear();
      Impl::setBeginPtrSmallData(this);
    } else {
      //Simply update pointers!
      std::swap(m_count,o.m_count);
      this->setBeginPtr( m_data.large.data = o.m_data.large.data);
      m_data.large.capacity = o.m_data.large.capacity;
      o.m_data.large.capacity = 0;//why did we leave out this line earlier?
      Impl::setBeginPtrSmallData(&o);
    }
    assert(o.empty());
    assert(!Impl::large(&o));
    return *this;
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline SmallVector<TValue,NSMALL,MODE>::SmallVector( std::initializer_list<TValue> l)
    : SmallVector()
  {
    if ( l.size() <= NSMALL ) {
      auto it = begin();
      for ( auto&& e : l )
        new(it++) TValue(std::move(e));
      m_count = l.size();
    } else {
      auto heap = Impl::createNewDetachedHeap( l.size() );
      for ( auto&& e : l )
        heap.emplace_back(std::move(e));
      Impl::adoptHeap(this,heap);
    }
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  template <class TIter>
  inline void SmallVector<TValue,NSMALL,MODE>::setByCopy( TIter it_begin, TIter it_end ) {
    clear();
    std::size_t n = std::distance(it_begin,it_end);
    if ( n <= NSMALL ) {
      auto it = begin();
      for ( ; it_begin!=it_end; ++it_begin)
        new(it++) TValue(*it_begin);
      m_count = n;
    } else {
      auto heap = Impl::createNewDetachedHeap( n );
      for ( ; it_begin!=it_end; ++it_begin)
        heap.emplace_back( *it_begin );
      Impl::adoptHeap(this,heap);
    }
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  template <class TIter>
  inline void SmallVector<TValue,NSMALL,MODE>::setByMove( TIter it_begin, TIter it_end ) {
    clear();
    std::size_t n = std::distance(it_begin,it_end);
    if ( n <= NSMALL ) {
      auto it = begin();
      for ( ; it_begin!=it_end; ++it_begin)
        new(it++) TValue(std::move(*it_begin));
      m_count = n;
    } else {
      auto heap = Impl::createNewDetachedHeap( n );
      for ( ; it_begin!=it_end; ++it_begin)
        heap.emplace_back( std::move(*it_begin) );
      Impl::adoptHeap(this,heap);
    }
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline constexpr SmallVector<TValue,NSMALL,MODE>::SmallVector() noexcept
    : Impl::SVBase(Impl::smallDataBegin(this))
  {
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  template<class TOther>
  inline SmallVector<TValue,NSMALL,MODE>::SmallVector( SVAllowCopy_t, const TOther& o )
    : SmallVector(SVAllowCopy,o.begin(),o.end())
  {
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  template <class TIter>
  SmallVector<TValue,NSMALL,MODE>::SmallVector( SVAllowCopy_t, TIter it_begin, TIter it_end )
    : SmallVector()
  {
    this->setByCopy(it_begin,it_end);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  ncnodiscard17 inline ncconstexpr17 TValue* SmallVector<TValue,NSMALL,MODE>::data() noexcept
  {
    if (detail::SVUseFast<TValue,MODE>())
      return static_cast<TValue*>(this->beginPtr());
    if ( Impl::small(this) ) nclikely
      return reinterpret_cast<TValue*>(&m_data.small_data[0]);
    else
      return m_data.large.data;
  }
  template<class TValue, std::size_t NSMALL, SVMode MODE>
  ncnodiscard17 inline ncconstexpr17 const TValue* SmallVector<TValue,NSMALL,MODE>::data() const noexcept {
    if (detail::SVUseFast<TValue,MODE>())
      return static_cast<const TValue*>(this->beginPtr());
    if ( Impl::small(this) ) nclikely
      return reinterpret_cast<const TValue*>(&m_data.small_data[0]);
    else
      return m_data.large.data;
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  ncnodiscard17 inline constexpr typename SmallVector<TValue,NSMALL,MODE>::size_type
  SmallVector<TValue,NSMALL,MODE>::size() const noexcept
  {
    return m_count;
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  ncnodiscard17 inline constexpr bool SmallVector<TValue,NSMALL,MODE>::empty() const noexcept
  {
    return m_count == 0;
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  ncnodiscard17 inline ncconstexpr17 const typename SmallVector<TValue,NSMALL,MODE>::value_type&
  SmallVector<TValue,NSMALL,MODE>::operator[]( size_type i ) const noexcept
  {
    assert(i<m_count);
    return *std::next(data(),i);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  ncnodiscard17 inline ncconstexpr17 typename SmallVector<TValue,NSMALL,MODE>::value_type&
  SmallVector<TValue,NSMALL,MODE>::operator[]( size_type i ) noexcept
  {
    assert(i<m_count);
    return *std::next(data(),i);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  ncnodiscard17 inline ncconstexpr17 const typename SmallVector<TValue,NSMALL,MODE>::value_type& SmallVector<TValue,NSMALL,MODE>::at( size_type i ) const
  {
    if ( i >= m_count )
      throw std::out_of_range("SmallVector::at(): index out of out of range");
    return *std::next(data(),i);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  ncnodiscard17 inline ncconstexpr17 typename SmallVector<TValue,NSMALL,MODE>::value_type& SmallVector<TValue,NSMALL,MODE>::at( size_type i )
  {
    if ( i >= m_count )
      throw std::out_of_range("SmallVector::at(): index out of out of range");
    return *std::next(data(),i);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  ncnodiscard17 inline ncconstexpr17 TValue* SmallVector<TValue,NSMALL,MODE>::begin() noexcept { return data(); }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  ncnodiscard17 inline ncconstexpr17 TValue* SmallVector<TValue,NSMALL,MODE>::end() noexcept { return data() + size(); }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  ncnodiscard17 inline constexpr const TValue* SmallVector<TValue,NSMALL,MODE>::begin() const noexcept { return data(); }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  ncnodiscard17 inline constexpr const TValue* SmallVector<TValue,NSMALL,MODE>::end() const noexcept { return data() + size(); }


  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline void SmallVector<TValue,NSMALL,MODE>::push_back( const TValue& value )
  {
    Impl::emplace_back(this,value);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline void SmallVector<TValue,NSMALL,MODE>::push_back( TValue&& value )
  {
    Impl::emplace_back(this,std::move(value));
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  template<typename ...Args>
  inline TValue& SmallVector<TValue,NSMALL,MODE>::emplace_back( Args&& ...args )
  {
    return Impl::emplace_back(this,std::forward<Args>(args)...);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline void SmallVector<TValue,NSMALL,MODE>::clear() noexcept
  {
    Impl::clear(this);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline constexpr bool SmallVector<TValue,NSMALL,MODE>::isLocalStorage() const noexcept
  {
    return Impl::small(this);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline constexpr bool SmallVector<TValue,NSMALL,MODE>::isFastAccess() noexcept
  {
    return detail::SVUseFast<TValue,MODE>();
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  ncnodiscard17 inline ncconstexpr17 typename SmallVector<TValue,NSMALL,MODE>::size_type
  SmallVector<TValue,NSMALL,MODE>::capacity() const noexcept
  {
    if ( Impl::small(this) ) nclikely
      return NSMALL ;
    else
      return m_data.large.capacity;
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline void SmallVector<TValue,NSMALL,MODE>::reserve_hint( size_type n ) {
    if ( Impl::small(this) || n <= capacity() )
      return;
    Impl::resizeLargeCapacity(this,n);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline void SmallVector<TValue,NSMALL,MODE>::shrink_to_fit() {
    if ( Impl::small(this) || m_count == m_data.large.capacity )
      return;
    Impl::resizeLargeCapacity( this, m_count );
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline SmallVector<TValue,NSMALL,MODE>::~SmallVector() noexcept
  {
    Impl::clear(this);
  }


  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline void SmallVector<TValue,NSMALL,MODE>::swap(SmallVector& o) noexcept
  {
    if ( Impl::large(this) && Impl::large(&o) ) {
      //Both use heap-storage, extremely fast:
      std::swap( m_count, o.m_count );
      std::swap( m_data.large.data, o.m_data.large.data );
      std::swap( m_data.large.capacity, o.m_data.large.capacity );
      o.setBeginPtr( o.m_data.large.data );
      this->setBeginPtr( m_data.large.data );
    } else {
      //Potentially complicated, just use move semantics:
      SmallVector<TValue,NSMALL,MODE> tmp(std::move(*this));
      *this = std::move(o);
      o = std::move(tmp);
    }
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  template<class U, typename>
  inline void SmallVector<TValue,NSMALL,MODE>::resize( size_type n,
                                                       detail::SV_Dummy * )
  {
    n == m_count+1 ? this->emplace_back() : this->resize( n, TValue() );
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  template<class U, typename>
  inline void SmallVector<TValue,NSMALL,MODE>::resize( size_type n )
  {
    static_assert( std::is_nothrow_default_constructible<TValue>::value, "" );
    if ( m_count >= n ) {
      Impl::resizeDown( this, n );
      return;
    }

    //Must add new default-constructed elements at the end.
    if ( n <= capacity() ) {
      //In existing capacity (since TValue is_nothrow_default_constructible,
      //we can do this directly without intermediate buffer):
      while ( m_count < n ) {
        new (end()) TValue();
        ++m_count;
      }
      assert( (size_type)std::distance(begin(),end()) == m_count );
      return;
    }

    //Needs more capacity - use new heap-storage into which we move over old
    //elements and append new default constructed elements as needed (again
    //relying on TValue to be is_nothrow_default_constructible):
    assert( n > NSMALL );
    auto heap = Impl::createNewDetachedHeap( n );
    for ( auto&&e : *this )
      heap.emplace_back(std::move(e));
    for ( size_type i = m_count; i < n; ++i )
      heap.emplace_back();//TValue() is noexcept
    assert( (size_type)std::distance(heap.begin(),heap.end()) == n );
    Impl::adoptHeap(this,heap);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  template<class U, typename>
  inline void SmallVector<TValue,NSMALL,MODE>::resize( size_type n,
                                                       const TValue& val_to_copy )
  {
    static_assert(std::is_nothrow_copy_constructible<TValue>::value,"");
    if ( m_count >= n ) {
      Impl::resizeDown( this, n );
      return;
    }

    //Must add new copy-constructed elements at the end.
    if ( n <= capacity() ) {
      //In existing capacity (since TValue is_nothrow_copy_constructible,
      //we can do this directly without intermediate buffer):
      while ( m_count < n ) {
        new (end()) TValue(val_to_copy);
        ++m_count;
      }
      assert( (size_type)std::distance(begin(),end()) == m_count );
      return;
    }

    //Needs more capacity - use new heap-storage into which we move over old
    //elements and append new copy constructed elements as needed (again
    //relying on TValue to be is_nothrow_copy_constructible):
    assert( n > NSMALL );
    auto heap = Impl::createNewDetachedHeap( n );
    for ( auto&&e : *this )
      heap.emplace_back(std::move(e));
    for ( size_type i = m_count; i < n; ++i )
      heap.emplace_back(val_to_copy);//TValue(val_to_copy) is noexcept
    assert( (size_type)std::distance(heap.begin(),heap.end()) == n );
    Impl::adoptHeap(this,heap);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline void SmallVector<TValue,NSMALL,MODE>::pop_back() noexcept
  {
    if ( m_count && m_count != (NSMALL+1) ) nclikely {
      //Simply destruct last entry and move count down.
      ( data() + --m_count)->~TValue();
      return;
    }
    if ( !m_count ) ncunlikely {
      //std::vector::pop_back on empty vector is UB. Here we simply choose the
      //easy (and well defined ) solution of asserting and otherwise doing
      //nothing.
      assert( false && "pop_back on empty vector is undefined behaviour");
      return;
    }

    ///////////////////////////////////////////////////////////////////
    //Pop_back causes transition from large to small.

    //First detach heap allocation clear our own state:
    auto old_heap = Impl::detachHeapDataAndClear(this);

    //Then regain by move the NSMALL values which are to be kept:
    assert( std::next(old_heap.begin(),NSMALL+1) == old_heap.end() );
    this->setByMove( old_heap.begin(), std::next(old_heap.begin(),NSMALL) );
    assert(m_count==NSMALL);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline SmallVector<TValue,NSMALL,MODE>::SmallVector( SVCountConstruct_t, size_type count, const TValue& value )
    : SmallVector()
  {
    if ( count <= NSMALL ) {
      if ( count == 0 )
        return;
      while (count--)
        this->push_back(value);
      return;
    } else {
      auto heap = Impl::createNewDetachedHeap( count );
      while (count--)
        heap.emplace_back(value);
      Impl::adoptHeap(this,heap);
    }
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline void swap(SmallVector<TValue,NSMALL,MODE>& a, SmallVector<TValue,NSMALL,MODE>& b) noexcept
  {
    a.swap(b);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline bool SmallVector<TValue,NSMALL,MODE>::operator==( const SmallVector& o ) const noexcept
  {
    if ( this->m_count != o.m_count )
      return false;
    if ( this == &o || this->m_count == 0 )
      return true;
    auto it = this->begin();
    auto itE = it + this->m_count;
    auto itO = o.begin();
    for ( ; it != itE; ++it, ++itO )
      if ( ! ( *it == *itO ) )
        return false;
    return true;
  }


  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline bool SmallVector<TValue,NSMALL,MODE>::operator<( const SmallVector& o ) const noexcept
  {
    if ( this->m_count != o.m_count )
      return this->m_count < o.m_count;
    if ( this == &o || this->m_count == 0 )
      return false;
    auto it = this->begin();
    auto itE = it + this->m_count;
    auto itO = o.begin();
    for ( ; it != itE; ++it, ++itO )
      if ( ! ( *it == *itO ) )
        return *it < *itO;
    return false;
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  SmallVector<TValue,NSMALL,MODE>& SmallVector_IC<TValue,NSMALL,MODE>::asBase() noexcept
  {
    return *static_cast<SmallVector<TValue,NSMALL,MODE>*>(this);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  const SmallVector<TValue,NSMALL,MODE>& SmallVector_IC<TValue,NSMALL,MODE>::asBase() const noexcept
  {
    return *static_cast<SmallVector<TValue,NSMALL,MODE>*>(this);
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline SmallVector_IC<TValue,NSMALL,MODE>::SmallVector_IC( SmallVector_IC&& o ) noexcept
    : SmallVector_IC()
  {
    this->asBase() = std::move(o.asBase());
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline SmallVector_IC<TValue,NSMALL,MODE>& SmallVector_IC<TValue,NSMALL,MODE>::operator=( SmallVector_IC&& o ) noexcept
  {
    this->asBase() = std::move(o.asBase());
    return *this;
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline SmallVector_IC<TValue,NSMALL,MODE>& SmallVector_IC<TValue,NSMALL,MODE>::operator=( const SmallVector_IC& o )
  {
    *this = SmallVector_IC(SVAllowCopy,o);
    return *this;
  }

  template<class TValue, std::size_t NSMALL, SVMode MODE>
  inline SmallVector_IC<TValue,NSMALL,MODE>::SmallVector_IC( const SmallVector_IC& o )
    : SmallVector_IC(SVAllowCopy,o)
  {
  }

  static_assert(std::is_copy_constructible<SmallVector_IC<int,1>>::value,"");
  static_assert(std::is_copy_assignable<SmallVector_IC<int,1>>::value,"");
  static_assert(!std::is_copy_constructible<SmallVector<int,1>>::value,"");
  static_assert(!std::is_copy_assignable<SmallVector<int,1>>::value,"");
  static_assert(std::is_nothrow_default_constructible<SmallVector<int,1>>::value,"");
  static_assert(std::is_nothrow_default_constructible<SmallVector_IC<int,1>>::value,"");
}

namespace std {
  //Would like to avoid injections in std::, but it seems to not work reliably without it:
  template<class TValue, std::size_t NSMALL, NCrystal::SVMode MODE>
  inline void swap(NCrystal::SmallVector<TValue,NSMALL,MODE>& a, NCrystal::SmallVector<TValue,NSMALL,MODE>& b)
  {
    a.swap(b);
  }
}

#endif

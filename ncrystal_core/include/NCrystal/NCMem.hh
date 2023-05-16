#ifndef NCrystal_Mem_hh
#define NCrystal_Mem_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2023 NCrystal developers                                   //
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

#include "NCrystal/NCException.hh"
#include <memory>
#include <type_traits>
#include <utility>//std::move
#include <cassert>
#include <functional>
#include <atomic>
#include <cstddef>
#include <new>

namespace NCrystal {

  //Attempt to clear all NCrystal caches (should be safe to call as it will not
  //clear data associated to active object for which client code has ownership):
  NCRYSTAL_API void clearCaches();

  //For internal NCrystal usage, registered functions will be invoked whenever
  //clearCaches() is called:
  NCRYSTAL_API void registerCacheCleanupFunction(std::function<void()>);

  //Type alias for std::shared_ptr which makes it clear when to use shared_obj
  //and when to use the nullable alternative:
  template <class T>
  using optional_shared_obj = std::shared_ptr<T>;

  template <class T>
  class shared_obj {

    // This thin wrapper around std::shared_ptr can never be initialised with a
    // null-ptr (constructors will throw if this is attempted), and it can be
    // passed directly to functions taking a reference or pointer to the
    // contained object type. It can also easily be converted to/from
    // std::shared_ptr as needed.
    //
    // Note that a moved-from shared_obj will still contain a nullptr (however,
    // one should never use a moved-from object for anything else than being
    // destructed).

    //For simplicity, explicitly disallow esoteric types:
    static_assert(!std::is_volatile<T>::value,"shared_obj class does not support volative types");
    static_assert(!std::is_reference<T>::value,"shared_obj class does not support reference types");
    static_assert(!std::is_member_pointer<T>::value,"shared_obj class does not support pointer types");
    static_assert(!std::is_member_object_pointer<T>::value,"shared_obj class does not support pointer types");
    static_assert(!std::is_member_function_pointer<T>::value,"shared_obj class does not support pointer types");
    static_assert(!std::is_pointer<T>::value,"shared_obj class does not support pointer types");
    static_assert(!std::is_array<T>::value,"shared_obj class does not support array types");

  public:

    using element_type = T;
    using element_type_nocv = typename std::remove_cv<T>::type;
    using shared_ptr_type = std::shared_ptr<T>;

    constexpr shared_obj() = delete;//Can not be default constructed
    explicit ncconstexpr17 shared_obj( T* t ) : m_shptr(t) { check_nonnull(); }//unsafe if not explict
    ncconstexpr17 shared_obj( std::shared_ptr<T> sp ) : m_shptr(std::move(sp)) { check_nonnull(); }
    ncconstexpr17 shared_obj( std::unique_ptr<T> up ) : m_shptr(std::move(up)) { check_nonnull(); }

    //If T is const, we can construct from non-const objects:
    template<typename U = T>
    ncconstexpr17 shared_obj( std::shared_ptr<element_type_nocv> sp, typename std::enable_if<std::is_const<U>::value>::type* = nullptr ) : m_shptr(std::move(sp)) { check_nonnull(); }
    template<typename U = T>
    ncconstexpr17 shared_obj( std::unique_ptr<element_type_nocv> up, typename std::enable_if<std::is_const<U>::value>::type* = nullptr ) : m_shptr(std::move(up)) { check_nonnull(); }
    template<typename U = T>
    ncconstexpr17 shared_obj( element_type_nocv* t, typename std::enable_if<std::is_const<U>::value>::type* = nullptr ) : m_shptr(t) { check_nonnull(); }
    template<typename U = T>
    constexpr shared_obj( const shared_obj<element_type_nocv>& so, typename std::enable_if<std::is_const<U>::value>::type* = nullptr ) noexcept : m_shptr(so) { }

    ncconstexpr17 shared_obj( const shared_obj& ) = default;
    ncconstexpr17 shared_obj& operator=( const shared_obj& ) = default;
    ncconstexpr17 shared_obj( shared_obj&& ) = default;
    ncconstexpr17 shared_obj& operator=( shared_obj&& ) = default;

    //Disallow initialisation/assignment from nullptr:
    shared_obj & operator= ( decltype(nullptr) ) = delete;
    shared_obj( decltype(nullptr) ) = delete;

    //Explicit getters:
    ncconstexpr17 T* get() noexcept { return m_shptr.get(); }
    constexpr const T* get() const noexcept { return m_shptr.get(); }
    constexpr const std::shared_ptr<T>& getsp() const noexcept { return m_shptr; }
    constexpr const optional_shared_obj<T>& optional() const noexcept { return m_shptr; }//same, different names

    //Object is always true in boolean contexts and not equal to nullptr:
    constexpr bool operator()() const noexcept { return true; }
    constexpr bool operator!() const noexcept { return false; }
    constexpr bool operator==( decltype(nullptr) ) const noexcept { return false; }
    constexpr bool operator!=( decltype(nullptr) ) const noexcept { return true; }

    //comparisons:
    constexpr bool operator==( const shared_obj& o ) const noexcept { return m_shptr == o.m_shptr; }
    constexpr bool operator!=( const shared_obj& o ) const noexcept { return m_shptr != o.m_shptr; }
    constexpr bool operator<( const shared_obj& o ) const noexcept { return m_shptr < o.m_shptr; }

    //deref:
    ncconstexpr17 T* operator->() noexcept { return m_shptr.get(); }
    constexpr const T* operator->() const noexcept { return m_shptr.get(); }
    ncconstexpr17 T& operator*() noexcept { return *m_shptr; }
    constexpr const T& operator*() const noexcept { return *m_shptr; }

    //Automatic conversions:
    constexpr operator const T&() const noexcept { return *m_shptr; }
    ncconstexpr17 operator T&() noexcept { return *m_shptr; }
    constexpr operator std::shared_ptr<const T>() const noexcept { return m_shptr; }
    ncconstexpr17 operator std::shared_ptr<T>() noexcept { return m_shptr; }
    constexpr operator std::weak_ptr<const T>() const noexcept { return m_shptr; }
    ncconstexpr17 operator std::weak_ptr<T>() noexcept { return m_shptr; }

    //Conversion to pointer is explicit to avoid footgun "const Info* info =
    //NC::createInfo(..)". The footgun is still available for references (const
    //Info& info = ...), but that is a necessary tradeoff since the implicit
    //conversion to a reference is so useful when calling functions:
    explicit constexpr operator const T*() const noexcept { return m_shptr.get(); }
    explicit ncconstexpr17 operator T*() noexcept { return m_shptr.get(); }

    //If T is non-const, we can convert to const object:
    template<class U = T>
    constexpr operator typename std::enable_if<std::negate<std::is_const<U>>::value, shared_obj<const T>>::type() const noexcept
    {
      return shared_obj<const T>(m_shptr);
    }

    //Move/copy/assign to base class pointer (in shared_obj or std::shared_ptr)
    //should be ok (only compiles if T* can be assigned from TOther*):
    template<class TOther>
    ncconstexpr17 shared_obj& operator=( shared_obj<TOther>&& o ) noexcept {
      m_shptr = std::move(o).detail_consume_shptr();
      check_nonnull();
      return *this;
    }

    template<class TOther>
    ncconstexpr17 shared_obj( shared_obj<TOther>&& o ) noexcept : m_shptr(std::move(o).detail_consume_shptr()) {}

    template<class TOther>
    ncconstexpr17 shared_obj& operator=( const shared_obj<TOther>& o ) noexcept
    {
      m_shptr = o.m_shptr;
      return *this;
    }

    template<class TOther>
    ncconstexpr17 shared_obj( std::shared_ptr<TOther>&& sp ) noexcept : m_shptr(std::move(sp)) { check_nonnull(); }

    template<class TOther>
    ncconstexpr17 shared_obj( const std::shared_ptr<TOther>& sp ) noexcept : m_shptr(sp) { check_nonnull(); }

    template<class TOther>
    ncconstexpr17 shared_obj& operator=( std::shared_ptr<TOther>&& sp ) noexcept
    {
      m_shptr = std::move(sp);
      check_nonnull();
      return *this;
    }

    template<class TOther>
    ncconstexpr17 shared_obj& operator=( const std::shared_ptr<TOther>& sp ) noexcept
    {
      m_shptr = sp;
      check_nonnull();
      return *this;
    }

    template<class TOther>
    ncconstexpr17 shared_obj( const shared_obj<TOther>& o ) noexcept : m_shptr(o.getsp()) {}

    template<class TOther>
    ncconstexpr17 operator std::shared_ptr<TOther>() noexcept { return m_shptr; }

    template<class TOther>
    constexpr operator std::shared_ptr<const TOther>() const noexcept { return m_shptr; }

    template<class TDerived>
    std::shared_ptr<TDerived> tryDynCast() { return std::dynamic_pointer_cast<TDerived>(this->m_shptr); }

    template<class TDerived>
    shared_obj<TDerived> dynCast()
    {
      auto sp = tryDynCast<TDerived>();
      if ( sp==nullptr )
        NCRYSTAL_THROW(LogicError,"shared_obj::dynCast ERROR: dynamic cast failed");
      return { sp };
    }

    //Efficient creation (same as global makeSO<T>(...) function):
    template<typename ...Args>
    static shared_obj<T> make( Args&& ...args )
    {
      return shared_obj<T>( std::make_shared<T>(std::forward<Args>(args)...), guaranteed_non_null_t() );
    }

  private:
    ncconstexpr17 void check_nonnull() const
    {
      if (!m_shptr)
        NCRYSTAL_THROW(BadInput,"Attempt to initialise shared_obj<T> object with null pointer is illegal");
    }
    std::shared_ptr<T> m_shptr;
    friend class shared_obj<element_type_nocv>;
    friend class shared_obj<const element_type_nocv>;
    //Efficient constructor for make(...) function:
    struct guaranteed_non_null_t {};
    constexpr shared_obj( std::shared_ptr<T>&& sp, guaranteed_non_null_t ) noexcept : m_shptr(std::move(sp)) {}
  public:
    //Don't use this! Only for casting functions above:
    ncconstexpr17 std::shared_ptr<T>&& detail_consume_shptr() && noexcept { return std::move(m_shptr); }
  };

  //Efficient/safe creation of shared_object's (similar to std::make_unique/std::make_shared):
  template<typename T, typename ...Args>
  inline shared_obj<T> makeSO( Args&& ...args )
  {
    return shared_obj<T>::make( std::forward<Args>(args)... );
  }

  template <class Derived>
  struct EnableSharedFromThis : std::enable_shared_from_this<Derived> {
    //Get smart pointer to current object (object MUST be already managed by
    //shared pointer, e.g. created with NCrystal::makeSO or std::make_shared):
    shared_obj<Derived> shared_obj_from_this() {
#  if __cplusplus >= 201703L
      auto shptr = this->weak_from_this().lock();
      if (!shptr)
        NCRYSTAL_THROW( LogicError, "shared_obj_from_this() does not work since this"
                        " instance is not managed by std::shared_ptr / NCrystal::shared_obj" );
      return shptr;
#else
      return this->shared_from_this();
#endif
    }
    shared_obj<const Derived> shared_obj_from_this() const {
#  if __cplusplus >= 201703L
      auto shptr = this->weak_from_this().lock();
      if (!shptr)
        NCRYSTAL_THROW( LogicError, "shared_obj_from_this() does not work since this"
                        " instance is not managed by std::shared_ptr / NCrystal::shared_obj" );
      return shptr;
#else
      return this->shared_from_this();
#endif
    }

  };

  //Aligned allocation suitable for any alignment size. Returned memory must be
  //eventually released by call to std::free. Throws std::bad_alloc in case the
  //allocation fails.
  void * alignedAlloc( std::size_t alignment, std::size_t size );

  //Type-safe version of alignedAlloc (returned memory must be free'd by
  //call to std::free, not delete):
  template<class TValue>
  inline TValue* alignedAlloc( std::size_t number_of_objects );

}

#if __cplusplus < 201402L
//Make sure we can use std::make_unique from C++14 even in C++11 code
namespace std {
  template<typename T, typename ...Args>
  inline std::unique_ptr<T> make_unique( Args&& ...args )
  {
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
  }
}
#endif


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCrystal {
  namespace detail {
    //For alignment > alignof(std::max_align_t):
    void * bigAlignedAlloc( std::size_t alignment, std::size_t size );
  }

  //Aligned allocation suitable for any alignment size. Returned memory must be
  //eventually released by call to std::free. Throws std::bad_alloc in case the
  //allocation fails.
  inline void * alignedAlloc( std::size_t alignment, std::size_t size ) {
    nc_assert( size>0 );//size 0 malloc is UB
    nc_assert( alignment>0 );//duh!
    nc_assert( (alignment & (alignment - 1)) == 0 );//checks if (non-zero) alignment is power of 2.

#if defined(__GNUC__) && (__GNUC__ > 11 )
    if ( false ) {//Very temporary workaround for gcc 12! This is of course pretty nasty!
#else

#  if defined(__GNUC__) && (__GNUC__*1000+__GNUC_MINOR__)<4009
    //gcc didn't add max_align_t to std:: until gcc 4.9
    constexpr auto alignof_max_align_t = alignof(max_align_t);
#  else
    constexpr auto alignof_max_align_t = alignof(std::max_align_t);
#  endif


    if ( alignment <= alignof_max_align_t ) {
#endif
      //std::malloc is supposed to work in this case (and std::aligned_alloc
      //might NOT work!):
      void * result = std::malloc(size);
      if ( result )
        return result;
      throw std::bad_alloc();
    } else {
      return detail::bigAlignedAlloc(alignment,size);
    }
  }

  template<class TValue>
  inline TValue* alignedAlloc( std::size_t number_of_objects )
  {
    return static_cast<TValue*>(alignedAlloc( alignof(TValue), number_of_objects * sizeof(TValue) ));
  }

}

#endif

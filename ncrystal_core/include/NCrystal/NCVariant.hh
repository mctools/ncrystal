#ifndef NCrystal_Variant_hh
#define NCrystal_Variant_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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

#include "NCrystal/NCDefs.hh"

namespace NCrystal {

  template <class T1, class T2>
  class Variant {
    //Very simplistic alternative to std::variant which only supports two
    //variant types and which optionally can hold no type.

    //Todo: some of the methods below could be noexcept depending on methods of T1 and T2
  public:

    Variant() noexcept;
    Variant( const T1& );
    Variant( T1&& );
    Variant( const T2& );
    Variant( T2&& );
    ~Variant();

    void clear();
    constexpr bool empty() const noexcept;

    template<class T> constexpr bool has_value() const noexcept;
    template<class T> T& get() ncnoexceptndebug;
    template<class T> const T& get() const ncnoexceptndebug;
    template<class T, typename ...Args> void emplace( Args&& ... );

    Variant( const Variant& );
    Variant& operator=( const Variant& );
    Variant( Variant&& );
    Variant& operator=( Variant&& );

    //Copy/move (NB: moved from Variant is always empty)
    Variant& operator=( const T1& );
    Variant& operator=( T1&& );
    Variant& operator=( const T2& );
    Variant& operator=( T2&& );

  private:
    static constexpr unsigned MAXSIZE = ( sizeof(T1) > sizeof(T2) ? sizeof(T1) : sizeof(T2) );
    //gcc 4.8 has problems with:
    //  static constexpr unsigned MAXALIGN = ( alignof(T1) > alignof(T2) ? alignof(T1) : alignof(T2) );
    //  alignas(MAXALIGN) uint8_t m_data[MAXSIZE];
    //So do it like:
    alignas(alignof(T1) > alignof(T2) ? alignof(T1) : alignof(T2)) uint8_t m_data[MAXSIZE];
    enum class Content { HasT1, HasT2, Empty };
    Content m_content = Content::Empty;
    struct Impl;
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCrystal {

  template <class T1, class T2>
  struct Variant<T1,T2>::Impl {

    static_assert(!std::is_same<T1,T2>::value,"Variant class must be given two distinct types");

    //For simplicity, explicitly disallow esoteric types:
    static_assert(!std::is_volatile<T1>::value,"Variant class does not support volative types");
    static_assert(!std::is_reference<T1>::value,"Variant class does not support reference types");
    static_assert(!std::is_member_pointer<T1>::value,"Variant class does not support pointer types");
    static_assert(!std::is_member_object_pointer<T1>::value,"Variant class does not support pointer types");
    static_assert(!std::is_member_function_pointer<T1>::value,"Variant class does not support pointer types");
    static_assert(!std::is_pointer<T1>::value,"Variant class does not support pointer types");
    static_assert(!std::is_array<T1>::value,"Variant class does not support array types");

    static_assert(!std::is_volatile<T2>::value,"Variant class does not support volative types");
    static_assert(!std::is_reference<T2>::value,"Variant class does not support reference types");
    static_assert(!std::is_member_pointer<T2>::value,"Variant class does not support pointer types");
    static_assert(!std::is_member_object_pointer<T2>::value,"Variant class does not support pointer types");
    static_assert(!std::is_member_function_pointer<T2>::value,"Variant class does not support pointer types");
    static_assert(!std::is_pointer<T2>::value,"Variant class does not support pointer types");
    static_assert(!std::is_array<T2>::value,"Variant class does not support array types");

    //gcc 4.8 -Wstrict-aliasing will falsely warn if we directly cast and
    //dereference &m_data[0]. Using these functions to first convert to a void
    //pointer seems to shut it up:

    static constexpr void * data( Variant* v ) noexcept
    {
      return reinterpret_cast<void*>(&(v->m_data[0]));
    }

    static constexpr const void * data( const Variant* v ) noexcept
    {
      return reinterpret_cast<const void*>(&(v->m_data[0]));
    }

  };


  template <class T1, class T2>
  template <class T>
  inline constexpr bool Variant<T1,T2>::has_value() const noexcept
  {
    static_assert(std::is_same<T,T1>::value||std::is_same<T,T2>::value,
                  "Provided type is not one of the variants supported for this class");
    return m_content == ( std::is_same<T,T1>::value ? Content::HasT1 : Content::HasT2 );
  }

  template <class T1, class T2>
  template<class T>
  inline T& Variant<T1,T2>::get() ncnoexceptndebug {
#ifndef NDEBUG
    if (!has_value<T>())
      NCRYSTAL_THROW(LogicError,"Invalid variant access (does not currently have the requested type)");
#endif
    return *reinterpret_cast<T*>(Impl::data(this));
  }

  template <class T1, class T2>
  template<class T>
  inline const T& Variant<T1,T2>::get() const ncnoexceptndebug {
#ifndef NDEBUG
    if (!has_value<T>())
      NCRYSTAL_THROW(LogicError,"Invalid variant access (does not currently have the requested type)");
#endif
    return *reinterpret_cast<const T*>(Impl::data(this));
  }

  template <class T1, class T2>
  template<class T, typename ...Args>
  inline void Variant<T1,T2>::emplace( Args&& ...args )
  {
    static_assert(std::is_same<T,T1>::value||std::is_same<T,T2>::value,
                  "Provided type is not one of the variants supported for this class");
    if (m_content!=Content::Empty)
      clear();
    new(Impl::data(this)) T(std::forward<Args>(args)...);
    m_content = ( std::is_same<T,T1>::value ? Content::HasT1 : Content::HasT2 );
  }

  template <class T1, class T2>
  inline Variant<T1,T2>::Variant( const Variant& o )
    : Variant()
  {
    *this = o;
  }

  template <class T1, class T2>
  inline Variant<T1,T2>& Variant<T1,T2>::operator=( const Variant& o )
  {
    clear();
    if ( o.has_value<T1>() ) {
      emplace<T1>( o.get<T1>() );
    } else if ( o.has_value<T2>() ) {
      emplace<T2>( o.get<T2>() );
    }
    return *this;
  }

  template <class T1, class T2>
  inline Variant<T1,T2>::Variant( Variant&& o )
    : Variant()
  {
    *this = std::move(o);
  }

  template <class T1, class T2>
  inline Variant<T1,T2>&  Variant<T1,T2>::operator=( Variant&& o )
  {
    clear();
    if ( o.has_value<T1>() ) {
      emplace<T1>( std::move(o.get<T1>()));
      o.clear();
    } else if ( o.has_value<T2>() ) {
      emplace<T2>( std::move(o.get<T2>()));
      o.clear();
    }
    return *this;
  }


  template <class T1, class T2>
  inline Variant<T1,T2>::Variant() noexcept : m_content(Content::Empty) {}

  template <class T1, class T2>
  inline Variant<T1,T2>::Variant( const T1& t1 ) : m_content(Content::HasT1) { new(Impl::data(this)) T1(t1); }

  template <class T1, class T2>
  inline Variant<T1,T2>::Variant( T1&& t1 ) : m_content(Content::HasT1) { new(Impl::data(this)) T1(std::move(t1)); }

  template <class T1, class T2>
  inline Variant<T1,T2>::Variant( const T2& t2 ) : m_content(Content::HasT2) { new(Impl::data(this)) T2(t2); }

  template <class T1, class T2>
  inline Variant<T1,T2>::Variant( T2&& t2 ) : m_content(Content::HasT2) { new(Impl::data(this)) T2(std::move(t2)); }

  template <class T1, class T2>
  inline constexpr bool Variant<T1,T2>::empty() const noexcept { return m_content == Content::Empty; }

  template <class T1, class T2>
  inline void Variant<T1,T2>::clear()
  {
    if ( m_content == Content::Empty )
      return;
    if ( m_content==Content::HasT1 )
      reinterpret_cast<T1*>(Impl::data(this))->~T1();
    else
      reinterpret_cast<T2*>(Impl::data(this))->~T2();
    m_content = Content::Empty;
  }

  template <class T1, class T2>
  inline Variant<T1,T2>::~Variant()
  {
    clear();
  }

  template <class T1, class T2>
  inline Variant<T1,T2>& Variant<T1,T2>::operator=( const T1& o )
  {
    emplace<T1>(o);
    return *this;
  }

  template <class T1, class T2>
  inline Variant<T1,T2>& Variant<T1,T2>::operator=( T1&& o )
  {
    emplace<T1>(std::move(o));
    return *this;
  }

  template <class T1, class T2>
  inline Variant<T1,T2>& Variant<T1,T2>::operator=( const T2& o )
  {
    emplace<T2>(o);
    return *this;
  }


  template <class T1, class T2>
  inline Variant<T1,T2>& Variant<T1,T2>::operator=( T2&& o )
  {
    emplace<T2>(std::move(o));
    return *this;
  }

}

#endif

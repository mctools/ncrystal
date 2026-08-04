#ifndef NCrystal_TinyVector_hh
#define NCrystal_TinyVector_hh

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

#ifndef NCrystal_Defs_hh
#  include "NCrystal/core/NCDefs.hh"
#endif
#include <initializer_list>
#include <type_traits>

namespace NCRYSTAL_NAMESPACE {

  //This file provides a TinyVector container class, which like SmallVector is
  //designed to be a drop-in for std::vector. Unlike SmallVector, however, it
  //ONLY has local storage and thus comes with a maximum size.

  template<class TValue, std::size_t NMAX>
  class TinyVector
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

    using size_type = decltype(NMAX);
    static constexpr size_type nmax = NMAX;

    ///////////////////////////////////////////////////////////////////////////
    //Construction:
    TinyVector() = default;
    ~TinyVector() noexcept;

    template <class TIter>
    TinyVector( TIter it_begin, TIter it_end )
      : TinyVector() { setByCopy(it_begin,it_end); }

    ///////////////////////////////////////////////////////////////////////////
    //Allow copy/move semantics (might be inefficient for large vectors!)
    TinyVector( const TinyVector& );
    TinyVector& operator=( const TinyVector& );
    TinyVector( TinyVector&& );
    TinyVector& operator=( TinyVector&& );

    ncconstexpr17 TValue* data() noexcept;
    ncconstexpr17 const TValue* data() const noexcept;
    ncconstexpr17 TValue* begin() noexcept;
    ncnodiscard17 ncconstexpr17 TValue* end() noexcept;
    constexpr const TValue* begin() const noexcept;
    ncnodiscard17 constexpr const TValue* end() const noexcept;
    constexpr const TValue* cbegin() const noexcept { return begin(); }
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
    ncnodiscard17 constexpr size_type size() const noexcept { return m_size; }
    ncnodiscard17 constexpr bool empty() const noexcept { return m_size==0; }
    constexpr size_type capacity() const noexcept { return NMAX; }
    void shrink_to_fit() const noexcept {}
    void reserve( size_type ) const noexcept {}
    void reserve_hint( size_type ) const noexcept {}

    //fixme: remove?
    constexpr bool isLocalStorage() const noexcept { return true; }
    static constexpr bool isFastAccess() noexcept { return true; }

    ///////////////////////////////////////////////////////////////////////////
    //Insert or remove elements:
    void push_back( const TValue& value );
    void push_back( TValue&& value );
    template<typename ...Args>
    TValue& emplace_back( Args&& ... );
    void pop_back() ncnoexceptndebug;
    void clear() noexcept;

    ///////////////////////////////////////////////////////////////////////////
    //Comparison (first on size, then element-wise - this is different than
    //std::vector):
    bool operator==( const TinyVector& ) const noexcept;
    bool operator<( const TinyVector& ) const noexcept;

    ///////////////////////////////////////////////////////////////////////////
    //Various ways to set contents:
    TinyVector( std::initializer_list<TValue> );
    template <class TIter>
    void setByCopy( TIter it_begin, TIter it_end );
    template <class TIter>
    void setByMove( TIter it_begin, TIter it_end );
    void swap( TinyVector& o );

    ///////////////////////////////////////////////////////////////////////////
    //TinyVector does not have constructors like std::vector(count) or
    //std::vector(count,value). Deleting these explicitly gives better error
    //diagnostics for the user (specifying various bit widths for robustness):
    TinyVector( std::uint8_t ) = delete;
    TinyVector( std::uint16_t ) = delete;
    TinyVector( std::uint32_t ) = delete;
    TinyVector( std::uint64_t ) = delete;
    TinyVector( std::int8_t ) = delete;
    TinyVector( std::int16_t ) = delete;
    TinyVector( std::int32_t ) = delete;
    TinyVector( std::int64_t ) = delete;

    //To make std::swap work (also has additional global functions inlined
    //below):
    friend void swap( TinyVector& a, TinyVector& b ) noexcept { a.swap(b); }

  private:
    alignas(TValue) uint8_t m_data[ NMAX * sizeof(TValue) ];
    size_type m_size = 0;
  };

}

////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {


  template<class TValue, std::size_t NMAX>
  inline TinyVector<TValue,NMAX>::~TinyVector() noexcept
  {
    clear();
  }

  template<class TValue, std::size_t NMAX>
  inline ncconstexpr17 TValue* TinyVector<TValue,NMAX>::data() noexcept
  { return reinterpret_cast<TValue*>(&m_data[0]); }

  template<class TValue, std::size_t NMAX>
  inline ncconstexpr17 const TValue* TinyVector<TValue,NMAX>::data() const noexcept
  { return reinterpret_cast<const TValue*>(&m_data[0]); }

  template<class TValue, std::size_t NMAX>
  inline ncconstexpr17 TValue* TinyVector<TValue,NMAX>::begin() noexcept { return data(); }

  template<class TValue, std::size_t NMAX>
  ncnodiscard17 inline ncconstexpr17 TValue* TinyVector<TValue,NMAX>::end() noexcept { return data() + size(); }

  template<class TValue, std::size_t NMAX>
  inline constexpr const TValue* TinyVector<TValue,NMAX>::begin() const noexcept { return data(); }

  template<class TValue, std::size_t NMAX>
  ncnodiscard17 inline constexpr const TValue* TinyVector<TValue,NMAX>::end() const noexcept { return data() + size(); }



  template<class TValue, std::size_t NMAX>
  ncnodiscard17 inline ncconstexpr17 const typename TinyVector<TValue,NMAX>::value_type&
  TinyVector<TValue,NMAX>::operator[]( size_type i ) const noexcept
  {
    assert(i<m_size);
    return *std::next(data(),i);
  }

  template<class TValue, std::size_t NMAX>
  ncnodiscard17 inline ncconstexpr17 typename TinyVector<TValue,NMAX>::value_type&
  TinyVector<TValue,NMAX>::operator[]( size_type i ) noexcept
  {
    assert(i<m_size);
    return *std::next(data(),i);
  }

  template<class TValue, std::size_t NMAX>
  ncnodiscard17 inline ncconstexpr17 const typename TinyVector<TValue,NMAX>::value_type&
  TinyVector<TValue,NMAX>::at( size_type i ) const
  {
    if ( i >= m_size )
      throw std::out_of_range("TinyVector::at(): index out of out of range");
    return *std::next(data(),i);
  }

  template<class TValue, std::size_t NMAX>
  ncnodiscard17 inline ncconstexpr17 typename TinyVector<TValue,NMAX>::value_type&
  TinyVector<TValue,NMAX>::at( size_type i )
  {
    if ( i >= m_size )
      throw std::out_of_range("TinyVector::at(): index out of out of range");
    return *std::next(data(),i);
  }

  template<class TValue, std::size_t NMAX>
  template<typename ...Args>
  inline TValue& TinyVector<TValue,NMAX>::emplace_back( Args&& ...args ) {
    nc_assert( size() < capacity() );
    TValue * newobjaddr = this->end();
    new( (void*)(newobjaddr) ) TValue(std::forward<Args>(args)...);
    ++(this->m_size);
    return *newobjaddr;
  }

  template<class TValue, std::size_t NMAX>
  inline void TinyVector<TValue,NMAX>::push_back( const TValue& value )
  {
    this->emplace_back(value);
  }

  template<class TValue, std::size_t NMAX>
  inline void TinyVector<TValue,NMAX>::push_back( TValue&& value )
  {
    this->emplace_back(std::move(value));
  }

  template<class TValue, std::size_t NMAX>
  inline void TinyVector<TValue,NMAX>::clear() noexcept
  {
    static_assert(std::is_nothrow_destructible<TValue>::value,
                  "TinyVector can only keep objects with noexcept destructors");
    if ( std::is_trivially_destructible<TValue>::value ) {
      m_size = 0;
      return;
    }
    if ( m_size == 0 )
      return;
    auto it = begin();
    auto itE = end();
    m_size = 0;
    for ( ; it!=itE; ++it )
      it->~TValue();
  }

  template<class TValue, std::size_t NMAX>
  inline void TinyVector<TValue,NMAX>::pop_back() ncnoexceptndebug
  {
    nc_assert( m_size > 0 );
    auto it = std::prev(end());
    static_assert(std::is_nothrow_destructible<TValue>::value,
                  "TinyVector can only keep objects with noexcept destructors");
    --m_size;
    if ( !std::is_trivially_destructible<TValue>::value )
      it->~TValue();
  }

  template<class TValue, std::size_t NMAX>
  inline TinyVector<TValue,NMAX>::TinyVector( std::initializer_list<TValue> l)
    : TinyVector()
  {
    nc_assert( ( l.size() <= capacity() )&&"insufficient TinyVector capacity" );
    auto it = begin();
    for ( auto&& e : l )
      new(it++) TValue(std::move(e));
    m_size = l.size();
  }

  template<class TValue, std::size_t NMAX>
  inline bool TinyVector<TValue,NMAX>::operator==( const TinyVector& o ) const noexcept
  {
    if ( this->m_size != o.m_size )
      return false;
    if ( this == &o || this->m_size == 0 )
      return true;
    auto it = this->begin();
    auto itE = it + this->m_size;
    auto itO = o.begin();
    for ( ; it != itE; ++it, ++itO )
      if ( ! ( *it == *itO ) )
        return false;
    return true;
  }

  template<class TValue, std::size_t NMAX>
  inline bool TinyVector<TValue,NMAX>::operator<( const TinyVector& o ) const noexcept
  {
    if ( this->m_size != o.m_size )
      return this->m_size < o.m_size;
    if ( this == &o || this->m_size == 0 )
      return false;
    auto it = this->begin();
    auto itE = it + this->m_size;
    auto itO = o.begin();
    for ( ; it != itE; ++it, ++itO )
      if ( ! ( *it == *itO ) )
        return *it < *itO;
    return false;
  }

  template<class TValue, std::size_t NMAX>
  inline TinyVector<TValue,NMAX>::TinyVector( TinyVector&& o )
    : TinyVector()
  {
    *this = std::move(o);
  }

  template<class TValue, std::size_t NMAX>
  inline TinyVector<TValue,NMAX>::TinyVector( const TinyVector& o )
    : TinyVector()
  {
    *this = o;
  }

  template<class TValue, std::size_t NMAX>
  inline TinyVector<TValue,NMAX>& TinyVector<TValue,NMAX>::operator=( const TinyVector& o )
  {
    if ( this == &o )
      return *this;
    clear();
    for ( const auto& e : o )
      push_back( e );
    return *this;
  }

  template<class TValue, std::size_t NMAX>
  inline TinyVector<TValue,NMAX>& TinyVector<TValue,NMAX>::operator=( TinyVector&& o )
  {
    if ( this == &o )
      return *this;
    clear();
    for ( auto& e : o )
      emplace_back( std::move(e) );
    o.clear();
    return *this;
  }

  template<class TValue, std::size_t NMAX>
  template <class TIter>
  inline void TinyVector<TValue,NMAX>::setByCopy( TIter it_begin, TIter it_end )
  {
    clear();
    for ( ; it_begin!=it_end; ++it_begin )
      push_back( *it_begin );
  }

  template<class TValue, std::size_t NMAX>
  template <class TIter>
  inline void TinyVector<TValue,NMAX>::setByMove( TIter it_begin, TIter it_end )
  {
    clear();
    for ( ; it_begin!=it_end; ++it_begin )
      emplace_back( std::move(*it_begin) );
  }

  template<class TValue, std::size_t NMAX>
  inline void TinyVector<TValue,NMAX>::swap( TinyVector& o )
  {
    //NB: Just a simple and non-optimal implementation for now:
    TinyVector tmp( std::move(o) );
    o = std::move(*this);
    *this = std::move(tmp);
  }

  static_assert(std::is_copy_constructible<TinyVector<int,1>>::value,"");
  static_assert(std::is_copy_assignable<TinyVector<int,1>>::value,"");
  static_assert(std::is_nothrow_default_constructible<TinyVector<int,1>>::value,"");
}

namespace std {
  //Would like to avoid injections in std::, but it seems to not work reliably without it:
  template<class TValue, std::size_t NMAX>
  inline void swap(NCrystal::TinyVector<TValue,NMAX>& a,
                   NCrystal::TinyVector<TValue,NMAX>& b)
  {
    a.swap(b);
  }
}

#endif

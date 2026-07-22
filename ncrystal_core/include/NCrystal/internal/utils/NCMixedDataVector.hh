#ifndef NCrystal_MixedDataVector_hh
#define NCrystal_MixedDataVector_hh

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

namespace NCRYSTAL_NAMESPACE {

  class MixedDataVector : private MoveOnly {
    std::vector<std::uint8_t> m_data;
  public:

  //////////////////////////////////////////////////////////////////////////
  // Low-level utility class for packing of heterogeneous trivial datatypes.

    MixedDataVector( std::size_t reserve_num_bytes = 0 );

    MixedDataVector( std::vector<std::uint8_t>&& data )
      : m_data( std::move(data) ) {}
    const std::vector<std::uint8_t>& data() const { return m_data; }
    std::vector<std::uint8_t> stealData() && { return std::move(m_data); }

    void clear() { m_data.clear(); }
    MixedDataVector clone() const;
    void reserve_bytes( std::size_t nbytes ) { m_data.reserve( nbytes ); }
    void shrink_to_fit() { m_data.shrink_to_fit(); }
    std::size_t size_bytes() const { return m_data.size(); }
    void swap( MixedDataVector& o ) { m_data.swap( o.m_data ); }

    //Appending and extracting data is implemented via std::memcpy, for
    //portability and alignment safety.

    template <class T>
    std::size_t append( const T& value )
    {
      static_assert(std::is_default_constructible<T>::value, "");
      static_assert(std::is_trivially_constructible<T>::value, "");
      static_assert(std::is_trivially_copyable<T>::value, "");
      static_assert(std::is_trivially_destructible<T>::value, "");
      const std::size_t pos = m_data.size();
      const std::size_t newsize = pos + sizeof(T);
      if ( m_data.capacity() < newsize )
        m_data.reserve( newsize*2 );
      m_data.resize( newsize );
      std::memcpy( m_data.data() + pos, &value, sizeof(T));
      return pos;
    }

    template <class T>
    T extract( std::size_t pos ) const
    {
      static_assert(std::is_default_constructible<T>::value, "");
      static_assert(std::is_trivially_constructible<T>::value, "");
      static_assert(std::is_trivially_copyable<T>::value, "");
      static_assert(std::is_trivially_destructible<T>::value, "");
      nc_assert( pos + sizeof(T) <= m_data.size() );
      T res;
      std::memcpy(&res, m_data.data() + pos, sizeof(T));
      return res;
    }
  };

}

////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::MixedDataVector NCrystal::MixedDataVector::clone() const
{
  MixedDataVector o;
  o.m_data = m_data;
  return o;
}

inline NCrystal::MixedDataVector::MixedDataVector( std::size_t reserve_num_bytes )
{
  if ( reserve_num_bytes )
    reserve_bytes( reserve_num_bytes );
}

#endif

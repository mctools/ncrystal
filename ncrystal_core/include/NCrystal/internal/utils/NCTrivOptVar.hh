#ifndef NCrystal_TrivOptVar_hh
#define NCrystal_TrivOptVar_hh

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

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Optional / Variant for trivial types that are themselves trivial. //
//                                                                    //
////////////////////////////////////////////////////////////////////////

//Fixme: Cleanup (move stuff to inline section) and unit test extensively.

namespace NCRYSTAL_NAMESPACE {

  namespace Utils {
    namespace detail {
      struct TV_ph2 {};
      struct TV_ph3 {};
    }

    template<class T1, class T2 = detail::TV_ph2, class T3 = detail::TV_ph3>
    class TrivialVariant {

      // A Variant class for N trivially destructible and copyable types, that
      // can be serialised with a simple bit-wise copy and do not need actual
      // destructors to run. It can be empty.
      //
      // For now we support N=1, N=2 and N=3, but it can easily be extended to
      // higher numbers.

      static_assert(std::is_trivially_destructible<T1>::value,"");
      static_assert(std::is_trivially_copyable<T1>::value,"");
      static_assert(std::is_trivially_destructible<T2>::value,"");
      static_assert(std::is_trivially_copyable<T2>::value,"");
      static_assert(std::is_trivially_destructible<T3>::value,"");
      static_assert(std::is_trivially_copyable<T3>::value,"");
      static_assert(!std::is_same<T1,T2>::value,"");
      static_assert(!std::is_same<T1,T3>::value,"");
      static_assert(!std::is_same<T2,T3>::value,"");

      constexpr static std::size_t maxsize = ncconstexpr_max( sizeof(T1),
                                                              sizeof(T2),
                                                              sizeof(T3) );
      constexpr static std::size_t maxalign = ncconstexpr_lcm( alignof(T1),
                                                               alignof(T2),
                                                               alignof(T3) );
      alignas(maxalign) char m_data[maxsize];
      int m_content = 0;//empty

    public:
      TrivialVariant() = default;
      ~TrivialVariant() = default;
      TrivialVariant( const TrivialVariant& ) = default;
      TrivialVariant& operator=( const TrivialVariant& ) = default;

      TrivialVariant( const T1& o ) { *this = o; }
      TrivialVariant( const T2& o ) { *this = o; }
      TrivialVariant( const T3& o ) { *this = o; }

      TrivialVariant( NullOptType ) {}
      TrivialVariant& operator=( const NullOptType& ) { m_content = 0; return *this; }

      void clear() { m_content = 0; }
      void reset() { clear(); }

      template<class Tn>
      bool has_value() const noexcept {
        static_assert( std::is_same<Tn,T1>::value
                       || std::is_same<Tn,T2>::value
                       || std::is_same<Tn,T3>::value, "" );
        constexpr int valid = ( std::is_same<Tn,T1>::value ? 1
                                : ( std::is_same<Tn,T2>::value ? 2 : 3 ) );
        return m_content == valid;
      }

      template<class Tn>
      const Tn& value() const
      {
        nc_assert_always(has_value<Tn>());//fixme _always
        return *reinterpret_cast<const Tn*>(m_data);
      }

      template<class Tn>
      Tn& value()
      {
        nc_assert_always(has_value<Tn>());//fixme _always
        return *reinterpret_cast<Tn*>(m_data);
      }

      bool empty() const { return m_content == 0; }

      template<class Tn>
      TrivialVariant& operator=( const Tn& o )
      {
        if ( m_data != (char*)(&o) )
          std::memcpy( m_data, (char*)&o, sizeof(Tn) );
        static_assert( std::is_same<Tn,T1>::value
                       || std::is_same<Tn,T2>::value
                       || std::is_same<Tn,T3>::value, "" );
        constexpr int valid = ( std::is_same<Tn,T1>::value ? 1
                                : ( std::is_same<Tn,T2>::value ? 2 : 3 ) );
        m_content = valid;
        return *this;
      }
    };

    static_assert(std::is_trivially_destructible<TrivialVariant<double,int,char>>::value,"");
    static_assert(std::is_trivially_copyable<TrivialVariant<double,int,char>>::value,"");

    template<class TData>
    class TrivialOptional {

      // An Optional class for trivially destructible and copyable types, that
      // can be serialised with a simple bit-wise copy and do not need actual
      // destructors to run..

      static_assert(std::is_trivially_destructible<TData>::value,"");
      static_assert(std::is_trivially_copyable<TData>::value,"");
      alignas(TData) char m_data[sizeof(TData)];
      bool m_hasData = false;
    public:

      using data_t = TData;

      TrivialOptional() = default;
      ~TrivialOptional() = default;
      TrivialOptional( const TrivialOptional& ) = default;
      TrivialOptional& operator=( const TrivialOptional& ) = default;

      TrivialOptional( const TData& o ) { *this = o; }
      TrivialOptional( NullOptType ) {}
      TrivialOptional& operator=( const NullOptType& ) { m_hasData = false; return *this; }

      void clear() { m_hasData = false; }
      void reset() { clear(); }
      bool has_value() const { return m_hasData; }
      const TData& value() const
      {
        nc_assert(has_value());
        return *reinterpret_cast<const TData*>(m_data);
      }
      TData& value()
      {
        nc_assert(has_value());
        return *reinterpret_cast<TData*>(m_data);
      }

      bool empty() const { return !m_hasData; }

      TrivialOptional& operator=( const TData& o )
      {
        if ( m_data != (char*)(&o) )
          std::memcpy( m_data, (char*)&o, sizeof(TData) );
        m_hasData = true;
        return *this;
      }

      template<typename... Args>
      void emplace( Args&& ...args )
      {
        new(&m_data) TData(std::forward<Args>(args)...);
        m_hasData = true;
      }

    };
    static_assert(std::is_trivially_destructible<TrivialOptional<double>>::value,"");
    static_assert(std::is_trivially_copyable<TrivialOptional<double>>::value,"");

  }

}

#endif

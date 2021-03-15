#ifndef NCrystal_LegacyMem_hh
#define NCrystal_LegacyMem_hh

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

/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// Legacy memory related classes and functions. Provided here for now for      //
// backwards compatibility with code using NCrystal. Everything is implemented //
// in the NCrystal::Legacy namespace, which is for now injected into the main  //
// NCrystal namespace. This will eventually change!                            //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/NCException.hh"
#include <memory>
#include <type_traits>
#include <utility>//std::move
#include <functional>
#include <atomic>

namespace NCrystal {

  namespace Legacy {

    //Base class for ref-counted objects, destructors of which should be protected
    //rather than public.

    class NCRYSTAL_API RCBase {
    public:
    unsigned refCount() const noexcept { return m_refCount; }
    void ref() const noexcept { ++m_refCount; }
    void unref() const noexcept
    {
      if (m_refCount.fetch_sub(1)==1)
        delete this;
    }

    //Careful with this one (only used by deprecated NCrystal factories!):
    void unrefNoDelete() const noexcept
    {
      --m_refCount;
    }

    //Monitor number of RCBase instances or enable dbg printouts. The lvl
    //parameter should be 0, 1 or 2. Default will be 0, or the value set in the
    //NCRYSTAL_DEBUGMEM env var:
    static long nInstances();
    static void enableMemDbg(int lvl);

    //No copy/move:
    RCBase( const RCBase& ) = delete;
    RCBase& operator=( const RCBase& ) = delete;
    RCBase( RCBase&& ) = delete;
    RCBase& operator=( RCBase&& ) = delete;

    protected:
    RCBase() noexcept;
    virtual ~RCBase();
    private:
    mutable std::atomic<unsigned> m_refCount;
    };

    template< class T >
    struct NCRYSTAL_API RCHolder {
      typedef T element_type;
      ncconstexpr17 RCHolder() = default;
      ncconstexpr17 explicit RCHolder( T* t ) : m_obj( t ) { if ( m_obj ) m_obj->ref(); }
      ~RCHolder() { if ( m_obj ) m_obj->unref(); }
      ncconstexpr17 RCHolder( const RCHolder & o ) : m_obj(o.m_obj) { if (m_obj) m_obj->ref(); }
      ncconstexpr17 T* obj() noexcept { return m_obj; }
      constexpr const T* obj() const noexcept  { return m_obj; }
      void clear() {
        const T* old = m_obj;
        m_obj = 0;//clear first in case unref throws
        if (old)
          old->unref();
      }
      RCHolder & operator= ( const RCHolder & o)
      {
        if (o.m_obj!=m_obj) {
          clear();
          m_obj = o.m_obj;
          if (m_obj)
            m_obj->ref();
        }
        return *this;
      }
      RCHolder & operator= ( T* newobj )
      {
        if (newobj!=m_obj) {
          clear();
          m_obj = newobj;
          if (m_obj)
            m_obj->ref();
        }
        return *this;
      }
      RCHolder( RCHolder&& o ) { std::swap(m_obj,o.m_obj); }
      RCHolder& operator=( RCHolder&& o ) { clear(); std::swap(m_obj,o.m_obj); return *this; }
      bool operator()() const { return m_obj!=0; }
      bool operator!() const { return m_obj==0; }

      //Assign, construct and compare with nullptr:
      RCHolder & operator= ( decltype(nullptr) ) { m_obj = nullptr; return *this; }
      RCHolder( decltype(nullptr) ) {}
      constexpr bool operator==( decltype(nullptr) ) const noexcept { return m_obj==nullptr; }
      constexpr bool operator!=( decltype(nullptr) ) const noexcept { return m_obj!=nullptr; }

      //Move/copy/assign to base class pointer should be ok (only compiles if T* can be assigned from TOther*):
      template<class TOther>
      ncconstexpr17 RCHolder& operator=( RCHolder<TOther>&& o ) noexcept {
        RCHolder<TOther> o_taken = std::move(o);
        m_obj = o_taken.obj();
        if (m_obj)
          m_obj->ref();
        return *this;
      }
      template<class TOther>
      ncconstexpr17 RCHolder( RCHolder<TOther>&& o ) noexcept {
        *this = std::move(o);
      }
      template<class TOther>
      ncconstexpr17 RCHolder& operator=( const RCHolder<TOther>& o ) noexcept {
        m_obj = o.obj();
        if (m_obj)
          m_obj->ref();
        return *this;
      }
      template<class TOther>
      ncconstexpr17 RCHolder( const RCHolder<TOther>& o ) noexcept {
        *this = o;
      }

      //Release obj without triggering deletion (careful, only used by deprecated NCrystal factories!):
      T * releaseNoDelete() noexcept { T* old = m_obj; if ( m_obj ) m_obj->unrefNoDelete(); m_obj = 0; return old; }

      template<class TOther>
      constexpr RCHolder<TOther> dyncast() const noexcept
      {
        return RCHolder<TOther>(dynamic_cast<TOther*>(m_obj));
      }

      ncconstexpr17 T* operator->() noexcept { return m_obj; }
      constexpr const T* operator->() const noexcept { return m_obj; }
      ncconstexpr17 T& operator*() ncnoexceptndebug { nc_assert(m_obj); return *m_obj; }
      ncconstexpr17 const T& operator*() const ncnoexceptndebug { nc_assert(m_obj); return *m_obj; }

    private:
    T* m_obj = nullptr;
    };

    typedef RCHolder<const RCBase> RCGuard;

    //get_pointer specialisation for RCHolder allows easy integration with boost:
    template< class T >
    T* get_pointer(const RCHolder<T>& r) { return const_cast<T*>(r.obj()); }

    //makeRC can be used for RCBase-derived objects, similarly to how one would
    //use make_unique/make_shared for std smart pointers:
    template<typename T, typename ...Args>
    inline RCHolder<T> makeRC( Args&& ...args )
    {
      return RCHolder<T>( new T( std::forward<Args>(args)... ) );
    }

    template<typename T, typename Tauto>
    inline RCHolder<T> castRC( RCHolder<Tauto> tin )
    {
      return RCHolder<T>(dynamic_cast<T*>(tin.obj()));
    }

    template<typename T, typename Tauto>
    inline RCHolder<T> static_castRC( RCHolder<Tauto> tin )
    {
      nc_assert(dynamic_cast<T*>(tin.obj()));
      return RCHolder<T>(static_cast<T*>(tin.obj()));
    }
  }
}

#endif

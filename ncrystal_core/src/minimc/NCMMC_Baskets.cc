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

#include "NCrystal/internal/minimc/NCMMC_Baskets.hh"
#include "NCMMC_BasketSrcFiller.hh"
#include "NCMMC_MemPool.hh"
#include "NCMMC_BasketUtils.hh"
#ifndef NCRYSTAL_DISABLE_THREADS
#  include <condition_variable>
#  include <mutex>
#endif

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    class detail::UBImpl {
    public:
      static void*& access_internal( Basket& b ) noexcept
      {
        return b.internal;
      }
      static const void* access_internalc( const Basket& b ) noexcept
      {
        return b.internal;
      }
      static bool& access_token_activeflag( WorkerToken& t ) noexcept
      {
        return t.m_active;
      }
      static Basket& access_token_basket_buf( WorkerToken& t ) noexcept
      {
        return t.m_basketbuf;
      }
      static Basket&
      token_basket_buf_refresh( WorkerToken& t,
                                BasketMgr& bmgr ) ncnoexceptndebug
      {
        Basket& b = t.m_basketbuf;
        if ( !b.valid() )
          b = bmgr.allocateBasket();
        if ( b.size()>0 )
          b.neutrons->nused = 0;
        nc_assert ( b.size() == 0 );
        return b;
      }
    };

    namespace {

      struct Basket_Basic final {
        NeutronBasket neutrons;
        BasketValBufInt nscat;
        BasketValBufInt nscat_inelas;
        BasketValBufDbl buf1;

        //Infrastructure:
        constexpr bool full() const noexcept { return neutrons.full(); }
        constexpr bool empty() const noexcept { return neutrons.empty(); }
        constexpr std::size_t size() const noexcept { return neutrons.size(); }

        void markAsMissedTarget( std::size_t i ) noexcept { nscat.data[i] = -1; }

        NeutronBasket& get_neutrons() { return neutrons; }
        const NeutronBasket& get_neutrons() const { return neutrons; }

        void init_extra( std::size_t i ) ncnoexceptndebug
        {
          nscat.data[i] = 0;
          nscat_inelas.data[i] = 0;
          //Although the engine is responsible for any contents in buf1, we
          //initialise it anyway to ensure no uninitialised memory is ever
          //copied around or assigned. We initialise it to 0.0
          buf1.data[i] = 0.0;
        }

        void validateIfDbg() const ncnoexceptndebug
        {
          BasketUtils::basket_validateIfDbg(neutrons);
        }

        void assignToUB( Basket& dst ) noexcept
        {
          dst.neutrons = &neutrons;
          dst.nscat = &nscat;
          dst.nscat_inelas = &nscat_inelas;
          dst.buf1 = &buf1;
        }

        void copy1( const Basket_Basic& o,
                    std::size_t i_o, std::size_t i ) ncnoexceptndebug
        {
          //NB: Does NOT update neutrons.nused!
          nc_assert(i_o < o.size());
          //can't check i<this->size(), since size might not have been updated
          //yet by caller.
          BasketUtils::basketfields_set( neutrons.fields,
                                         o.neutrons.fields, i_o, i );
          nscat.data[i] = o.nscat[i_o];
          nscat_inelas.data[i] = o.nscat_inelas[i_o];
          buf1.data[i] = o.buf1[i_o];
        }

        void appendN( const Basket_Basic& o,
                      std::size_t i_o, std::size_t n ) ncnoexceptndebug
        {
          nc_assert( n>0 );
          nc_assert( this != &o );
          nc_assert( n>=1 );
          nc_assert( this->size() + n <= basket_N );
          nc_assert( i_o + n <= o.size() );
          const std::size_t i = this->size();

          //copy all data:
          BasketUtils::memcpydata<int>( nscat.data + i,o.nscat.data + i_o,n );
          BasketUtils::memcpydata<int>( nscat_inelas.data + i,
                                        o.nscat_inelas.data + i_o, n );
          BasketUtils::memcpydata<double>( buf1.data + i,
                                           o.buf1.data + i_o, n );
          nc_assert( neutrons.nused + n <= basket_N );
          nc_assert( i_o + n <= o.neutrons.nused );
          BasketUtils::basketfields_setrange( neutrons.fields,
                                              o.neutrons.fields,
                                              i_o, neutrons.nused, n );
          nc_assert( neutrons.fields.x[i] == o.neutrons.fields.x[i_o] );
          nc_assert( neutrons.fields.y[i] == o.neutrons.fields.y[i_o] );
          //finally update size:
          neutrons.nused += n;
        }

        std::size_t append1( const Basket_Basic& o,
                             std::size_t i_o ) ncnoexceptndebug
        {
          nc_assert(!neutrons.full());
          std::size_t i = neutrons.nused++;
          copy1( o, i_o, i );
          return i;
        }
      };

      struct Basket_Extended final {

        //Basic info:
        Basket_Basic basic;

        //Extended info:
        NeutronBasketFields neutrons_initial;

        //Infrastructure:
        constexpr bool full() const noexcept { return basic.neutrons.full(); }
        constexpr bool empty() const noexcept { return basic.neutrons.empty(); }
        constexpr std::size_t size() const noexcept { return basic.neutrons.size(); }

        void markAsMissedTarget( std::size_t i ) noexcept { basic.markAsMissedTarget(i); }

        NeutronBasket& get_neutrons() { return basic.get_neutrons(); }
        const NeutronBasket& get_neutrons() const { return basic.get_neutrons(); }

        void init_extra( std::size_t i ) ncnoexceptndebug
        {
          basic.init_extra(i);
          //The neutrons_initial parameters are just a copy of the initial
          //values of the neutrons:
          BasketUtils::basketfields_set( neutrons_initial,
                                         basic.neutrons.fields, i, i );
          nc_assert( neutrons_initial.x[i]
                     == basic.neutrons.fields.x[i] );
          nc_assert( neutrons_initial.ekin[i]
                     == basic.neutrons.fields.ekin[i] );
        }

        void validateIfDbg() const ncnoexceptndebug
        {
          basic.validateIfDbg();
          BasketUtils::basketfields_validateIfDbg( neutrons_initial,
                                                   basic.size() );
        }

        void assignToUB( Basket& dst ) noexcept
        {
          basic.assignToUB( dst );
          dst.neutrons_initial = &neutrons_initial;
        }

        void copy1( const Basket_Extended& o,
                    std::size_t i_o, std::size_t i ) ncnoexceptndebug
        {
          //NB: Does NOT update neutrons.nused!
          basic.copy1( o.basic, i_o, i );
          BasketUtils::basketfields_set( neutrons_initial,
                                         o.neutrons_initial, i_o, i );
        }

        void appendN( const Basket_Extended& o,
                      std::size_t i_o, std::size_t n ) ncnoexceptndebug
        {
          const std::size_t this_size = basic.neutrons.nused;
          nc_assert( n > 0 );
          nc_assert( this_size + n <= basket_N );
          basic.appendN( o.basic, i_o, n );
          BasketUtils::basketfields_setrange( neutrons_initial,
                                              o.neutrons_initial,
                                              i_o, this_size, n );
          nc_assert( neutrons_initial.x[this_size] == o.neutrons_initial.x[i_o] );
          nc_assert( basic.neutrons.fields.x[this_size] == o.basic.neutrons.fields.x[i_o] );
        }

        std::size_t append1( const Basket_Extended& o,
                             std::size_t i_o ) ncnoexceptndebug
        {
          nc_assert(!full());
          std::size_t i = basic.neutrons.nused;
          copy1( o, i_o, i );
          ++basic.neutrons.nused;
          return i;
        }

      };

      static_assert( std::is_standard_layout<Basket_Basic>::value, "" );
      static_assert( std::is_standard_layout<Basket_Extended>::value, "" );
      static_assert( std::is_trivially_destructible<Basket_Basic>::value, "" );
      static_assert( std::is_trivially_destructible<Basket_Extended>::value, "" );

      inline void *& ub_internal( Basket& b ) noexcept
      {
        return detail::UBImpl::access_internal(b);
      }

      template<class TBasket>
      struct UBHImpl{
        //typedef templated internal data structures:
        using basket_t = TBasket;
        using basket_src_filler_t = BasketSrcFiller<basket_t>;
        using heapmem_t = HeapMem<alignof(TBasket),sizeof(TBasket)>;
        using heapmempool_t = HeapMemPool<heapmem_t,16>;

        static basket_t& real_basket( Basket& b )
        {
          auto rbptr = reinterpret_cast<basket_t*>
            ( detail::UBImpl::access_internal(b) );
          nc_assert( rbptr != nullptr );
          return *rbptr;
        }

        static const basket_t& real_basket( const Basket& b )
        {
          auto rbptr = reinterpret_cast<const basket_t*>
            ( detail::UBImpl::access_internalc(b) );
          nc_assert( rbptr != nullptr );
          return *rbptr;
        }

        static void assign_basket_parameters( Basket& dst,
                                              const basket_t& src_c ) noexcept
        {
          auto& src = *const_cast<basket_t*>(&src_c);
          src.assignToUB( dst );
        }

        class ViewBasketAsUB final {
          //Non-owning way to view a basket_t as a Basket. The main point
          //is to reset the internal pointer of the Basket object in the
          //destructor, before it is itself destructed.
          Basket m_ub;
        public:
          ViewBasketAsUB( const basket_t& b ) noexcept
          {
            detail::UBImpl::access_internal(m_ub)
              = const_cast<void*>(static_cast<const void*>(&b));
            assign_basket_parameters( m_ub, b );
          }
          ~ViewBasketAsUB()
          {
            detail::UBImpl::access_internal(m_ub) = nullptr;
          }
          const Basket& view() noexcept { return m_ub; }
        };


        static std::size_t append1( Basket& dst, const Basket& src,
                                    std::size_t idx_src )
        {
          return real_basket(dst).append1( real_basket(src), idx_src );
        }

        static heapmem_t releaseMem( Basket&& b )
        {
          nc_assert(ub_internal(b)!=nullptr);
          Basket tmp;
          tmp.swap(b);
          heapmem_t mem;
          ub_internal(tmp) = nullptr;
          mem.set_raw(ub_internal(tmp));
          return mem;
        }

        static void emergencyCleanupMemory( Basket&& b )
        {
          releaseMem( std::move(b) ).deallocate();
        }

        static Basket constructBasketFromHeapMem( heapmem_t mem )
        {
          nc_assert(mem.data()!=nullptr);
          //Construct a basket object in the heap memory and move ownership to
          //the result basket. It is important that the result basket is fully
          //initialised before any exception might be thrown:
          static_assert( std::is_nothrow_constructible<basket_t>::value, "" );
          basket_t * b  = new(mem.release_raw()) basket_t;
          Basket res;
          ub_internal(res) = (void*)b;
          assign_basket_parameters( res, *b );
          nc_assert( res.neutrons != nullptr && res.neutrons->nused == 0);
          return res;
        }

        class BskMgr final : public BasketMgr {
        public:
          BskMgr( GeometryPtr geom, SourcePtr src )
            : m_srcfiller( std::move(geom), std::move(src) )
          {
          }

          void haltSource() override
          {
#ifndef NCRYSTAL_DISABLE_THREADS
            std::unique_lock<std::mutex> lock(m_mutex);
#endif
            m_data.src_ran_out = true;
          }

          void haltError() override
          {
#ifndef NCRYSTAL_DISABLE_THREADS
            std::unique_lock<std::mutex> lock(m_mutex);
#endif
            m_data.src_ran_out = true;
            m_data.is_error = true;
#ifndef NCRYSTAL_DISABLE_THREADS
            m_condvar.notify_all();
#endif
          }

          Basket allocateBasketNoLock()
          {
            //call only when lock is already acquired
            return constructBasketFromHeapMem( m_data.mempool.allocate() );
          }

          Basket allocateBasket() override
          {
#ifndef NCRYSTAL_DISABLE_THREADS
            std::unique_lock<std::mutex> lock(m_mutex);
#endif
            return allocateBasketNoLock();
          }

          void deallocateBasketNoLock( Basket&& b )
          {
            //call only when lock is already acquired
            if ( b.valid() ) {
              auto mem = releaseMem(std::move(b));
              m_data.mempool.deallocate( std::move(mem) );
            }
          }

          void deallocateBasket( Basket&& b ) override
          {
            if ( b.valid() ) {
#ifndef NCRYSTAL_DISABLE_THREADS
              std::unique_lock<std::mutex> lock(m_mutex);
#endif
              deallocateBasketNoLock( std::move(b) );
            }
          }

          void addPendingBasket( Basket&& b ) override
          {
            nc_assert( b.valid() );
#ifndef NCRYSTAL_DISABLE_THREADS
            std::unique_lock<std::mutex> lock(m_mutex);
#endif
            if ( m_data.is_error || b.empty() ) {
              auto mem = releaseMem(std::move(b));
              m_data.mempool.deallocate( std::move(mem) );
              return;
            }
            m_data.pending.emplace_back( std::move(b) );
#ifndef NCRYSTAL_DISABLE_THREADS
            m_condvar.notify_one();
#endif
          }

          Basket getInputBasket( WorkerToken& token,
                                 RNG& rng,
                                 const TallyFct& tallyfct,
                                 ParticleCountSum& missCount ) override
          {
            Basket result;
            std::function<void(const basket_t&)> resultfct;
            if ( tallyfct != nullptr ) {
              resultfct = [&tallyfct]( const basket_t& b )
              {
                tallyfct( ViewBasketAsUB( b ).view() );
              };
            }
#ifdef NCRYSTAL_DISABLE_THREADS
            if ( !m_data.pending.empty() ) {
              result = std::move( m_data.pending.back() );
              m_data.pending.pop_back();
            } else {
              result = allocateBasketNoLock();
            }
            nc_assert(result.valid());
            if ( result.size() > basket_N_almost_Full )
              return result;
            if ( !m_data.src_ran_out ) {
              Basket& tk_buf
                = detail::UBImpl::token_basket_buf_refresh(token,*this);
              nc_assert( tk_buf.valid() && tk_buf.size() == 0 );
              if (!m_srcfiller.fillFromSource( real_basket( result ),
                                               rng,
                                               resultfct,
                                               missCount,
                                               real_basket(tk_buf))) {
                m_data.src_ran_out = true;
              }
            }
            if ( result.size() == 0 )
              deallocateBasketNoLock( std::move(result) );//done!
            return result;
#else
            bool use_source = false;
            bool& tk_active = detail::UBImpl::access_token_activeflag( token );
            std::unique_lock<std::mutex> lock(m_mutex);
            if ( !tk_active ) {
              tk_active = true;
              ++m_data.ntokens;
            }
            ++m_data.ntokenswaiting;
            auto& data = m_data;
            m_condvar.wait( lock, [&data, &use_source, &result ]()
            {
              if ( data.is_error ) {
                use_source = false;
                return true;//done waiting: general error
              }
              use_source = !data.src_ran_out;
              if ( !data.pending.empty() ) {
                result = std::move( data.pending.back() );
                data.pending.pop_back();
                if ( result.size() > basket_N_almost_Full )
                  use_source = false;
                return true;//done waiting: got something
              }
              //Nothing pending, check the source:
              if ( use_source )
                return true;//done waiting: must fill new basket from source!

              //Nothing pending, source ran out, and there was no error. If all
              //worker threads are waiting for something, this is a sign that
              //the simulation is actually complete.
              if ( data.ntokens > data.ntokenswaiting ) {
                return false;//keep waiting: some workers are still working and
                             //might produce new pending baskets.
              }
              return true;//done waiting: no more work to do
            });

            nc_assert( bool(lock) );//we still have the lock
            --m_data.ntokenswaiting;

            if ( use_source && !result.valid() ) {
              //allocate a new basket while we still have the lock
              result = constructBasketFromHeapMem( m_data.mempool.allocate() );
              nc_assert(result.size()==0);
            }

            lock.unlock();//release lock, src filling can happen concurrently

            if ( use_source ) {
              nc_assert( result.valid() );
              if ( result.size() < basket_N_almost_Full ) {
                //Ensure we have an extra basket buffer, needed for the
                //fillFromSource call:
                Basket& tk_buf
                  = detail::UBImpl::token_basket_buf_refresh(token,*this);
                nc_assert( tk_buf.valid() && tk_buf.size() == 0 );
                if (!m_srcfiller.fillFromSource( real_basket( result ),
                                                 rng,
                                                 resultfct,
                                                 missCount,
                                                 real_basket(tk_buf))) {
                  //source ran out, must lock again to note this down:
                  lock.lock();
                  m_data.src_ran_out = true;
                  lock.unlock();
                }
              }
              nc_assert(!bool(lock));
              if ( result.valid() && result.empty() ) {
                //Unlikely to happen, but just in case:
                deallocateBasket( std::move(result) );
              }
            }
            return result;
#endif
          }

          void tokenDispose( WorkerToken& tk ) override
          {
            bool& tk_active = detail::UBImpl::access_token_activeflag( tk );
            Basket tk_basket = std::move( detail::UBImpl
                                          ::access_token_basket_buf(tk) );
            if ( !tk_active && !tk_basket.valid() )
              return;
#ifndef NCRYSTAL_DISABLE_THREADS
            std::unique_lock<std::mutex> lock(m_mutex);//NB: not noexcept!
#endif
            if ( tk_basket.valid() )
              deallocateBasketNoLock( std::move(tk_basket) );
            if ( tk_active ) {
              tk_active = false;
#ifndef NCRYSTAL_DISABLE_THREADS
              nc_assert( m_data.ntokens > 0 );
              --m_data.ntokens;
              m_condvar.notify_all();
#endif
            }
          }

        private:
#ifndef NCRYSTAL_DISABLE_THREADS
          std::mutex m_mutex;
          std::condition_variable m_condvar;
#endif
          //The condvar and data on the following struct is protected by the
          //mutex access:
          struct Data {
#ifndef NCRYSTAL_DISABLE_THREADS
            std::size_t ntokens = 0;
            std::size_t ntokenswaiting = 0;
#endif
            std::vector<Basket> pending;
            heapmempool_t mempool;
            bool src_ran_out = false;
            bool is_error = false;
          };
          Data m_data;
          //The srcfiller can be used concurrently, without locking m_mutex:
          basket_src_filler_t m_srcfiller;
        };
      };
    }
  }
}

std::size_t
NCMMC::Basket::append1( const Basket& other, std::size_t idx_other )
{
  nc_assert(valid()&&other.valid());
  nc_assert(basketType() == other.basketType());
  nc_assert( size() < basket_N );
  if ( neutrons_initial ) {
    nc_assert(basketType()==BasketType::Extended);
    return UBHImpl<Basket_Extended>::append1( *this, other, idx_other );
  } else {
    nc_assert(basketType()==BasketType::Basic);
    return UBHImpl<Basket_Basic>::append1( *this, other, idx_other );
  }
}

void NCMMC::Basket::dealloc_warn() noexcept
{
  //Called from destructor in the rare case the basket was not returned
  //explicitly to a memory pool before it went out of scope. This could happen
  //in case of an exception elsewhere, so we want to prevent a memory leak in
  //this case but we emit warning so we can catch regular coding errors in the
  //simulation engine.
  //
  //We can NOT leak from this as we are called from a destructor.
  try {
    {
#ifndef NDEBUG
      //only in debug builds, this is for developers not when user callbacks
      //throws an exception.
      NCRYSTAL_WARN("MiniMC Basket went out of scope without"
                    " being handed to the manager.");
#endif
      nc_assert_always(internal);
      auto bt = this->basketType();
      Basket tmp;
      tmp.swap(*this);
      if ( bt == BasketType::Basic ) {
        UBHImpl<Basket_Basic>::emergencyCleanupMemory(std::move(tmp));
      } else if ( bt == BasketType::Extended ) {
        UBHImpl<Basket_Extended>::emergencyCleanupMemory(std::move(tmp));
      } else {
        nc_assert_always(false);
      }
    }
  } catch ( std::exception& e ) {
    //try something without exceptions:
    const char * msg = e.what();
    std::printf("NCrystal ERROR:: unexpected exception during"
                " Basket emergency cleanup: %s\n",
                (msg?msg:"<unknown>"));
  }
}

NC::shared_obj<NCMMC::BasketMgr>
NCMMC::createBasketMgr( GeometryPtr geom,
                        SourcePtr src,
                        BasketType bt )
{
  if ( bt == BasketType::Basic ) {
    return makeSO<UBHImpl<Basket_Basic>::BskMgr>
      ( std::move(geom), std::move(src) );
  } else {
    if ( bt != BasketType::Extended )
      NCRYSTAL_THROW(BadInput,"Invalid basket type.");
    return makeSO<UBHImpl<Basket_Extended>::BskMgr>
      ( std::move(geom), std::move(src) );
  }
}

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

#include "NCrystal/internal/minimc/NCMMC_UBView.hh"
#include "NCMMC_BasketSrcFiller.hh"
#include "NCMMC_BasketUtils.hh"

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    class detail::UBImpl {
    public:
      static void*& access_internal( UniversalBasket& b ) noexcept
      {
        return b.internal;
      }
      static const void* access_internalc( const UniversalBasket& b ) noexcept
      {
        return b.internal;
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

        void assignToUB( UniversalBasket& dst ) noexcept
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

        void assignToUB( UniversalBasket& dst ) noexcept
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

      inline void *& ub_internal( UniversalBasket& b )
      {
        return detail::UBImpl::access_internal(b);
      }

      template<class TBasket>
      struct UBHImpl{
        //typedef templated internal data structures:
        using basket_t = TBasket;
        using basket_src_filler_t = BasketSrcFiller<basket_t>;
        using basketmgr_t = BasketMgr<basket_t>;
        using basket_holder_t = typename basketmgr_t::basket_holder_t;
        static constexpr int local_basket_max_poolsize = 4;
        using heapmem_t = typename basket_holder_t::heapmem_t;
        using heapmempool_t = HeapMemPool<heapmem_t,local_basket_max_poolsize>;

        static basket_t& real_basket( UniversalBasket& b )
        {
          auto rbptr = reinterpret_cast<basket_t*>
            ( detail::UBImpl::access_internal(b) );
          nc_assert( rbptr != nullptr );
          return *rbptr;
        }

        static const basket_t& real_basket( const UniversalBasket& b )
        {
          auto rbptr = reinterpret_cast<const basket_t*>
            ( detail::UBImpl::access_internalc(b) );
          nc_assert( rbptr != nullptr );
          return *rbptr;
        }

        static void assign_basket_parameters( UniversalBasket& dst,
                                              const basket_t& src_c ) noexcept
        {
          auto& src = *const_cast<basket_t*>(&src_c);
          src.assignToUB( dst );
        }

        static UniversalBasket toUniversalBasket( basket_holder_t&& bh ) {
          //NB: Map invalid basket to invalid basket.
          UniversalBasket res;
          nc_assert( !res.valid() );
          if ( bh.valid() ) {
            assign_basket_parameters( res, bh.basket() );
            ub_internal(res) = bh.stealMemory().release_raw();
            nc_assert( res.valid() );
          }
          return res;
        }

        class ViewBasketAsUB final {
          //Non-owning way to view a basket_t as a UniversalBasket. The main point
          //is to reset the internal pointer of the UniversalBasket object in the
          //destructor, before it is itself destructed.
          UniversalBasket m_ub;
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
          const UniversalBasket& view() noexcept { return m_ub; }
        };


        static std::size_t append1( UniversalBasket& dst,
                                    const UniversalBasket& src,
                                    std::size_t idx_src )
        {
          return real_basket(dst).append1( real_basket(src), idx_src );
        }

        static basket_holder_t moveToBasketHolder( UniversalBasket&& b )
        {
          nc_assert_always(ub_internal(b)!=nullptr);
          UniversalBasket tmp;
          tmp.swap(b);
          heapmem_t mem;
          mem.set_raw(ub_internal(tmp));
          ub_internal(tmp) = nullptr;
          basket_holder_t bh{ no_init };
          bh.set_raw_from_active_heapmem( std::move(mem) );
          return bh;
        }

        class UBMgr final : public UniversalBasketMgr {
        public:
          //NOTE: basketmgr_t is MT safe, but UBMgr is NOT, needing one cloned
          //instance per thread.
          UBMgr( std::shared_ptr<basketmgr_t> bmgr = nullptr )
            : m_mgr( [&bmgr]()
            { return ( bmgr==nullptr
                       ? std::make_shared<basketmgr_t>()
                       : std::move(bmgr) ); }() )
          {
          }

          UniversalBasket allocateBasket() override
          {
            UniversalBasket res;
            //Get the heap memory:
            heapmem_t mem;
            if ( !m_localMemPool.empty() )
              mem = m_localMemPool.allocate();
            else
              mem = m_mgr->allocateBasket().stealMemory();
            nc_assert(mem.data()!=nullptr);
            //Construct a basket object in the heap memory and move ownership to
            //the result basket. It is important that the result basket is fully
            //initialised before any exception might be thrown:
            static_assert( std::is_nothrow_constructible<basket_t>::value, "" );
            basket_t * b  = new(mem.release_raw()) basket_t;
            ub_internal(res) = (void*)b;
            assign_basket_parameters( res, *b );
            return res;
          }

          void deallocateBasket( UniversalBasket&& b ) override
          {
            nc_assert( b.valid() );
            auto mem = moveToBasketHolder( std::move(b) ).stealMemory();
            nc_assert( mem.data() != nullptr );
            //return mem to local mempool or overflow to global mgr:
            nc_assert( m_localMemPool.size() <= local_basket_max_poolsize );
            if ( m_localMemPool.size() == local_basket_max_poolsize )
              m_mgr->deallocateHeapMem( std::move(mem) );
            else
              m_localMemPool.deallocate( std::move(mem) );
          }

          void addPendingBasket( UniversalBasket&& b ) override
          {
            nc_assert( b.valid() );
            if ( b.empty() )
              deallocateBasket( std::move(b) );
            else
              m_mgr->addPendingBasket( moveToBasketHolder( std::move(b) ) );
          }

          shared_obj<UniversalBasketMgr> cloneMgrForThread() override
          {
            return makeSO<UBMgr>( m_mgr );
          }

          shared_obj<basketmgr_t> realManager() { return m_mgr; }

        private:
          shared_obj<basketmgr_t> m_mgr;//the actual manager implementation, one
                                        //for all threads.
          heapmempool_t m_localMemPool;//A small thread-local memory pool.
        };

        class InBskProv final : public InputBasketProvider {
          basket_src_filler_t m_srcfiller;
          ThreadCount m_nthreads;

          static ThreadedUsage determineThreadedUsage(ThreadCount nthreads)
          {
#ifndef NCRYSTAL_DISABLE_THREADS
            return ( nthreads.get() > 1
                     ? ThreadedUsage::Multi
                     : ThreadedUsage::Single );
#else
            (void)nthreads;
            return ThreadedUsage::Single;
#endif
          }

        public:
          InBskProv( GeometryPtr geom,
                     SourcePtr src,
                     shared_obj<basketmgr_t> bmgr,
                     ThreadCount nthreads )
            : m_srcfiller( std::move(geom),
                           std::move(src),
                           std::move(bmgr),
                           determineThreadedUsage(nthreads) ),
              m_nthreads( nthreads )
          {
          }

          void haltSource() override
          {
            m_srcfiller.haltSource();
          }

          UniversalBasket getInputBasket( RNG& rng,
                                          const TallyFct& tallyfct,
                                          ParticleCountSum& missCount ) override
          {
            std::function<void(const basket_t&)> resultfct;
            if ( tallyfct != nullptr ) {
              resultfct = [&tallyfct]( const basket_t& b )
              {
                tallyfct( ViewBasketAsUB( b ).view() );
              };
            }
            auto bh = m_srcfiller.getPendingBasket( m_nthreads, rng,
                                                    resultfct, missCount );
            return toUniversalBasket( std::move(bh) );
          }
        };
      };
    }
  }
}

std::size_t
NCMMC::UniversalBasket::append1( const UniversalBasket& other,
                                 std::size_t idx_other )
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

void NCMMC::UniversalBasket::dealloc_warn() noexcept
{
  //Called from destructor in the rare case the basket was not returned
  //explicitly to a memory pool before it went out of scope. This could happen
  //in case of an exception elsewhere, so we want to prevent a memory leak in
  //this case but we emit warning so we can catch regular coding errors in the
  //simulation engine.
  //
  //We can NOT leak from this as we are called from a destructor.
  //
  //Fixme: try to verify that the logic below actually works, doesn't trigger
  //new issues, etc. At least do it once manually.

  try {
    {
      NCRYSTAL_WARN("MiniMC UniversalBasket went out of scope without"
                    " being handed to the manager.");
      nc_assert_always(!internal);
      auto bt = this->basketType();
      UniversalBasket tmp;
      tmp.swap(*this);
      if ( bt == BasketType::Basic ) {
        UBHImpl<Basket_Basic>::moveToBasketHolder(std::move(tmp));
      } else if ( bt == BasketType::Extended ) {
        UBHImpl<Basket_Extended>::moveToBasketHolder(std::move(tmp));
      } else {
        nc_assert_always(false);
      }
    }
  } catch ( std::exception& e ) {
    //try something without exceptions:
    const char * msg = e.what();
    std::printf("NCrystal ERROR:: unexpected exception during"
                " UniversalBasket emergency cleanup: %s\n",
                (msg?msg:"<unknown>"));
  }
}

NCMMC::BasketManagementPair
NCMMC::createBasketManagement( GeometryPtr geom,
                               SourcePtr src,
                               BasketType bt,
                               ThreadCount nthreads )
{
  if ( bt == BasketType::Basic ) {
    auto ubmgr = makeSO< UBHImpl<Basket_Basic>::UBMgr >();
    auto prov
      = makeSO<UBHImpl<Basket_Basic>::InBskProv>
      ( std::move(geom), std::move(src), ubmgr->realManager(), nthreads );
    return { std::move(ubmgr), std::move(prov) };
  } else {
    if ( bt != BasketType::Extended )
      NCRYSTAL_THROW(BadInput,"Invalid basket type.");
    auto ubmgr = makeSO< UBHImpl<Basket_Extended>::UBMgr >();
    auto prov
      = makeSO<UBHImpl<Basket_Extended>::InBskProv>
      ( std::move(geom), std::move(src), ubmgr->realManager(), nthreads );
    return { std::move(ubmgr), std::move(prov) };
  }
}

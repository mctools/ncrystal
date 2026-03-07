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

#include "NCrystal/internal/minimc/NCMMC_CBMgr.hh"
#include "NCrystal/internal/minimc/NCMMC_Baskets.hh"
#include "NCMMC_BasketUtils.hh"
#include "NCMMC_ParseCfg.hh"

namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

NCMMC::CB::DataArea::DataArea( std::size_t capacity )
  : m_capacity(capacity)
{
  //fixme: unit test this next exception:
  for ( auto i : ncrange(nfieldsmax) )
    m_datacache[i] = nullptr;
  if ( capacity > CBMgrInput::cachelen_max )
    NCRYSTAL_THROW2(BadInput,"Due to memory concerns, callback functions"
                    " are limited to providing at most "
                    <<CBMgrInput::cachelen_max<<" neutrons at a time");
}

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace CB {

      class DataArea::Mutable final : NoCopyMove {
      public:
        Mutable() = delete;
        static inline std::size_t& size(DataArea&da) noexcept
        {
          return da.m_size;
        }

        static inline decltype(DataArea::m_datacache)&
        datacache(DataArea&da) noexcept
        {
          return da.m_datacache;
        }

        static inline std::unique_ptr<double[]>& memholder(DataArea&da) noexcept
        {
          return da.m_memholder;
        }

        static inline BasketType& basketType(DataArea&da) noexcept
        {
          return da.m_basketType;
        }
      };

      namespace {

        void data_append_neutronfields( double** dst,
                                        unsigned offset,
                                        std::size_t this_size,
                                        const NeutronBasketFields& b,
                                        std::size_t b_size )
        {
          //order: x,y,z,ux,uy,uz,ekin,w
          if (!b_size)
            return;
          auto aa = [&dst,&offset,this_size,b_size]( const BasketValBufDbl& ov )
          {
            double * arrB = dst[offset++];
            BasketUtils::memcpydata<double>( arrB + this_size,
                                             ov.data, b_size );
          };
          aa(b.x);
          aa(b.y);
          aa(b.z);
          aa(b.ux);
          aa(b.uy);
          aa(b.uz);
          aa(b.ekin);
          aa(b.w);
        }

        void data_init( DataArea& da, const Basket& b )
        {
          auto& da_basketType = DataArea::Mutable::basketType(da);
          auto& da_datacache = DataArea::Mutable::datacache(da);
          auto& da_memholder = DataArea::Mutable::memholder(da);

          //Initialise based on the first basketview we see:
          BasketType bt = b.basketType();
          nc_assert( bt != BasketType::Invalid );
          nc_assert( bt == BasketType::Basic || bt == BasketType::Extended );
          nc_assert( b.nscat != nullptr );
          nc_assert( b.nscat_inelas != nullptr );
          nc_assert( (b.nscat_inelas == nullptr) == ( bt == BasketType::Basic) );
          const std::size_t nfields = ( bt == BasketType::Basic ? 10 : 18 );

          //Based on this, figure out number of fields and setup memory:
          nc_assert_always(da.size() == 0);
          nc_assert_always( da_basketType == BasketType::Invalid);//fixme _always
          nc_assert_always( nfields > 0 && nfields <= DataArea::nfieldsmax );
          nc_assert_always( da.capacity() > 0
                            && da.capacity() <= CBMgrInput::cachelen_max );
          const std::size_t nvalues = nfields * da.capacity();
          nc_assert_always( da_memholder == nullptr );
          da_memholder = ncmake_unique_array<double>( nvalues );
          double * it = da_memholder.get();
          for ( auto i : ncrange( nfields ) ) {
            da_datacache[i] = it;
            it += da.capacity();
          }
          da_basketType = bt;
          nc_assert_always( nfields == da.nFields() );
        }

        void data_append( DataArea& da, const Basket& b )
        {
          nc_assert(b.valid()&&b.neutrons);
          const std::size_t this_size = da.size();
          const std::size_t o_size = b.size();
          nc_assert_always( o_size <= basket_N );//fixme _always
          nc_assert_always( this_size + o_size <= da.capacity() );

          const auto& da_basketType = DataArea::Mutable::basketType(da);

          if ( da_basketType == BasketType::Invalid )
            data_init( da, b );
          nc_assert_always( da_basketType == BasketType::Basic
                            || da_basketType == BasketType::Extended );

          auto& da_datacache = DataArea::Mutable::datacache(da);
          data_append_neutronfields( da_datacache, 0, da.size(),
                                     b.neutrons->fields, b.size() );

          //Convert nscat and nscat_inelas to doubles:
          {
            auto b_nscat = b.nscat;
            nc_assert_always( b_nscat != nullptr );
            double * dst = da_datacache[8] + this_size;
            double * dstE = dst + o_size;
            const int * src = b_nscat->data;
            for ( ; dst!=dstE; ++dst, ++src )
              *dst = double(*src);
          }

          {
            auto b_nscat_inelas = b.nscat_inelas;
            nc_assert_always( b_nscat_inelas != nullptr );
            double * dst = da_datacache[9] + this_size;
            double * dstE = dst + o_size;
            const int * src = b_nscat_inelas->data;
            for ( ; dst != dstE; ++dst, ++src )
              *dst = double(*src);
          }

          if ( da_basketType == BasketType::Extended ) {
            auto nb_orig = b.neutrons_initial;
            nc_assert_always( nb_orig != nullptr );
            data_append_neutronfields( da_datacache, 10, da.size(),
                                       *nb_orig, b.size() );
          }

          DataArea::Mutable::size(da) += o_size;
        }

#ifndef NCRYSTAL_DISABLE_THREADS
        void data_append_other( DataArea& dst, DataArea& o )
        {
          nc_assert_always( &dst != &o );
          if ( dst.size() >= dst.capacity() || o.size() == 0 )
            return;//we are full, or the other area is empty.

          const std::size_t freecap = dst.capacity() - dst.size();
          const std::size_t ntransfer = ( o.size() <= freecap
                                          ? o.size() : freecap );
          nc_assert_always( ntransfer >= 1 );
          nc_assert_always( ntransfer <= o.size() );
          nc_assert_always( ntransfer <= freecap );
          nc_assert_always( dst.size()+ntransfer <= dst.capacity() );
          const std::size_t new_osize = o.size() - ntransfer;

          //Carry out the transfer:
          auto dst_data = DataArea::Mutable::datacache(dst);
          auto o_data = o.view_data();
          const std::size_t nfields = dst.nFields();
          nc_assert_always( nfields == o.nFields() );
          for ( auto ifield : ncrange(nfields) )
            BasketUtils::memcpydata<double>( dst_data[ifield] + dst.size(),
                                             o_data[ifield] + new_osize,
                                             ntransfer );
          DataArea::Mutable::size(o) = new_osize;
          DataArea::Mutable::size(dst) += ntransfer;
          nc_assert_always(dst.size()<=dst.capacity());
          nc_assert_always(o.size()<=o.capacity());
        }
#endif
      }
    }
  }
}

unsigned NCMMC::CB::DataArea::nFields() const
{
  if ( m_basketType == BasketType::Basic )
    return 10;
  nc_assert_always( m_basketType == BasketType::Extended );
  return 18;
}

void NCMMC::CB::CBMgr::registerData( const Basket& b )
{
#ifndef NCRYSTAL_DISABLE_THREADS
  auto datauptr = threadAcquireCache();
  nc_assert_always(datauptr!=nullptr);
  data_append( *datauptr, b );
  threadReturnCache( std::move(datauptr) );
#else
  if ( m_cache == nullptr ) {
    m_cache = ncmake_unique<DataArea>(m_nmax);
    nc_assert(!(m_cache->size() + b->size() > m_nmax));
  } else if ( m_cache->size() + b->size() > m_nmax ) {
    flush();
  }
  nc_assert_always(m_cache!=nullptr);
  data_append( *m_cache, b );
#endif
}

void NCMMC::CB::CBMgr::flush()
{
#ifndef NCRYSTAL_DISABLE_THREADS
  NCRYSTAL_LOCK_GUARD(m_cachemtx);//likely not needed
  //Merge data before firing it off, to potentially avoid multiple callbacks
  //with a small number of neutrons in each. Since we can not iterate over a
  //std::queue, we first move everything to a vector:
  if ( m_caches.empty() )
    return;

  SmallVector<DataAreaPtr,64> cachevector;
  while( !m_caches.empty() ) {
    cachevector.emplace_back(std::move(m_caches.front()));
    m_caches.pop();
    nc_assert_always(cachevector.back()!=nullptr);
  }

  const std::size_t capacity1 = cachevector.front()->capacity();
  std::size_t capacity_total = capacity1 * cachevector.size();
  std::size_t size_total = 0;
  for ( auto& e : cachevector )
    size_total += e->size();

  auto merge_makes_sense = [&cachevector,&capacity_total,size_total,capacity1]()
  {
    //Can we remove a basket and have room for all data?
    return cachevector.size() >= 2 && size_total+capacity1 <= capacity_total;
  };

  //We have to move less neutrons if we sort the cachevector by size, shortest
  //at the end:
  auto sort_cachevector = [&cachevector]()
  {
    std::stable_sort( cachevector.begin(), cachevector.end(),
                      [](const DataAreaPtr& a, const DataAreaPtr&b)
                      {
                        nc_assert_always(a!=nullptr&&b!=nullptr);
                        return a->size() > b->size();
                      } );
  };

  while ( merge_makes_sense() ) {
    nc_assert_always(cachevector.size()>=2);
    sort_cachevector();
    auto itBack = std::prev(cachevector.end());
    nc_assert_always(itBack->get()!=nullptr);
    nc_assert_always(std::prev(itBack)->get()!=nullptr);
    nc_assert_always(std::prev(itBack)>=cachevector.begin());
    DataArea& src_area = *(itBack->get());
    DataArea& tgt_area = *(std::prev(itBack)->get());
    nc_assert_always(tgt_area.size()<tgt_area.capacity());//should not be able
                                                          //to happen if
                                                          //merge_makes_sense()
    if ( src_area.size() > 0 )
      data_append_other(tgt_area, src_area );
    if ( src_area.size() == 0 ) {
      cachevector.pop_back();//fully absorbed, discard
      nc_assert_always(capacity_total>capacity1);
      capacity_total -= capacity1;
      nc_assert_always(capacity_total>size_total);
    }
  }

  //Fire the callbacks, most to least data (just seems more in line with what
  //users might expect: potentially one not very full area at the very end).
  sort_cachevector();
  std::reverse( cachevector.begin(), cachevector.end() );
  while( !cachevector.empty() ) {
    DataAreaPtr data = std::move(cachevector.back());
    nc_assert_always(data!=nullptr);
    cachevector.pop_back();
    if ( data != nullptr && data->size() > 0 )
      fireCallback( *data );
  }
#else
  if ( m_cache!=nullptr && m_cache->size() > 0 ) {
    nc_assert( m_callback != nullptr );
    m_callback( *m_cache );
    DataArea::Mutable::size(*m_cache) = 0;
    //m_cache->clear();
  }
#endif
}

NCMMC::CB::CBMgr::CBMgr( CBMgrInput input )
  : m_nmax( input.cachelen > basket_N ? input.cachelen : basket_N ),
#ifndef NCRYSTAL_DISABLE_THREADS
    m_nmax_caches( input.ncaches > 1 ? input.ncaches : 1 ),
#endif
    m_callback(std::move(input.callbackfct))
{
  nc_assert_always(m_callback!=nullptr);

  //1e9 is high, memorable, fits in any kind of integer of at least 32 bits.
  if ( input.cachelen > CBMgrInput::cachelen_max )
    NCRYSTAL_THROW2(BadInput,"cachelen can never exceed "
                    <<CBMgrInput::cachelen_max)

  if ( input.cachelen < m_nmax ) {
    NCRYSTAL_WARN("cachelen was increased to"
                  " the minimum value of "<<m_nmax);
  } else if ( input.cachelen
              > static_cast<std::size_t>(5e7) ) {
    NCRYSTAL_WARN("Value of cachelen is extremely"
                  " high. This might lead to performance issues or"
                  " crashes.");
  }
}

NCMMC::CB::CBMgr::~CBMgr()
{
  //flush() should have been invoked before destruction (NB: no exceptions from
  //destructor!)
#ifndef NCRYSTAL_DISABLE_THREADS
  assert( m_caches.empty() );
#else
  assert( m_cache == nullptr || m_cache->size() == 0 );
#endif
}


#ifndef NCRYSTAL_DISABLE_THREADS

NCMMC::CB::CBMgr::DataAreaPtr NCMMC::CB::CBMgr::threadAcquireCache()
{
  std::unique_lock<std::mutex> lock(m_cachemtx);
  DataAreaPtr res;
  m_cachecondvar.wait( lock,[&res,this]()
  {
    if ( !m_caches.empty() ) {
      res = std::move(m_caches.front());
      m_caches.pop();
    }
    if ( res == nullptr && m_ncaches_created < m_nmax_caches ) {
      res = ncmake_unique<DataArea>(m_nmax);
      ++m_ncaches_created;
    }
    return res != nullptr;
  });
  nc_assert_always( res != nullptr );
  return res;
}

void NCMMC::CB::CBMgr::threadReturnCache( DataAreaPtr data )
{
  nc_assert_always(data!=nullptr);
  if ( data->size()+basket_N > m_nmax ) {
    //Not room for another basket. Deliver + clear:
    fireCallback( *data );
    DataArea::Mutable::size(*data) = 0;
    nc_assert_always( data->size() == 0 );
  }
  //Return to queue:
  std::unique_lock<std::mutex> lock(m_cachemtx);
  m_caches.push( std::move(data) );
  m_cachecondvar.notify_one();
}

void NCMMC::CB::CBMgr::fireCallback( const DataArea& data )
{
  nc_assert_always(data.size()>0);
  NCRYSTAL_LOCK_GUARD( m_cbmtx );
  nc_assert( m_callback != nullptr );
  m_callback( data );
}

NCMMC::CB::CBMgrInput NCMMC::CB::decodeCBMgrInput( StrView raw_str )
{
  CBMgrInput res;
  namespace PMC = parseMMCCfg;
  auto tokeninfo = PMC::tokenize( raw_str );
  if ( tokeninfo.mainName.has_value() )
    NCRYSTAL_THROW2(BadInput,"Invalid callback cfg: \""<<raw_str<<"\"");
  auto& tokens = tokeninfo.tokens;
  PMC::checkNoUnknown(tokens,"basket_type;cachelen;ncaches","callback");

  {
    auto val = PMC::getValue_str_allowempty(tokens,"basket_type");
    if ( val=="basic" )
      res.basketType = BasketType::Basic;
    else if ( val=="extended" )
      res.basketType = BasketType::Extended;
    else if ( !val.empty() )
      NCRYSTAL_THROW2(BadInput,"Invalid basket_type \""<<val
                      <<"\" (must be \"basic\" or \"extended\").");
  }
  if ( PMC::hasValue(tokens,"cachelen") ) {
    auto val = PMC::getValue_sizet(tokens,"cachelen");
    if ( val > CBMgrInput::cachelen_max )
      NCRYSTAL_THROW2(BadInput,"Invalid cachelen \""<<val
                      <<"\" (must be less than "
                      <<CBMgrInput::cachelen_max<<")");
    if ( val < basket_N )
      val = basket_N;
    res.cachelen = val;
  }
  if ( PMC::hasValue(tokens,"ncaches") ) {
    auto val = PMC::getValue_sizet(tokens,"ncaches");
    res.ncaches = ( val>=1 ? val : 1 );
  }

  return res;
}

#endif

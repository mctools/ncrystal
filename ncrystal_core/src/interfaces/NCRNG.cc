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

#include "NCrystal/interfaces/NCRNG.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"

#ifndef NCRYSTAL_DISABLE_THREADS
#  include <thread>
namespace NCRYSTAL_NAMESPACE { namespace { using ThreadID = std::thread::id; } }
#else
namespace NCRYSTAL_NAMESPACE { namespace { using ThreadID = int; } }
#endif


namespace NC = NCrystal;

bool NC::RNGStream::coinflip()
{
  return generate() > 0.5;
}

uint64_t NC::RNGStream::generate64RndmBits()
{
  //As generate() returns a double with only 53 random bits, we need two
  //calls. Most straight forward to do this by combining two 32 bit integers:
  const uint32_t g1 = generate32RndmBits();
  const uint32_t g2 = generate32RndmBits();
  return ( uint64_t{ g1 } << 32 ) | g2;
}

uint32_t NC::RNGStream::generate32RndmBits()
{
  return static_cast<uint32_t>(generate()*std::numeric_limits<uint32_t>::max());
}

uint32_t NC::RNGStream::stateTypeUID() const noexcept
{
  return 0;
}

void NC::RNGStream::actualSetState( std::vector<uint8_t>&& )
{
}

std::vector<uint8_t> NC::RNGStream::actualGetState() const
{
  return {};
}

NC::shared_obj<NC::RNGStream> NC::RNGStream::actualCloneWithNewState( std::vector<uint8_t>&& ) const
{
  return std::shared_ptr<RNGStream>{nullptr};//dummy (will trigger exception immediately)
}

NC::RNGStreamState NC::RNGStream::getState() const
{
  uint32_t stateuid = stateTypeUID();
  if (!stateuid)
    NCRYSTAL_THROW(LogicError,"RNG::getState should never be called without first checking supportsStateManipulation().");
  std::vector<uint8_t> v = actualGetState();
  nc_assert_always(!v.empty());
  v.reserve(v.size()+sizeof(stateuid));
  appendToStateVector(v, stateuid);
  return RNGStreamState{ bytes2hexstr(v) };
}

namespace NCRYSTAL_NAMESPACE {
  namespace RNGStream_detail {

    constexpr uint32_t builtinRNGStateTypeUID = 0xb067bd44;//randomly generated

    uint32_t extractStateUID( const char * fullfct, const std::string& state )
    {
      std::vector<uint8_t> v = hexstr2bytes(state);
      if ( ! ( v.size() > sizeof(uint32_t) ) )
        NCRYSTAL_THROW2(BadInput,fullfct<<" got too short state.");
      return RNGStream::popFromStateVector<uint32_t>(v);
    }

    std::vector<uint8_t> extractStateBytes( const char * fct, const std::string& state, uint32_t expected_stateuid )
    {
      if (!expected_stateuid)
        NCRYSTAL_THROW2(LogicError,"RNGStream::"<<fct<<" should never be called without first checking supportsStateManipulation().");

      std::vector<uint8_t> v = hexstr2bytes(state);
      if ( ! ( v.size() > sizeof(expected_stateuid) ) )
        NCRYSTAL_THROW2(BadInput,"RNGStream::"<<fct<<" got too short state.");

      if ( RNGStream::popFromStateVector<uint32_t>(v) != expected_stateuid )
        NCRYSTAL_THROW2(BadInput,"RNGStream::"<<fct<<" got invalid state (or state originating in different RNG implementation).");

      return v;
    }
  }
}

void NC::RNGStream::setState( const RNGStreamState& state)
{
  actualSetState( RNGStream_detail::extractStateBytes("setState",state.get(),stateTypeUID()) );
  nc_assert( getState() == state );
}

NC::shared_obj<NC::RNGStream> NC::RNGStream::cloneWithNewState( const RNGStreamState& state ) const
{
  auto rng = actualCloneWithNewState( RNGStream_detail::extractStateBytes("cloneWithNewState",state.get(),stateTypeUID()) );
  nc_assert( rng->getState() == state );
  return rng;
}

NC::shared_obj<NC::RNGStream> NC::RNGStream::createJumped() const
{
  NCRYSTAL_THROW(LogicError,"createJumped() is not supported by this RNG stream (check isJumpCapable() before calling).");
  return optional_shared_obj<RNGStream>{nullptr};//can't just return nullptr when return type is shared_obj
}

bool NC::stateIsFromBuiltinRNG( const RNGStreamState& state )
{
  return RNGStream_detail::extractStateUID( "NCrystal::stateIsFromBuiltinRNG", state.get() ) == RNGStream_detail::builtinRNGStateTypeUID;
}

namespace NCRYSTAL_NAMESPACE {
  class RNG_XRSR final : public RNGStream {
  public:

    RNG_XRSR( uint64_t seed ) : m_impl{seed} {}
    RNG_XRSR( RandXRSRImpl&& impl ) : m_impl(std::move(impl)) {}
    RNG_XRSR( no_init_t ) : m_impl(no_init) {}

    bool coinflip() override { return m_impl.coinflip(); }
    uint64_t generate64RndmBits() override { return m_impl.genUInt64(); }
    uint32_t generate32RndmBits() override { return m_impl.genUInt32(); }

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
    void generateMany( std::size_t n, double* tgt ) override
    {
      return m_impl.generateMany(n,tgt);
    }
#endif

  protected:

    double actualGenerate() override { return m_impl.generate(); }

    uint32_t stateTypeUID() const noexcept override {
      return RNGStream_detail::builtinRNGStateTypeUID;
    }

    static RandXRSRImpl::state_t detail_convstate(std::vector<uint8_t>&& v)
    {
      RandXRSRImpl::state_t newstate;
      nc_assert_always( v.size() == 2*sizeof(uint64_t) );
      newstate[1] = popFromStateVector<uint64_t>(v);
      nc_assert( v.size() == sizeof(uint64_t) );
      newstate[0] = popFromStateVector<uint64_t>(v);
      nc_assert(v.empty());
      return newstate;
    }
    void actualSetState( std::vector<uint8_t>&& v ) override
    {
      m_impl.state() = detail_convstate(std::move(v));
    }
    std::vector<uint8_t> actualGetState() const override
    {
      std::vector<uint8_t> v;
      v.reserve( 2*sizeof(uint64_t) );
      appendToStateVector<uint64_t>(v,m_impl.state()[0]);
      appendToStateVector<uint64_t>(v,m_impl.state()[1]);
      nc_assert( v.size() == 2*sizeof(uint64_t) );
      return v;
    }

    shared_obj<RNGStream> actualCloneWithNewState( std::vector<uint8_t>&& v ) const override
    {
      return makeSO<RNG_XRSR>( RandXRSRImpl( detail_convstate(std::move(v)) ) );
    }

    bool isJumpCapable() const override
    {
      return true;
    }

    shared_obj<RNGStream> createJumped() const override
    {
      auto clone = makeSO<RNG_XRSR>(RandXRSRImpl(m_impl.state()));
      clone->m_impl.jump();
      return clone;
    }

  private:
    RandXRSRImpl m_impl;
  };

  class RNG_OneFctForAllThreads final : public RNGStream {
  public:
    RNG_OneFctForAllThreads( std::function<double()> fct ) : m_fct{std::move(fct)} {}
    bool useInAllThreads() const override { return true; }
  protected:
    double actualGenerate() override { return m_fct(); }
  private:
    std::function<double()> m_fct;
  };
}

NC::shared_obj<NC::RNGStream> NC::createBuiltinRNG( uint64_t seed )
{
  return makeSO<RNG_XRSR>(seed);
}

NC::shared_obj<NC::RNGStream> NC::createBuiltinRNG( const RNGStreamState& state )
{
  auto rng = makeSO<RNG_XRSR>( no_init );
  rng->setState(state);
  return rng;
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    struct DefRNGProd {
      std::mutex mtx;
      optional_shared_obj<RNGProducer> producer;
    };
    DefRNGProd& defRNGProdDB() {
      static DefRNGProd s_rngprod;
      return s_rngprod;
    }
  }
}

NC::shared_obj<NC::RNGProducer> NC::getDefaultRNGProducer()
{
  auto& d = defRNGProdDB();
  NCRYSTAL_LOCK_GUARD(d.mtx);
  if ( d.producer == nullptr )
    d.producer = makeSO<RNGProducer>( createBuiltinRNG() );
  return d.producer;
}

void NC::setDefaultRNG( NC::shared_obj<NC::RNGStream> rng )
{
  auto newprod = makeSO<RNGProducer>( std::move(rng) );
  auto& d = defRNGProdDB();
  NCRYSTAL_LOCK_GUARD(d.mtx);
  d.producer = newprod;
}

void NC::clearDefaultRNG()
{
  auto& d = defRNGProdDB();
  NCRYSTAL_LOCK_GUARD(d.mtx);
  d.producer = nullptr;
}

void NC::setDefaultRNGFctForAllThreads( std::function<double()> fct )
{
  setDefaultRNG(makeSO<RNG_OneFctForAllThreads>(fct));
}

namespace NCRYSTAL_NAMESPACE {
  struct RNGProducer::Impl final {
    Impl( shared_obj<RNGStream> rng ) : m_nextproduct( std::move(rng) ) {}
    Impl( no_init_t ) {}
    optional_shared_obj<RNGStream> m_nextproduct;
    optional_shared_obj<RNGStream> m_nextnextproduct;
    std::map<RNGStreamIndex,optional_shared_obj<RNGStream>> m_idxdb;
    std::map<ThreadID,optional_shared_obj<RNGStream>> m_thread_idxdb;
    std::mutex m_mtx;
    void jumpFillNextNextIfAppropriate();
    shared_obj<RNGStream> produceUnlocked();
    shared_obj<RNGStream> produceByIdxUnlocked( RNGStreamIndex );
    shared_obj<RNGStream> produceByThreadIdxUnlocked( ThreadID );
    static uint64_t currentThreadID();
  };
}
void NC::RNGProducer::Impl::jumpFillNextNextIfAppropriate()
{
  nc_assert_always(m_nextnextproduct==nullptr);
  if ( m_nextproduct == nullptr || m_nextproduct->useInAllThreads() || !m_nextproduct->isJumpCapable() )
    return;
  m_nextnextproduct = m_nextproduct->createJumped();
  if ( m_nextnextproduct == nullptr )
    NCRYSTAL_THROW(LogicError,"RNG stream claimed to be jump capable but a call to produce() returned nullptr!");
}

NC::shared_obj<NC::RNGStream> NC::RNGProducer::Impl::produceUnlocked()
{
  if ( !m_nextproduct )
    NCRYSTAL_THROW(CalcError,"Can not produce more independent RNG streams.");
  if ( m_nextproduct->useInAllThreads() )
    return m_nextproduct;//special case, just keep producing the same
  auto result = std::move(m_nextproduct);
  m_nextproduct = std::move(m_nextnextproduct);
  //Create jumped state immediately if possible (before anyone consumes
  //numbers from m_nextproduct, thereby modifying the state):
  jumpFillNextNextIfAppropriate();
  return result;
}

NC::shared_obj<NC::RNGStream> NC::RNGProducer::Impl::produceByIdxUnlocked( RNGStreamIndex idx )
{
  optional_shared_obj<RNGStream>& entry = m_idxdb[idx];
  if ( entry == nullptr )
    entry = produceUnlocked();
  return entry;
}

NC::shared_obj<NC::RNGStream> NC::RNGProducer::Impl::produceByThreadIdxUnlocked( ThreadID idx )
{
  optional_shared_obj<RNGStream>& entry = m_thread_idxdb[idx];
  if ( entry == nullptr )
    entry = produceUnlocked();
  return entry;
}

NC::shared_obj<NC::RNGStream> NC::RNGProducer::produce()
{
  NCRYSTAL_LOCK_GUARD(m_impl->m_mtx);
  return m_impl->produceUnlocked();
}

NC::shared_obj<NC::RNGStream> NC::RNGProducer::produceByIdx( RNGStreamIndex idx )
{
  NCRYSTAL_LOCK_GUARD(m_impl->m_mtx);
  return m_impl->produceByIdxUnlocked(idx);
}

NC::shared_obj<NC::RNGStream> NC::RNGProducer::produceForCurrentThread()
{
#ifndef NCRYSTAL_DISABLE_THREADS
  ThreadID thread_id = std::this_thread::get_id();
#else
  ThreadID thread_id = 1;//always the same!
#endif
  NCRYSTAL_LOCK_GUARD(m_impl->m_mtx);
  return m_impl->produceByThreadIdxUnlocked(thread_id);
}

NC::RNGProducer::RNGProducer( no_init_t )
  : m_impl( no_init )
{
  //sterile null producer
}

NC::shared_obj<NC::RNGProducer> NC::RNGProducer::getNullProducer()
{
  static shared_obj<RNGProducer> nullproducer = NC::makeSO<RNGProducer>(no_init);
  return nullproducer;
}

NC::RNGProducer::RNGProducer( shared_obj<RNGStream> rng, SkipOriginal skip_orig )
  : m_impl( std::move(rng) )
{
  //Create jumped state immediately if possible (before anyone consumes
  //numbers from m_nextproduct, thereby modifying the state):
  m_impl->jumpFillNextNextIfAppropriate();
  nc_assert_always( m_impl->m_nextproduct != nullptr );

  if ( skip_orig == SkipOriginal::True ) {
    //Always provides the original rng first, skip over it as requested (so it
    //won't be provided to anyone else):
    m_impl->produceUnlocked();
  }

  //Immediately set aside (as promised in class description) the first RNGStream
  //for the current thread. That way, registering a given RNG first and then
  //calling the global getRNG() function in the same thread, will give the
  //original stream back - which is probably what most people will expect in a
  //single threaded application:
  produceForCurrentThread();
}

NC::RNGProducer::RNGProducer( RNGProducer&& ) noexcept = default;
NC::RNGProducer& NC::RNGProducer::operator=( RNGProducer&& ) noexcept = default;
NC::RNGProducer::~RNGProducer() = default;

NC::shared_obj<NC::RNGStream> NC::getRNG()
{
  return getDefaultRNGProducer()->produceForCurrentThread();
}

NC::shared_obj<NC::RNGStream> NC::getIndependentRNG()
{
  return getDefaultRNGProducer()->produce();
}

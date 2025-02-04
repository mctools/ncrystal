#ifndef NCrystal_RNG_hh
#define NCrystal_RNG_hh

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

#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  //RNG stream and producer classes (defined further down):
  class RNGStream;
  class RNGProducer;

  //Modify the RNG default streams used by NCrystal:
  NCRYSTAL_API void setDefaultRNG( shared_obj<RNGStream> );
  NCRYSTAL_API void setDefaultRNGFctForAllThreads( std::function<double()> );
  NCRYSTAL_API void clearDefaultRNG();

  //For some applications it might be desirable to access the RNG streams
  //through these global methods:
  NCRYSTAL_API shared_obj<RNGStream> getRNG();//returned stream is unique to current thread
  NCRYSTAL_API shared_obj<RNGStream> getIndependentRNG();//returned stream is entirely independent
  NCRYSTAL_API shared_obj<RNGProducer> getDefaultRNGProducer();

  //It is also possible to simply create a completely new instance of the
  //builtin RNG. WARNING: This is not usually what you should do, if you care
  //about unbiased random numbers (or multi-thread safety).
  NCRYSTAL_API shared_obj<RNGStream> createBuiltinRNG( uint64_t seed = 0 );
  NCRYSTAL_API shared_obj<RNGStream> createBuiltinRNG( const RNGStreamState& state );

  //Check whether a given RNG state is from the builtin RNG:
  NCRYSTAL_API bool stateIsFromBuiltinRNG(const RNGStreamState&);

  // The RNGStream and RNGProducer classes are most likely to be of interest
  // to people registering their own RNG. For simply consuming random numbers,
  // refer to the base RNG interface in NCDefs.hh.

  class NCRYSTAL_API RNGStream : public RNG {
  public:

    ////////////////////////////////////////////////////////////////////////////
    // Class which extends the basic random number provision features of the  //
    // RNG class with additional interfaces for controlling the state etc. of //
    // a RNG stream. In general, actual RNG implementations which are to be   //
    // registered and used with NCrystal should inherit from RNGStream.       //
    ////////////////////////////////////////////////////////////////////////////

    //Some RNG's supports getting/setting state. For convenience this is is
    //exposed as printable strings only containing alphanumeric
    //characters. Classes should reimplement all 4 of the following functions
    //if any are reimplemented:
    bool supportsStateManipulation() const noexcept;
    RNGStreamState getState() const;
    void setState( const RNGStreamState& );
    shared_obj<RNGStream> cloneWithNewState( const RNGStreamState& ) const;

    //Some RNG's support "jumping" to a new independent sequence of random
    //numbers. This is ideally implemented by clever algorithms which can
    //cheaply bring the RNG to the same state that would otherwise have been
    //reached by a VERY high number of calls to generate() (e.g. 2e64). This
    //will be used to create on-demand multi-thread safe streams.
    virtual bool isJumpCapable() const { return false; }
    virtual shared_obj<RNGStream> createJumped() const;

    //Override and return true if the same RNG stream should be used even if
    //objects are cloned for different threads (this supports wrapping of
    //global RNG sources which already takes care of concurrency by other
    //means - e.g. in Geant4).
    virtual bool useInAllThreads() const { return false; }

    //Default implementations (can be reimplemented if more efficient methods
    //are available from a given generator):
    bool coinflip() override;
    uint64_t generate64RndmBits() override;
    uint32_t generate32RndmBits() override;

  protected:
    //If supporting state manipulation, please implement the following three
    //functions (type UID should be non-zero number unique to the RNG type and
    //version):
    virtual uint32_t stateTypeUID() const noexcept;
    virtual void actualSetState( std::vector<uint8_t>&& );
    virtual std::vector<uint8_t> actualGetState() const;
    virtual shared_obj<RNGStream> actualCloneWithNewState( std::vector<uint8_t>&& ) const;

  public:
    //Utilities for implementing actualSetState/actualGetState:
    template <class TInteger>
      static void appendToStateVector( std::vector<uint8_t>&, TInteger );
    template <class TInteger>
      static TInteger popFromStateVector( std::vector<uint8_t>& );
  };

  class NCRYSTAL_API RNGProducer final : private MoveOnly {
  public:

    ///////////////////////////////////////////////////////////////////////
    // Class which, based on the capabilities of a given RNGStream, can  //
    // produce new independent streams.                                  //
    //                                                                   //
    // It always first sets aside the initial rng stream for the current //
    // thread, so it can be "produced" by subsequently calling           //
    // produceForCurrentThread().                                        //
    ///////////////////////////////////////////////////////////////////////

    enum class SkipOriginal { True, False };
    explicit RNGProducer( shared_obj<RNGStream> rng,
                          SkipOriginal = SkipOriginal::False );

    //Produce new independent stream:
    shared_obj<RNGStream> produce();

    //Produce independent stream by index. Calling later with the same index
    //will return the same stream:
    shared_obj<RNGStream> produceByIdx( RNGStreamIndex );

    //Produce independent stream for current thread. All calls within a given
    //thread will return the same stream, but will give independent streams
    //for different threads:
    shared_obj<RNGStream> produceForCurrentThread();

    RNGProducer( RNGProducer&& ) noexcept;
    RNGProducer& operator=( RNGProducer&& ) noexcept;
    ~RNGProducer();

    //Reinit this producer:
    void reinit( shared_obj<RNGStream> rng, SkipOriginal so ) { *this = RNGProducer(std::move(rng),so); }

    //Create null producer, not capable of producing anything:
    explicit RNGProducer( no_init_t );

    //Shared instance of null producer:
    static shared_obj<RNGProducer> getNullProducer();

  private:
    struct Impl;
    Pimpl<Impl> m_impl;
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////


namespace NCRYSTAL_NAMESPACE {

  inline bool RNGStream::supportsStateManipulation() const noexcept
  {
    return stateTypeUID()!=0;
  }

  template <class TUInteger>
  inline void RNGStream::appendToStateVector( std::vector<uint8_t>& v, TUInteger val )
  {
    constexpr auto nbytes = sizeof(TUInteger);
    static_assert(nbytes>0,"");
    static_assert(std::is_unsigned<TUInteger>::value,"only works with unsigned integers");
    for ( decltype(sizeof(TUInteger)) i = 0; i < nbytes; ++i )
      v.push_back( static_cast<uint8_t>((val>>8*((nbytes-1)-i)) & 0xFF) );
  }

  template <class TUInteger>
  inline TUInteger RNGStream::popFromStateVector( std::vector<uint8_t>& v )
  {
    constexpr auto nbytes = sizeof(TUInteger);
    static_assert(nbytes>0,"");
    static_assert(std::is_unsigned<TUInteger>::value,"only works with unsigned integers");
    TUInteger val{0};
    nc_assert_always( v.size() >= nbytes );
    unsigned offset = 8*(nbytes-1);
    for ( auto it = std::prev(v.end(),nbytes); it!=v.end(); ++it, offset-=8 )
      val |= ( TUInteger{*it} << offset );
    v.resize( v.size() - nbytes );
    return val;
  }

}

#endif

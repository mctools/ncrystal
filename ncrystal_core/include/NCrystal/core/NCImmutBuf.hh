#ifndef NCrystal_ImmutBuf_hh
#define NCrystal_ImmutBuf_hh

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

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  /////////////////////////////////////////////////////////////////////////////////
  //                                                                             //
  // Immutable buffer with custom alignment and SBO (small buffer                //
  // optimisation). If needed, remote storage is kept in a shared ptr, so        //
  // copying of buffer objects will always happen safely and malloc's are        //
  // reduced. Note that the data length is NOT stored on the object, so the code //
  // using ImmutableBuffer's are responsible for not using data beyond the       //
  // length. It is optionally possible to store a trivial meta data object (such //
  // as an integer) on the object, which is both for convenience and to ensure   //
  // any padding between the buffer and the metadata object are used to increase //
  // the size of the local buffer. Minimum local buffer size and required buffer //
  // alignment should be specified as template parameters.  Note that a single   //
  // char is used internally keep track of local vs. remote mode, so optimal     //
  // LOCALBUF_MINSIZE will value will typically be 1 less than a power of two.   //
  //                                                                             //
  // NB: On most platforms, only powers of two (1, 2, 4, 8, 16, ...) are valid   //
  // BUF_ALIGNMENT choices, but on some platforms e.g. long double apparently    //
  // has alignment of 10, so we don't static assert on it. It is always          //
  // recommended to set BUF_ALIGNMENT to the alignment of the data to be stored, //
  // extracted using alignof(..). UPDATE 2023: The underlying allignedAlloc code //
  // now static_asserts on powers of two, so the code will need to be updated if //
  // we encounter such platforms.                                                //
  //                                                                             //
  /////////////////////////////////////////////////////////////////////////////////

  template< std::size_t LOCALBUF_MINSIZE = 15, std::size_t BUF_ALIGNMENT=1, class TMetaData = NullOptType >
  class ImmutableBuffer final {
    class RemoteBuf;
    using TRemotePtr = std::shared_ptr<const RemoteBuf>;
  public:

    using size_type = std::size_t;

    //Use this constructor if TMetaData is NOT provided:
    ImmutableBuffer( const char * src, size_type len );

    //Use this constructor if TMetaData IS provided:
    ImmutableBuffer( const char * src, size_type len, TMetaData md );

    //Use these constructors if it is needed to construct unset objects ( this
    //is normally only sensible if the object is soon after filled by move or
    //copy assignment). Any MetaData object will be default constructed.:
    ImmutableBuffer( NullOptType );

    //Access data:
    const char * data() const ncnoexceptndebug;
    const TMetaData& metaData() const ncnoexceptndebug;

    //When the calling code is certain that the data fits in the local buffer
    //(i.e. len provided to constructor was <= buffer_local_size ), the
    //following method is slightly faster than data(). In debug builds, wrong
    //usage will trigger a runtime exception:
    const char * dataAssertLocal() const ncnoexceptndebug;

    //Destructor (does no work unless in remote mode):
    ~ImmutableBuffer();

    //Both copy and move semantics are supported. When the buffer is local, this
    //is usually implemented via a std::memcpy, otherwise it is more like a
    //copy/move of a std::shared_ptr.
    ImmutableBuffer( const ImmutableBuffer& ) noexcept;
    ImmutableBuffer& operator=( const ImmutableBuffer& ) noexcept;
    ImmutableBuffer( ImmutableBuffer&& ) noexcept;
    ImmutableBuffer& operator=( ImmutableBuffer&& ) noexcept;

    using metadata_type = TMetaData;
    static constexpr bool has_metadata = ! std::is_same<TMetaData,NullOptType>::value;
    static constexpr size_type buffer_alignment = BUF_ALIGNMENT;
    static constexpr size_type object_alignment = ncconstexpr_lcm<size_type>( BUF_ALIGNMENT,
                                                                              alignof(TRemotePtr),
                                                                              ( has_metadata ? alignof(TMetaData) : 1 ) );
    static constexpr size_type metadata_size = ( has_metadata
                                                 ? ncconstexpr_roundToNextMultipleOf<size_type>( sizeof(TMetaData), alignof(TMetaData) )
                                                 : 0 );
    static constexpr size_type unused_trailing_bytes = ( metadata_size > sizeof(TMetaData)
                                                         ? metadata_size - sizeof(TMetaData)
                                                         : 0 );
    static constexpr size_type object_size = ncconstexpr_roundToNextMultipleOf<size_type>( 1 + metadata_size +
                                                                                           ncconstexpr_max<size_type>( sizeof(TRemotePtr),
                                                                                                                       LOCALBUF_MINSIZE ),
                                                                                           object_alignment );
    static constexpr size_type buffer_local_size = object_size - metadata_size - 1;
    static constexpr size_type buffer_local_requested_size = LOCALBUF_MINSIZE;

  private:
    static constexpr size_type _detail_pos_md = object_size - metadata_size;
    static constexpr size_type _detail_pos_mode = buffer_local_size;

    //Implemenation notes: In remote mode, the local storage buffer is used to
    //keep the TRemotePtr. Thus, local storage has size and alignment suitable
    //for such a pointer, as well as for the user data in local mode. Finally,
    //the char just after the local buffer is used to keep the status flag
    //indicating remote/local/unset mode, and the object size will always be as
    //least as big as its alignment (to avoid padding the objects are kept in an
    //array).

    //All data in a single char array to avoid unexpected padding:
    alignas(object_alignment) char m_data[object_size];

    void initBuffer( const char *, size_type );
    constexpr bool isRemote() const noexcept;
    constexpr bool isEmpty() const noexcept;
    const TRemotePtr& remotePtr() const noexcept;
    TRemotePtr& remotePtr() noexcept;
    TMetaData& mdRef() noexcept;
    void unsetIfRemote() noexcept;
    const TMetaData& detail_metaData_nocheck() const noexcept;
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  template<std::size_t LOCALBUF_MINSIZE, std::size_t BUF_ALIGNMENT, class TMetaData>
  inline ImmutableBuffer<LOCALBUF_MINSIZE,BUF_ALIGNMENT,TMetaData>::~ImmutableBuffer()
  {
    static_assert( LOCALBUF_MINSIZE >= 1, "" );
    static_assert( LOCALBUF_MINSIZE <= 1024*1024, "" );
    static_assert( BUF_ALIGNMENT >= 1, "" );
    static_assert( BUF_ALIGNMENT <= 1024*1024, "" );
    static_assert( std::is_trivial<TMetaData>::value, "");
    static_assert( std::is_trivially_destructible<TMetaData>::value, "");
    static_assert( std::is_trivially_copyable<TMetaData>::value, "");
    static_assert( std::is_nothrow_destructible<TMetaData>::value, "");
    static_assert( std::is_nothrow_default_constructible<TMetaData>::value, "");
    static_assert( std::is_nothrow_copy_constructible<TMetaData>::value, "");
    static_assert( std::is_nothrow_copy_assignable<TMetaData>::value, "");
    static_assert( std::is_nothrow_move_constructible<TMetaData>::value, "");
    static_assert( std::is_nothrow_move_assignable<TMetaData>::value, "");
    static_assert( _detail_pos_md > 1, "" );
    static_assert( object_size == buffer_local_size + 1 + metadata_size, "" );
    static_assert( LOCALBUF_MINSIZE <= buffer_local_size, "" );
    static_assert(sizeof(char)==sizeof(uint8_t), "");

    //Release shared reference if remote mode, otherwise do nothing:
    unsetIfRemote();
  }

  template<std::size_t LBMS, std::size_t BA, class TMetaData>
  inline ImmutableBuffer<LBMS,BA,TMetaData>::ImmutableBuffer( NullOptType )
  {
    m_data[_detail_pos_mode] = 0;
    if ( has_metadata )
      *reinterpret_cast<TMetaData*>(&m_data[_detail_pos_md]) = TMetaData();
  }

  template<std::size_t LBMS, std::size_t BA, class MD>
  inline ImmutableBuffer<LBMS,BA,MD>::ImmutableBuffer( const char * src, std::size_t len )
  {
    static_assert( !has_metadata, "MetaData enabled and must be provided in constructor." );
    initBuffer(src,len);
  }

  template<std::size_t LOCALBUF_MINSIZE, std::size_t BUF_ALIGNMENT, class TMetaData>
  class ImmutableBuffer<LOCALBUF_MINSIZE,BUF_ALIGNMENT,TMetaData>::RemoteBuf : private NoCopyMove {
    char * m_data;
  public:
    RemoteBuf(const char * src, size_type len ) : m_data((char*)AlignedAlloc::alignedAllocFixedAlign<BUF_ALIGNMENT>(len)) { std::memcpy(m_data,src,len); }
    ~RemoteBuf() { AlignedAlloc::freeAlignedAllocFixedAlign<BUF_ALIGNMENT>(static_cast<void*>(m_data)); }
    const char * data() const noexcept { return m_data; }
  };

  template<std::size_t LBMS, std::size_t BA, class MD>
  inline void ImmutableBuffer<LBMS,BA,MD>::unsetIfRemote() noexcept {
    if (isRemote()) {
      TRemotePtr* rp = &remotePtr();
      m_data[_detail_pos_mode] = 0;
      rp->~TRemotePtr();
    }
  }

  template<std::size_t LBMS, std::size_t BA, class TMetaData>
  inline ImmutableBuffer<LBMS,BA,TMetaData>::ImmutableBuffer( const char * src, std::size_t len, TMetaData md )
  {
    initBuffer(src,len);
    *reinterpret_cast<TMetaData*>(&m_data[_detail_pos_md]) = std::move(md);
  }

  template<std::size_t LBMS, std::size_t BA, class MD>
  inline void ImmutableBuffer<LBMS,BA,MD>::initBuffer( const char * src, std::size_t len )
  {
    nc_assert( src != nullptr );//sanity
    nc_assert( len >= 1 );//sanity
    nc_assert( len < 1000000000 );//we take >1GB is sign of invalid len parameter (e.g. result of integer overflow).
    if ( len > buffer_local_size ) {
      static_assert(sizeof(TRemotePtr)+1<=sizeof(m_data), "");
      new(&m_data[0]) TRemotePtr(std::move(std::make_shared<const RemoteBuf>(src,len)));
      m_data[_detail_pos_mode] = (char)1;//remote mode
    } else {
      std::memcpy(&m_data[0],src,len);
      m_data[_detail_pos_mode] = (char)2;//local mode
    }
  }

  template<std::size_t LBMS, std::size_t BA, class MD>
  inline const char * ImmutableBuffer<LBMS,BA,MD>::data() const ncnoexceptndebug
  {
#ifndef NDEBUG
    if ( isEmpty() )
      NCRYSTAL_THROW(LogicError,"Trying to access data() on moved-from ImmutableBuffer object!");
#endif
    return ( isRemote() ? remotePtr()->data() : &m_data[0] );
  }

  template<std::size_t LBMS, std::size_t BA, class MD>
  inline const char * ImmutableBuffer<LBMS,BA,MD>::dataAssertLocal() const ncnoexceptndebug
  {
#ifndef NDEBUG
    if ( isEmpty() )
      NCRYSTAL_THROW(LogicError,"Trying to access dataAssertLocal() on moved-from ImmutableBuffer object!");
#endif
#ifndef NDEBUG
    if ( isRemote() )
      NCRYSTAL_THROW(LogicError,"Trying to access dataAssertLocal() on ImmutableBuffer object which is using remote buffer!");
#endif
    return &m_data[0];
  }

  template<std::size_t LBMS, std::size_t BA, class TMetaData>
  inline const TMetaData& ImmutableBuffer<LBMS,BA,TMetaData>::metaData() const ncnoexceptndebug
  {
#ifndef NDEBUG
    if ( isEmpty() )
      NCRYSTAL_THROW(LogicError,"Trying to access metaData() on moved-from ImmutableBuffer object!");
#endif
    return *reinterpret_cast<const TMetaData*>(&m_data[_detail_pos_md]);
  }

  template<std::size_t LBMS, std::size_t BA, class TMetaData>
  inline const TMetaData& ImmutableBuffer<LBMS,BA,TMetaData>::detail_metaData_nocheck() const noexcept
  {
    return *reinterpret_cast<const TMetaData*>(&m_data[_detail_pos_md]);
  }

  template<std::size_t LBMS, std::size_t BA, class MD>
  inline ImmutableBuffer<LBMS,BA,MD>::ImmutableBuffer( const ImmutableBuffer& o ) noexcept
  {
    m_data[_detail_pos_mode] = (char)0;//unset
    *this = o;//copy via assignment
  }

  template<std::size_t LBMS, std::size_t BA, class MD>
  inline constexpr bool ImmutableBuffer<LBMS,BA,MD>::isRemote() const noexcept { return m_data[_detail_pos_mode] == 1; }
  template<std::size_t LBMS, std::size_t BA, class MD>
  inline constexpr bool ImmutableBuffer<LBMS,BA,MD>::isEmpty() const noexcept { return m_data[_detail_pos_mode] == 0; }
  template<std::size_t LBMS, std::size_t BA, class MD>
  inline const typename ImmutableBuffer<LBMS,BA,MD>::TRemotePtr& ImmutableBuffer<LBMS,BA,MD>::remotePtr() const noexcept { /*?nc_assert(isRemote());*/ return *reinterpret_cast<const TRemotePtr*>(&m_data[0]); }
  template<std::size_t LBMS, std::size_t BA, class MD>
  inline typename ImmutableBuffer<LBMS,BA,MD>::TRemotePtr& ImmutableBuffer<LBMS,BA,MD>::remotePtr() noexcept { /*nc_assert(isRemote());*/ return *reinterpret_cast<TRemotePtr*>(&m_data[0]); }
  template<std::size_t LBMS, std::size_t BA, class TMetaData>
  inline TMetaData& ImmutableBuffer<LBMS,BA,TMetaData>::mdRef() noexcept { return *reinterpret_cast<TMetaData*>(&m_data[_detail_pos_md]); }

  template<std::size_t LBMS, std::size_t BA, class TMetaData>
  inline auto ImmutableBuffer<LBMS,BA,TMetaData>::operator=( const ImmutableBuffer& o ) noexcept -> ImmutableBuffer&
  {
    if ( !o.isRemote() ) {
      //optimise for local case. Entire object (including mode char and metadata
      //object) can then be trivially copied.
      unsetIfRemote();//(but first make ourselves non-remote if needed)
      std::memcpy(&m_data[0],&o.m_data[0], object_size );
    } else {
      //o is remote, we must copy its remote ptr..
      if (!isRemote()) {
        new(&m_data[0]) TRemotePtr();
        m_data[_detail_pos_mode] = (char)1;//remote mode
      }
      remotePtr() = o.remotePtr();
      if ( has_metadata )
        mdRef() = o.detail_metaData_nocheck();
    }
    return* this;
  }

  template<std::size_t LBMS, std::size_t BA, class MD>
  inline ImmutableBuffer<LBMS,BA,MD>::ImmutableBuffer( ImmutableBuffer&& o ) noexcept
  {
    m_data[_detail_pos_mode] = (char)0;//unset
    *this = std::move(o);//move via assignment
  }

  template<std::size_t LBMS, std::size_t BA, class TMetaData>
  inline auto ImmutableBuffer<LBMS,BA,TMetaData>::operator=( ImmutableBuffer&& o ) noexcept -> ImmutableBuffer&
  {
    //NB: always leave moved-from object unset for consistency
    if ( !o.isRemote() ) {
      //Easy, just copy over entire object byte by byte
      unsetIfRemote();//(but unset our own remote ptr first if set)
      std::memcpy(&m_data[0],&o.m_data[0],object_size);
      //always leave o unset for consistency (o was not in remote mode so this
      //is easy):
      o.m_data[_detail_pos_mode] = (char)0;
    } else {
      //o is remote. We move the shared ptr in o into our state (which might
      //already contain a remoteptr).
      if ( !isRemote() ) {
        new(&m_data[0]) TRemotePtr();
        m_data[_detail_pos_mode] = (char)1;//remote mode
      }
      remotePtr() = std::move(o.remotePtr());
      if ( has_metadata )
        mdRef() = std::move(o.mdRef());
      o.unsetIfRemote();
    }
    return *this;
  }

}

#endif

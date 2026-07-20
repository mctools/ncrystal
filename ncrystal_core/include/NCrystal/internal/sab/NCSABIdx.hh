#ifndef NCrystal_SABIdx_hh
#define NCrystal_SABIdx_hh

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

  //Fixme: major cleanup and inlining needed.

  namespace SABIdx {

    /////////////////////////////////////
    // Utilities for SAB Cell indexing //
    /////////////////////////////////////

    // A tabulated S(alpha,beta) table like the one in a SABData object requires
    // several kinds of indexing. A few complications and considerations:
    //
    //  * We have nalpha*nbeta grid points, but only (nalpha-1)*(nbeta-1) cells.
    //  * For some usage, the indices along one axis, ialpha and ibeta, are
    //    needed, and for other usage the indices into the arrays of S values
    //    ("sab" or "S", or even "logS") is needed.
    //  * For some usage, the most efficient choice would be to keep indices in
    //    std::size_t (for e.g. usage with the std::vector interface), while
    //    uint32_t has a more tight memory layout (which of course impose range
    //    limitations).
    //  * We recall that the SABData convention is (nalpha=alphagrid.size()):
    //    S(alphaGrid()[ia],betaGrid()[ib])=sab()[ib*nalpha+ia]
    //
    // For packing (ia,ib) into a single uint, we could either use
    // i=ib*nalpha+ia, which is exactly the index needed to access the S value
    // in the sab() array. Or we could simply pack ia and ib into their own half
    // of the bits (e.g. "ia&(ib<<16)" ). This imposes a more restrictive
    // maximum supported grid size, but makes it much more efficient to extract
    // ia and ib, which would otherwise need an expensive integer division by
    // nalpha. Given that final usage for e.g. looking up all data of a given
    // cell is likely to need ia and ib (to get cell alpha and beta values), we
    // will prefer the latter. We will also prioritise memory efficiency, so we
    // will require nalpha and nbeta to fit in std::uint16_t (~65k), and store
    // both into an uint32_t via a simply bitshift.


    //Basic integer types, not strongly typed for purpose:
    using raw_idx_t = std::uint32_t;
    using raw_size_t = std::size_t;
    using raw_sab_idx_t = raw_size_t;//fixme: this presumes that we do not
                                     //actually store it, it should actually be
                                     //safe to store in uint32_t

    //Fwd declarations:
    class PackedIndex;//packing (ialpha,ibeta) into a single integer
    class NAlpha;//strongly typed "nalpha"
    class NAlphaCells;//strongly typed "nalpha-1"
    class SABIndex;//Strongly typed SAB index ( idx = ib*nalpha+ia )
    class SABCellIndex;//Strongly typed SAB index ( idx = ib*(nalpha-1)+ia )

    class NAlpha final {
    public:
      //Strongly typed NAlpha ("alphaGrid.size()")
      using size_type = raw_size_t;
      NAlpha( const VectD& alphaGrid ) ncnoexceptndebug;
      explicit constexpr NAlpha( size_type ) ncnoexceptndebug;
      explicit constexpr NAlpha( NAlphaCells ) ncnoexceptndebug;
      constexpr size_type value() const noexcept { return m_value; }
    private:
      size_type m_value;
    };

    class NAlphaCells final {
    public:
      //Strongly typed NAlphaCells ("alphaGrid.size()-1")
      using size_type = raw_size_t;
      NAlphaCells( const VectD& alphaGrid ) ncnoexceptndebug;
      explicit constexpr NAlphaCells( size_type ) ncnoexceptndebug;
      explicit constexpr NAlphaCells( NAlpha ) ncnoexceptndebug;
      constexpr size_type value() const noexcept { return m_value; }
    private:
      size_type m_value;
    };

    namespace detail {

      template<class TTwo, class TOne>
      struct UIntPacker final {
        //pack two uint16 into one uint32, or two uint32 into one uint64.
        static constexpr int shift = ( std::is_same<TOne,std::uint16_t>::value
                                       ? 16 : 32 );
        static constexpr TTwo mask = ( std::is_same<TOne,std::uint16_t>::value
                                       ? 0xFFFFu : 0xFFFFFFFFu );
        static constexpr std::size_t max();
        static TTwo pack( TOne val1, TOne val2 ) noexcept
        {
          return ( static_cast<TTwo>(val1) << shift ) | static_cast<TTwo>(val2);
        }
        static constexpr TTwo unpack1( TTwo p ) noexcept { return p >> shift; }
        static constexpr TTwo unpack2( TTwo p ) noexcept { return p & mask; }

        //For packing to/from std::size_t:
        static TTwo pack_sizet( std::size_t, std::size_t ) ncnoexceptndebug;
        static constexpr std::size_t unpack1_sizet( TTwo ) noexcept;
        static constexpr std::size_t unpack2_sizet( TTwo ) noexcept;
        static constexpr std::size_t max_sizet() noexcept;
      };
    }

    class PackedIndex {
      using packer_t = detail::UIntPacker<raw_idx_t,std::uint16_t>;
    public:
      //Cell indices are in unpacked form (ialpha,ibeta) where the indices are
      //the ones of the lower grid point. Thus, if the grid has nalpha*nbeta
      //grid points, there are only (nalpha-1)*(nbeta-1) cells. We thus define
      //nalphacells=nalpha-1 and ncellbeta=nbeta-1, and let ialpha run from
      //0..(nalphacells-1), and similarly for ibeta.

      using index_t = raw_idx_t;

      static constexpr index_t invalid = std::numeric_limits<index_t>::max();
      constexpr PackedIndex() noexcept : m_idx(invalid) {}
      constexpr PackedIndex( NullOptType ) noexcept : m_idx(invalid) {}
      constexpr PackedIndex( raw_size_t ialpha,
                             raw_size_t ibeta ) ncnoexceptndebug;

      std::size_t idxAlpha() const
      {
        return packer_t::unpack1_sizet( m_idx );
      }

      std::size_t idxBeta() const
      {
        return packer_t::unpack2_sizet( m_idx );
      }

      //      CellIndex unpack( NAlphaCells ) const ncnoexceptndebug;
      constexpr bool isValid() const noexcept { return m_idx != invalid; }
      constexpr index_t value() const noexcept { return m_idx; }
      //static index_t nRawMax( NAlphaCells, std::size_t nbeta );
      struct from_raw_t {};
      constexpr PackedIndex( from_raw_t, index_t raw ) noexcept : m_idx(raw) { }
      //constexpr index_t rawValue() const noexcept { return m_idx; }
    private:
      index_t m_idx;
    };

    template<class NAlphaT>
    class SABIdxT {
    public:
      using size_type = raw_sab_idx_t;
      SABIdxT( NAlpha, size_type idxAlpha, size_type idxBeta ) ncnoexceptndebug;
      SABIdxT( NAlpha na, PackedIndex p ) ncnoexceptndebug
        : SABIdxT( na, p.idxAlpha(), p.idxBeta() )
      {
        nc_assert( p.isValid() );
      }
      explicit constexpr SABIdxT( size_type raw ) noexcept : m_value(raw)
      {
        nc_assert( m_value < static_cast<std::size_t>
                   (std::numeric_limits<std::uint32_t>::max()) );
      }
      constexpr size_type value() const noexcept { return m_value; }
    private:
      size_type m_value;
    };

    using SABIdx = SABIdxT<NAlpha>;
    using SABCellIdx = SABIdxT<NAlphaCells>;

  }
}

////////////////////////////
// Inline implementations //
////////////////////////////


namespace NCRYSTAL_NAMESPACE {
  namespace SABIdx {
    inline NAlphaCells::NAlphaCells( const VectD& alphaGrid ) ncnoexceptndebug
      : NAlphaCells( static_cast<size_type>(alphaGrid.size()-1) )
    {
      nc_assert( alphaGrid.size()>1);
    }
    inline constexpr NAlphaCells::NAlphaCells( size_type value ) ncnoexceptndebug
      : m_value( value )
    {
      nc_assert( (m_value+1)<std::numeric_limits<size_type>::max() );
      nc_assert( value>=1 );
    }
    inline constexpr NAlphaCells::NAlphaCells( NAlpha na ) ncnoexceptndebug
      : NAlphaCells( na.value()-1 )
    {
      nc_assert( na.value()>1);
    }
    inline NAlpha::NAlpha( const VectD& alphaGrid ) ncnoexceptndebug
      : NAlpha( static_cast<size_type>(alphaGrid.size()) )
    {
      nc_assert( alphaGrid.size()>1);
    }
    inline constexpr NAlpha::NAlpha( size_type value ) ncnoexceptndebug
      : m_value( value )
    {
      nc_assert( m_value<std::numeric_limits<size_type>::max() );
      nc_assert( value > 1 );
    }
    inline constexpr NAlpha::NAlpha( NAlphaCells nac ) ncnoexceptndebug
      : NAlpha( nac.value()+1 )
    {
      nc_assert( nac.value()>=1);
    }

    template<class TTwo, class TOne>
    inline constexpr std::size_t detail::UIntPacker<TTwo,TOne>::max()
    {
      static_assert( sizeof(TTwo) == 2*sizeof(TOne), "" );
      static_assert( std::is_same<TTwo,std::uint32_t>::value
                     || std::is_same<TTwo,std::uint64_t>::value, "" );
      static_assert( std::is_same<TOne,std::uint16_t>::value
                     || std::is_same<TOne,std::uint32_t>::value, "" );
      static_assert( sizeof(std::size_t) >= sizeof( TTwo ) );
      return static_cast<std::size_t>(std::numeric_limits<TOne>::max());
    }

    template<class TTwo, class TOne>
    inline TTwo detail::UIntPacker<TTwo,TOne>
    ::pack_sizet( std::size_t val1, std::size_t val2 ) ncnoexceptndebug
    {
      nc_assert( val1 <= max() );
      nc_assert( val2 <= max() );
      return pack( static_cast<TOne>(val1), static_cast<TOne>(val2) );
    }

    template<class TTwo, class TOne>
    inline constexpr std::size_t
    detail::UIntPacker<TTwo,TOne>::unpack1_sizet( TTwo p ) noexcept
    {
      return static_cast<std::size_t>(unpack1(p));
    }

    template<class TTwo, class TOne>
    inline constexpr std::size_t
    detail::UIntPacker<TTwo,TOne>::unpack2_sizet( TTwo p ) noexcept
    {
      return static_cast<std::size_t>(unpack2(p));
    }

    template<class TTwo, class TOne>
    inline constexpr std::size_t
    detail::UIntPacker<TTwo,TOne>::max_sizet() noexcept
    {
      return static_cast<std::size_t>(max());
    }

    inline constexpr PackedIndex
    ::PackedIndex( raw_size_t ialpha, raw_size_t ibeta ) ncnoexceptndebug
      : m_idx( packer_t::pack( ialpha, ibeta ) )
    {
      //disallow ialpha=ibeta=uint16_t::max, reserving it for "invalid":
      static_assert( packer_t::unpack1(invalid) == packer_t::max(), "" );
      static_assert( packer_t::unpack2(invalid) == packer_t::max(), "" );
      nc_assert( ialpha < packer_t::max_sizet() );
      nc_assert( ibeta < packer_t::max_sizet() );
      nc_assert( isValid() );
    }

    template<class NAlphaT>
    inline SABIdxT<NAlphaT>::SABIdxT( NAlpha na,
                                      size_type idxAlpha,
                                      size_type idxBeta ) ncnoexceptndebug
      : m_value( idxBeta * na.value() + idxAlpha )
    {
      nc_assert( idxAlpha < static_cast<std::size_t>
                 (std::numeric_limits<std::uint16_t>::max()) );
      nc_assert( idxBeta < static_cast<std::size_t>
                 (std::numeric_limits<std::uint16_t>::max()) );
      nc_assert( m_value < static_cast<std::size_t>
                 (std::numeric_limits<std::uint32_t>::max()) );
    }
  }
}
#endif

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
    // i=ib*nalpha+ia, which is also exactly the index needed to access the S
    // value in the sab() array. Or we could simply pack ia and ib into their
    // own half of the bits (e.g. "ia&(ib<<16)" ). This imposes a more
    // restrictive maximum supported grid size, but makes it much more efficient
    // to extract ia and ib, which would otherwise need an expensive integer
    // division by nalpha. Given that final usage for e.g. looking up all data
    // of a given cell is likely to need ia and ib (to get cell alpha and beta
    // values), we will prefer the latter. We will also prioritise memory
    // efficiency, so we will require nalpha and nbeta to fit in std::uint16_t
    // (~65k), and store both into an uint32_t via a simple bitshift.

    ////////////////////////////////////////////////////////////////
    //Call to trigger exception if out of range:
    void verifyGridSizes( std::size_t nalpha, std::size_t nbeta );

    ////////////////////////////////////////////////////////////////
    // Strongly typed SAB cell or grid point index. Deliberately
    // kept simple to be trivially constructible. This is
    // essentially just an uint32_t with utilities for very fast
    // packing/unpacking of (ialpha,ibeta), using just bit shifts
    // and masks:
    struct PackedIndex final {
      using index_t = std::uint32_t;
      index_t val;
      ncconstexprndebug std::size_t unpackAlphaIdx() const ncnoexceptndebug;
      ncconstexprndebug std::size_t unpackBetaIdx() const ncnoexceptndebug;

      static ncconstexprndebug
      PackedIndex fromAlphaIdxBetaIdx( std::size_t ialpha,
                                       std::size_t ibeta ) ncnoexceptndebug;
      static constexpr std::size_t maxAlphaIdx() noexcept;
      static constexpr std::size_t maxBetaIdx() noexcept;
    };

    ////////////////////////////////////////////////////////////////
    // Strongly typed number of grid points along alpha, which is
    // essentially the same as "alphaGrid.size()":
    class NAlphaCells;
    class NAlpha final {
    public:
      using size_type = std::size_t;
      NAlpha( const VectD& alphaGrid ) ncnoexceptndebug;
      explicit ncconstexprndebug NAlpha( size_type ) ncnoexceptndebug;
      explicit ncconstexprndebug NAlpha( NAlphaCells ) ncnoexceptndebug;
      ncconstexprndebug size_type value() const noexcept { return m_value; }
    private:
      size_type m_value;
    };

    ////////////////////////////////////////////////////////////////
    // Strongly typed number of cells along alpha, which is
    // essentially the same as "alphaGrid.size()-1":
    class NAlphaCells final {
    public:
      //Strongly typed NAlphaCells ("alphaGrid.size()-1")
      using size_type = std::size_t;
      NAlphaCells( const VectD& alphaGrid ) ncnoexceptndebug;
      explicit ncconstexprndebug NAlphaCells( size_type ) ncnoexceptndebug;
      explicit ncconstexprndebug NAlphaCells( NAlpha ) ncnoexceptndebug;
      ncconstexprndebug size_type value() const noexcept { return m_value; }
    private:
      size_type m_value;
    };

    ////////////////////////////////////////////////////////////////
    // Strongly typed index into a 1D vector of information about grid
    // points or cells. Examples of such vectors would be the S values
    // themselves at each grid point, or the integral of S over each
    // cell. This is essentially calculated as betaIdx * N + alphaIdx,
    // where N is either NAlpha or NAlphaCells:
    template<class NAlphaT>
    class SABIdxT {
    public:
      using size_type = std::size_t;
      ncconstexprndebug SABIdxT( NAlphaT,
                                 size_type idxAlpha,
                                 size_type idxBeta ) ncnoexceptndebug;
      ncconstexprndebug SABIdxT( NAlphaT na,
                                 PackedIndex p ) ncnoexceptndebug;
      ncconstexprndebug size_type value() const noexcept { return m_value; }
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
    inline ncconstexprndebug NAlphaCells::NAlphaCells( size_type value ) ncnoexceptndebug
      : m_value( value )
    {
#ifndef NDEBUG
      nc_assert( (m_value+1)<std::numeric_limits<size_type>::max() );
      nc_assert( value>=1 );
#endif
    }
    inline ncconstexprndebug NAlphaCells::NAlphaCells( NAlpha na ) ncnoexceptndebug
      : NAlphaCells( na.value()-1 )
    {
#ifndef NDEBUG
      nc_assert( na.value()>1);
#endif
    }
    inline NAlpha::NAlpha( const VectD& alphaGrid ) ncnoexceptndebug
      : NAlpha( static_cast<size_type>(alphaGrid.size()) )
    {
#ifndef NDEBUG
      nc_assert( alphaGrid.size()>1);
#endif
    }
    inline ncconstexprndebug NAlpha::NAlpha( size_type value ) ncnoexceptndebug
      : m_value( value )
    {
#ifndef NDEBUG
      nc_assert( m_value<std::numeric_limits<size_type>::max() );
      nc_assert( value > 1 );
#endif
    }
    inline ncconstexprndebug NAlpha::NAlpha( NAlphaCells nac ) ncnoexceptndebug
      : NAlpha( nac.value()+1 )
    {
#ifndef NDEBUG
      nc_assert( nac.value()>=1);
#endif
    }

    namespace detail {

      template<class TTwo, class TOne>
      struct UIntPacker final {
        //pack two uint16 into one uint32, or two uint32 into one uint64.
        static constexpr int shift = ( std::is_same<TOne,std::uint16_t>::value
                                       ? 16 : 32 );
        static constexpr TTwo mask = ( std::is_same<TOne,std::uint16_t>::value
                                       ? 0xFFFFu : 0xFFFFFFFFu );
        static constexpr TOne component_max();
        using packed_type = TTwo;
        using component_type = TOne;
        static constexpr TTwo pack( TOne val1, TOne val2 ) noexcept;
        static constexpr TTwo unpack1( TTwo p ) noexcept { return p >> shift; }
        static constexpr TTwo unpack2( TTwo p ) noexcept { return p & mask; }
      };
      using packer32_t = UIntPacker<std::uint32_t,std::uint16_t>;
    }

    template<class TTwo, class TOne>
    inline constexpr TTwo
    detail::UIntPacker<TTwo,TOne>::pack( TOne val1, TOne val2 ) noexcept
    {
      static_assert( sizeof(TTwo) == 2*sizeof(TOne), "" );
      static_assert( std::is_same<TTwo,std::uint32_t>::value
                     || std::is_same<TTwo,std::uint64_t>::value, "" );
      static_assert( std::is_same<TOne,std::uint16_t>::value
                     || std::is_same<TOne,std::uint32_t>::value, "" );
      static_assert( sizeof(std::size_t) >= sizeof( TTwo ), "" );
      return ( static_cast<TTwo>(val1) << shift ) | static_cast<TTwo>(val2);
    }

    template<class TTwo, class TOne>
    inline constexpr TOne detail::UIntPacker<TTwo,TOne>::component_max()
    {
      return std::numeric_limits<TOne>::max();
    }

    inline constexpr std::size_t PackedIndex::maxAlphaIdx() noexcept
    {
      //-1 to potentially leave room for an "invalid" state:
      return static_cast<std::size_t>(detail::packer32_t::component_max()-1);
    }

    inline constexpr std::size_t PackedIndex::maxBetaIdx() noexcept
    {
      return maxAlphaIdx();
    }

    inline ncconstexprndebug std::size_t
    PackedIndex::unpackAlphaIdx() const ncnoexceptndebug
    {
#ifndef NDEBUG
      nc_assert( static_cast<std::size_t>(detail::packer32_t::unpack1( val ))
                 <= maxAlphaIdx() );
#endif
      return static_cast<std::size_t>(detail::packer32_t::unpack1( val ));
    }

    inline ncconstexprndebug std::size_t
    PackedIndex::unpackBetaIdx() const ncnoexceptndebug
    {
#ifndef NDEBUG
      nc_assert( static_cast<std::size_t>(detail::packer32_t::unpack2( val ))
                 <= maxBetaIdx() );
#endif
      return static_cast<std::size_t>(detail::packer32_t::unpack2( val ));
    }

    inline ncconstexprndebug PackedIndex
    PackedIndex::fromAlphaIdxBetaIdx( std::size_t ialpha,
                                      std::size_t ibeta ) ncnoexceptndebug
    {
      static_assert( std::is_same<detail::packer32_t::packed_type,
                     PackedIndex::index_t>::value, "" );
#ifndef NDEBUG
      nc_assert( ialpha <= maxAlphaIdx() );
      nc_assert( ibeta <= maxBetaIdx() );
#endif
      return PackedIndex{ detail::packer32_t::pack( ialpha, ibeta ) };
    }

    inline void verifyGridSizes( std::size_t nalpha, std::size_t nbeta )
    {
      //Fixme: better location for this? In SABData class? Or
      //TransformKnlToStdFormat? Also, not inlined...
      static_assert(std::is_same<PackedIndex::index_t,std::uint32_t>::value,"");
      static_assert( std::numeric_limits<std::uint16_t>::max()==65535, "" );
      static_assert( PackedIndex::maxAlphaIdx()==65534, "" );
      static_assert( PackedIndex::maxBetaIdx()==65534, "" );
      if ( !( nalpha <= 65534 && nbeta <= 65534 ) )
        NCRYSTAL_THROW2(BadInput,"SAB grid too large"
                        " (max size is 65534x65534)");
      if ( !( nalpha >= 2 && nbeta >= 2 ) )
        NCRYSTAL_THROW2(BadInput,"SAB grid too small (min size is 2x2)");
    }

    template<class NAlphaT>
    inline ncconstexprndebug
    SABIdxT<NAlphaT>::SABIdxT( NAlphaT na,
                               size_type idxAlpha,
                               size_type idxBeta ) ncnoexceptndebug
      : m_value( idxBeta * na.value() + idxAlpha )
    {
#ifndef NDEBUG
      nc_assert( idxAlpha < static_cast<std::size_t>
                 (std::numeric_limits<std::uint16_t>::max()) );
      nc_assert( idxBeta < static_cast<std::size_t>
                 (std::numeric_limits<std::uint16_t>::max()) );
      nc_assert( m_value < static_cast<std::size_t>
                 (std::numeric_limits<std::uint32_t>::max()) );
#endif
    }

    template<class NAlphaT>
    inline ncconstexprndebug
    SABIdxT<NAlphaT>::SABIdxT( NAlphaT na,
                               PackedIndex p ) ncnoexceptndebug
      : SABIdxT( na, p.unpackAlphaIdx(), p.unpackBetaIdx() )
    {
    }
  }
}
#endif

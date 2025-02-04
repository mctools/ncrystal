#ifndef NCrystal_StrView_hh
#define NCrystal_StrView_hh

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

#include "NCrystal/core/NCSmallVector.hh"

namespace NCRYSTAL_NAMESPACE {

  class StrView;
  class StrView
  {
    // Simple string view class similar to std::string_view (which is
    // unavailable for C++11 code). Note however that this is an immutable
    // class, and provides a constant-view only. Also note that the StrView can
    // be without value, which is a different state than simply having a viewing
    // of an empty (size=0) string (this behavour is different than a regular
    // std::string_view).
    //
    // Note the usual caveats for a string view: The view does not participate
    // in life-time management of the pointed-too data, and calling any method
    // on a StrView whose pointed-too data is expired or modified, is undefined
    // behaviour (with the exception of StrView assignment operators or
    // destructor which does not actually access the data). Also be aware that
    // the data it points to is not necessarily null-terminated

  public:
    using size_type = std::size_t;
    static constexpr size_type npos = std::numeric_limits<size_type>::max();

    //Constructing from specific data range is simply setting the two internal
    //data members:
    constexpr StrView( const char * cdata, size_type len) noexcept;
    constexpr StrView( const StrView& ) noexcept;
    constexpr StrView( StrView&& ) noexcept;
    StrView( const std::string& ) noexcept;
    StrView( std::string&& ) = delete;//unsafe!

    //Constructing from a char-array does not need a call to strlen:
    template<size_type N>
    constexpr StrView( const char (&cdata)[N] ) noexcept
      : StrView( &cdata[0], ( cdata[N-1]=='\0' ? N-1 : N ) ) {}

    //Construction from a null-terminated string of unknown size needs a strlen
    //calculation. To avoid string literals selecting this constructor rather
    //than the templated array one, we protect with enable_if
    //(cf. https://stackoverflow.com/a/28243509).
    template <typename T>
    StrView(T cdata, typename std::enable_if<std::is_convertible<T, char const*>{}>::type* = nullptr) noexcept
      : m_data(cdata), m_size(cdata?std::strlen(cdata):0) {}

    //If needed in a constexpr context, this needs a special implementation
    //which is inefficient at run-time, hence we protect the constexpr version
    //with a special type signature (NB: the make(..) fct below might be more
    //convenient for string literals).
    struct constexpr_t {};
    constexpr StrView( constexpr_t, const char * cdata ) noexcept;

    //The following can be used to create constexpr strviews from string
    //literals (i.e. constexpr auto sv = StrView::make("foo")):
    template<size_type N>
    static constexpr StrView make( const char (&cdata)[N] ) noexcept
    {
      //NB: Defining directly in-class due to apple clang 11.
      return StrView( &cdata[0], ( cdata[N-1]=='\0' ? N-1 : N ) );
    }

    //Default constructing, or construction from nullptr/NullOpt, leads to a
    //StrView without value (which is NOT the same as an empty string).
    StrView() = default;
    constexpr StrView( NullOptType ) noexcept {}
    constexpr StrView( std::nullptr_t ) noexcept {}

    //Only mutable operation is move/copy assignment (so the objects can be
    //meaningfully kept in containers):
    ncconstexpr17 StrView& operator=( const StrView& o ) noexcept;
    ncconstexpr17 StrView& operator=( StrView&& o ) noexcept;

    //Check if the view has a value or not. If it does not, data() will return a
    //nullptr, and no methods should be used except the assignment operators.
    constexpr bool has_value() const noexcept;

    //Access data, and size (empty() is just a check for non-zero size). Note
    //that these three methods are safe to call even when has_value() is false,
    //(returning a nullptr, size 0, and empty()=true):
    constexpr const char * data() const noexcept { return m_data; }//Data. Warning: This is NOT necessarily null-terminated!!
    constexpr size_type size() const noexcept { return m_size; }
    constexpr bool empty() const noexcept { return m_size==0; }

    //Access data:
    constexpr const char * begin() const ncnoexceptndebug { return nc_assert_rv(m_data), m_data; }
    constexpr const char * end() const ncnoexceptndebug { return nc_assert_rv(m_data), m_data + m_size; }
    constexpr char operator[](size_type idx) const noexcept { return nc_assert_rv(m_data), m_data[idx]; }
    char at(size_type idx) const ncnoexceptndebug { return nc_assert_rv(m_data&&idx<m_size), m_data[idx]; }
    constexpr char front() const ncnoexceptndebug { return nc_assert_rv(m_data&&m_size>0), m_data[0]; }
    constexpr char back() const ncnoexceptndebug { return nc_assert_rv(m_data&&m_size>0), m_data[m_size - 1]; }

    //Create trimmed view with any leading or trailing whitespace removed:
    StrView trimmed() const;
    //left/right versions:
    StrView ltrimmed() const;
    StrView rtrimmed() const;

    //Substring or a full std::string:
    StrView substr( size_type first, size_type len = npos ) const noexcept;
    std::string to_string() const;

    void appendToString( std::string& ss ) const { ss.append( m_data, m_size ); }

    //To values:
    Optional<double> toDbl() const;
    Optional<int64_t> toInt() const { return toInt64(); }
    Optional<int32_t> toInt32() const;
    Optional<int64_t> toInt64() const;

    //String splitting:

    template< unsigned NPARTS_PREALLOC = 8>
    SmallVector<StrView,NPARTS_PREALLOC> split() const;//on normalised whitespace, no empty parts returned

    enum class SplitKeepEmpty { Yes, No };
    enum class SplitTrimParts { Yes, No };

    template< unsigned NPARTS_PREALLOC = 8,
              SplitKeepEmpty keep_empty = SplitKeepEmpty::Yes,
              SplitTrimParts trim_parts = SplitTrimParts::No >
    SmallVector<StrView,NPARTS_PREALLOC> split( char sep ) const;

    template< unsigned NPARTS_PREALLOC = 8,
              SplitKeepEmpty keep_empty = SplitKeepEmpty::Yes,
              SplitTrimParts trim_parts = SplitTrimParts::No >
    SmallVector<StrView,NPARTS_PREALLOC> split_any( const char * separators ) const;

    template< unsigned NPARTS_PREALLOC = 8,
              SplitKeepEmpty keep_empty = SplitKeepEmpty::Yes >
    SmallVector<StrView,NPARTS_PREALLOC> splitTrimmed( char sep ) const;

    template< unsigned NPARTS_PREALLOC = 8 >
    SmallVector<StrView,NPARTS_PREALLOC> splitTrimmedNoEmpty( char sep ) const;

    //startswith/endswith:
    bool startswith( StrView ) const;
    bool endswith( StrView ) const;
    bool startswith( char ) const;
    bool endswith( char ) const;

    //find/contains/contains_any functionality. Note that contains_any or
    //find_first_of are not as efficient as for std::string since we have to
    //implement it without the otherwise efficient std::strpbrk, which expects a
    //null terminated string.
    bool contains( char ) const noexcept;
    bool contains( StrView ) const;
    bool contains_any( const char * ) const noexcept;
    bool contains_any( StrView ) const noexcept;
    bool contains_only( StrView ) const;//only contains chars from given pattern

    size_type find( StrView ) const;
    size_type find( char ) const;
    size_type find_first_of( StrView ) const;
    size_type find_first_of( const char * ) const;

    enum class AllowTabs { Yes, No };
    enum class AllowNewLine { Yes, No };
    bool isSimpleASCII( AllowTabs = AllowTabs::Yes, AllowNewLine = AllowNewLine::Yes ) const;

    //Comparison operators (lexical sort) implemented with calls to std::strncmp:
    bool operator==(const StrView& o) const noexcept;
    bool operator!=(const StrView& o) const noexcept;
    bool operator<(const StrView& o) const noexcept;

    //Equivalent constexpr comparisons have implementations not suitable for
    //runtime evaluation, so we keep them in these separate methods:
    constexpr bool constexprEq(const StrView& o) const noexcept;
    constexpr bool constexprLessThan(const StrView& o) const noexcept;

  private:
    const char * m_data = nullptr;
    size_type m_size = 0;
  };

  //Miscellaneous:
  std::ostream& operator<<(std::ostream&, const StrView&);
  constexpr int64_t constexpr_strcmp( StrView sv1, StrView sv2 ) noexcept;

  class WordIterator {
  public:
    //Iterate over words in input, ignoring whitespace. End of text is reached
    //when next() returns an empty string.
    WordIterator( StrView text, StrView whitespace  = StrView::make(" \n\r\t") );
    StrView next();
  private:
    StrView m_text, m_ws;
  };

  ///////////////////////////////////////////////////////
  // Word wrapping:
  //
  // Normalise whitespace in text and stream word-wrapped Break up text on
  // whitespace and stream word-wrapped (inserting ' ' and '\n' between words
  // as appropriate). Note: Text is currently assumed to be ASCII (as UTF8
  // support would need us to handle multibyte chars).
  struct WordWrapCfg {
    //Column to wrap at:
    std::size_t colwidth = 80;
    //If set, assumed initial_offset.value() chars have already been streamed
    //on the first line. The first line will thus have less available
    //characters before wrapping takes place.
    Optional<std::size_t> initial_offset;
    //If all lines must be prefixed by a given string, provide it here:
    StrView prefix = StrView::make("");
    //If a single word is too long to fit in the available space, it will be
    //printed anyway unless the following is set:
    bool overflow_is_error = false;
    //If must always end with newline:
    bool ensure_final_newline = true;
    //Whitespace chars:
    StrView whitespace = StrView::make(" \n\r\t");
  };
  void streamWrappedText( std::ostream& os, StrView text, const WordWrapCfg& cfg = WordWrapCfg() );

}

#include "NCrystal/internal/utils/NCString.hh"

////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline constexpr StrView::StrView( const char * cdata, StrView::size_type len) noexcept : m_data(cdata), m_size(len) {}
  inline StrView::StrView( const std::string& stdstr ) noexcept : m_data(stdstr.c_str()), m_size(stdstr.size()) {}
  inline constexpr StrView::StrView( StrView::constexpr_t, const char * cdata ) noexcept : m_data(cdata), m_size(constexpr_strlen(cdata)) {}
  inline constexpr bool StrView::has_value() const noexcept { return m_data!=nullptr; }
  inline constexpr StrView::StrView( const StrView& o ) noexcept : m_data(o.m_data), m_size(o.m_size) {}
  inline constexpr StrView::StrView( StrView&& o ) noexcept : m_data(o.m_data), m_size(o.m_size) {}
  inline ncconstexpr17 StrView& StrView::operator=( const StrView& o ) noexcept { return m_data=o.m_data, m_size=o.m_size, *this; }
  inline ncconstexpr17 StrView& StrView::operator=( StrView&& o ) noexcept { return m_data=o.m_data, m_size=o.m_size, *this; }
  inline std::string StrView::to_string() const { nc_assert(m_data!=nullptr); return std::string(m_data,m_size); }

  inline constexpr int64_t constexpr_strcmp( StrView sv1, StrView sv2 ) noexcept
  {
    return constexpr_strcmp( sv1.data(), sv1.size(), sv2.data(), sv2.size() );
  }

  inline std::ostream& operator<<(std::ostream& os, const StrView& str)
  {
    os.write(str.data(),str.size());
    return os;
  }

  inline StrView StrView::substr( size_type first, size_type len ) const noexcept
  {
    //nc_assert(has_value());
    if ( first >= m_size || len == 0 )
      return StrView( m_data, 0 );//empty
    return StrView( m_data + first, std::min<size_type>(m_size-first,len) );
  }

  inline bool StrView::operator==(const StrView& o) const noexcept
  {
    return m_size == o.size() && 0 == std::strncmp( m_data, o.m_data, m_size );
  }

  inline bool StrView::operator!=(const StrView& o) const noexcept
  {
    return m_size != o.size() || 0 != std::strncmp( m_data, o.m_data, m_size );
  }

  inline bool StrView::operator<(const StrView& o) const noexcept
  {
    //Use std::strncmp to compare characters up to the length of the smallest
    //string. If not identical (cmp!=0), this already gives us the answer. If
    //identical (cmp==0), the shortest of the two strings is considered smaller.
    auto cr = std::strncmp( m_data, o.m_data,std::min<size_type>(m_size,o.m_size) );
    return ( cr == 0 ? m_size < o.m_size : cr < 0 );
  }

  inline constexpr bool StrView::constexprEq(const StrView& o) const noexcept
  {
    return constexpr_strcmp( *this, o );
  }

  inline constexpr bool StrView::constexprLessThan(const StrView& o) const noexcept
  {
    return constexpr_strcmp( m_data, m_size, o.m_data, o.m_size ) < 0;
  }

  inline bool StrView::contains( StrView pattern ) const
  {
    return find(pattern) != npos;
  }

  namespace detail {
    //strstr which does not assume null termination:
    const char * strstr_nonullterm( const char * cstr, std::size_t nc,
                                    const char * pattern, std::size_t np );
  }

  inline StrView::size_type StrView::find( StrView pattern ) const
  {
    nc_assert( has_value() );
    nc_assert( pattern.has_value() );
    nc_assert( !pattern.empty() );
    auto it = detail::strstr_nonullterm( m_data, m_size, pattern.m_data, pattern.m_size );
    return it ? std::distance( m_data, it ) : npos;
  }


  inline bool StrView::contains( char ch ) const noexcept
  {
    return m_size && std::memchr( m_data, ch, m_size ) != nullptr;
  }

  inline bool StrView::contains_any( const char * ss ) const noexcept
  {
    //Unfortunately we can not just call the otherwise efficient std::strpbrk,
    //since that expects a null terminated string.
    auto it = ss;
    while (*it)
      if (contains(*it++))
        return true;
    return false;
  }

  inline bool StrView::contains_any( StrView ss ) const noexcept
  {
    //Unfortunately we can not just call the otherwise efficient std::strpbrk,
    //since that expects a null terminated string.
    if ( m_size > 1024 && ss.size() > 1 && substr(0,128).contains_any(ss) )
      return true;//possible optimisation
    for ( char ch : ss )
      if (contains(ch))
        return true;
    return false;
  }

  inline StrView::size_type StrView::find( char ch ) const
  {
    auto it = (const char*)std::memchr( m_data, ch, m_size );
    return it ? std::distance(m_data,it): npos;
  }

  inline bool StrView::startswith( StrView o ) const
  {
    nc_assert(m_data);
    nc_assert(o.m_data);
    nc_assert(o.m_size);
    return o.m_size <= m_size && std::memcmp( m_data, o.m_data, o.m_size ) == 0;
  }

  inline bool StrView::endswith( StrView o ) const
  {
    nc_assert(m_data);
    nc_assert(o.m_data);
    nc_assert(o.m_size);
    return o.m_size <= m_size && std::memcmp( m_data + (m_size - o.m_size), o.m_data, o.m_size ) == 0;
  }

  inline bool StrView::startswith( char ch ) const
  {
    return m_size && m_data[0] == ch;
  }

  inline bool StrView::endswith( char ch ) const
  {
    return m_size && m_data[m_size-1] == ch;
  }

  inline Optional<double> StrView::toDbl() const { double val; if ( safe_str2dbl( *this, val ) ) return val; return NullOpt; }
  inline Optional<int32_t> StrView::toInt32() const { int32_t val; if ( safe_str2int( *this, val ) ) return val; return NullOpt; }
  inline Optional<int64_t> StrView::toInt64() const { int64_t val; if ( safe_str2int( *this, val ) ) return val; return NullOpt; }

  template< unsigned NP, StrView::SplitKeepEmpty KE, StrView::SplitTrimParts TP>
  inline SmallVector<StrView,NP> StrView::split( char sep ) const
  {
    constexpr bool do_keep_empty = (KE == SplitKeepEmpty::Yes);
    constexpr bool do_trim       = (TP == SplitTrimParts::Yes);
    SmallVector<StrView,NP> res;
    StrView sv(*this);
    nc_assert(sv.has_value());
    do {
      auto i = sv.find(sep);
      auto sv_part = sv.substr(0,i);
      if (do_trim)
        sv_part = sv_part.trimmed();
      if ( do_keep_empty || !sv_part.empty() )
        res.emplace_back(sv_part);
      sv = ( i==StrView::npos ? StrView() : sv.substr(i+1) );
    } while ( sv.has_value() );
    return res;
  }

  template< unsigned NP, StrView::SplitKeepEmpty KE, StrView::SplitTrimParts TP>
  inline SmallVector<StrView,NP> StrView::split_any( const char * separators ) const
  {
    nc_assert(separators&&*separators);
    constexpr bool do_keep_empty = (KE == SplitKeepEmpty::Yes);
    constexpr bool do_trim       = (TP == SplitTrimParts::Yes);
    SmallVector<StrView,NP> res;
    StrView sv(*this);
    nc_assert(sv.has_value());
    do {
      auto i = sv.find_first_of(separators);
      auto sv_part = sv.substr(0,i);
      if (do_trim)
        sv_part = sv_part.trimmed();
      if ( do_keep_empty || !sv_part.empty() )
        res.emplace_back(sv_part);
      sv = ( i==StrView::npos ? StrView() : sv.substr(i+1) );
    } while ( sv.has_value() );
    return res;
  }

  template<unsigned NP>
  inline SmallVector<StrView,NP> StrView::split() const
  {
    return split_any<NP,SplitKeepEmpty::No,SplitTrimParts::Yes>(" \t\n\r");
  }

  template<unsigned NP, StrView::SplitKeepEmpty KE>
  SmallVector<StrView,NP> StrView::splitTrimmed( char sep ) const
  {
    return split<NP,KE,SplitTrimParts::Yes>(sep);
  }

  template<unsigned NP>
  SmallVector<StrView,NP> StrView::splitTrimmedNoEmpty( char sep ) const
  {
    return split<NP,SplitKeepEmpty::No,SplitTrimParts::Yes>(sep);
  }

  inline bool StrView::isSimpleASCII( AllowTabs allow_tab, AllowNewLine allow_nl ) const
  {
    nc_assert(has_value());
    //NB: Keep in synch with NCString.icc
    //If not in 32..126, it must be a control char or char with high bit set
    //(e.g. extended ascii or utf-8 multibyte char). Of these, we only
    //(optionally) allow select whitespace.
    static_assert('\t'<32&&'\n'<32&&'\r'<32,"");
    auto it = m_data;
    auto itE = m_data+m_size;
    for ( ; it != itE; ++it ) {
      const char ch = *it;
      if ( ch > '\x1f' && ch < '\x7f' )//0x1f=31, 0x7f=127
        continue;
      if ( ch == '\t' && allow_tab == AllowTabs::No )
        return false;
      if ( ( ch == '\n' || ch == '\r' ) && allow_nl == AllowNewLine::No )
        return false;
    }
    return true;
  }

  inline WordIterator::WordIterator( StrView text, StrView whitespace )
    : m_text(text),
      m_ws(whitespace)
  {
    nc_assert_always(text.has_value());
    nc_assert_always(whitespace.has_value());
    nc_assert_always(!whitespace.empty());
  }

  //StrView-related ShortStr methods are defined here for dependency reasons:
  template <unsigned NMAX>
  inline ShortStr<NMAX>::ShortStr( StrView sv ) : ShortStr( sv.data(), sv.size() ) {}
  template <unsigned NMAX>
  inline constexpr StrView ShortStr<NMAX>::to_view() const noexcept { return StrView( data(), m_size ); }

}

#endif

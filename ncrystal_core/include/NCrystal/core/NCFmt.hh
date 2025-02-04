#ifndef NCrystal_Fmt_hh
#define NCrystal_Fmt_hh

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

// Utilities for formatting floating point values as strings. In particular
// efficient (malloc-free) std::ostream-compatbile formatting of floating point
// numbers. To support this, a short-string class is also provided.

#ifndef NCrystal_Defs_hh
#  include "NCrystal/core/NCDefs.hh"
#endif

namespace NCRYSTAL_NAMESPACE {

  class StrView;//Internal class (cf. NCStrView.hh)

  template <unsigned NMAX>
  class ShortStr {
    // Simple short-string class with a given maximum string length (including
    // null char), and no memory allocations.
  public:
    using size_type = unsigned;
    static_assert(NMAX>=1,"");
    static constexpr size_type max_strlen = static_cast<size_type>(NMAX - 1);
    static constexpr size_type bufsize = static_cast<size_type>(NMAX);

    ShortStr( const char * cstr );
    ShortStr( const char * cstr, unsigned str_len );
    ShortStr( StrView );

    //Empty string:
    ShortStr( NullOptType ) : m_size(0) { m_data[0] = '\0'; }
    ShortStr() : ShortStr( NullOpt ) {}

    //Special (and rather unsafe!!) constructor which allows direct access to
    //internal data members. This allows code to directly write to the ShortStr
    //buffers, but is not generally recommended:
    struct NCRYSTAL_API delayed_init_t {};
    ShortStr( delayed_init_t, char** c_handle, size_type** s_handle  )
    {
      *c_handle = &m_data[0];
      *s_handle = &m_size;
    }

    constexpr const char * data() const noexcept;
    constexpr size_type size() const noexcept;
    constexpr bool empty() const noexcept;
    constexpr StrView to_view() const noexcept;
    std::string to_string() const noexcept;

    ShortStr( const ShortStr& ) = default;
    ShortStr& operator=( const ShortStr& ) = default;
    ShortStr( ShortStr&& ) = default;
    ShortStr& operator=( ShortStr&& ) = default;

    //The following can be used to create constexpr objects from string literals
    //(i.e. constexpr auto sv = ShortStr::make("foo")): For now it is NOT
    //constexpr, due to limitations in C++11. This could be revisited if needed.
    template<size_type N>
    static ShortStr make( const char (&cdata)[N] )
    {
      return ShortStr( &cdata[0], ( cdata[N-1]=='\0' ? N-1 : N ) );
    }

  private:
    using TData = std::array<char,NMAX>;
    size_type m_size;
    TData m_data;
  };

  template <unsigned NMAX>
  std::ostream& operator<<(std::ostream&, const ShortStr<NMAX>&);

  ///////////////////////////////////////////////////////
  // String <-> FP conversions

  namespace detail {
    Optional<double> raw_str2dbl( const char *, std::size_t );
    Optional<std::int64_t> raw_str2int64( const char *, std::size_t );
  }

  //Buffer length needed to store any FP converted to string (NOT including any null terminator):
  namespace detail {
    NCRYSTAL_API inline constexpr unsigned ncconstexpr_log10ceil( unsigned val )
    {
      return val < 10u ? 1u : 1u + ncconstexpr_log10ceil( val / 10u );
    }

    template<typename T>
    inline constexpr unsigned ncconstexpr_max_size_fp_strbuf()
    {
      return 4u + std::numeric_limits<T>::max_digits10
        + ncconstexpr_max<unsigned>(2u, ncconstexpr_log10ceil(std::numeric_limits<T>::max_exponent10));
    }
    static constexpr unsigned max_size_strbuf_double = ncconstexpr_max_size_fp_strbuf<double>();//Normally 24
    static constexpr unsigned max_size_strbuf_float = ncconstexpr_max_size_fp_strbuf<float>();//Normally 15
  }

  //Short strings which are just large enough to hold non-lossy encoded floats:
  using ShortStrDbl = ShortStr< detail::max_size_strbuf_double + 1 >;
  using ShortStrFlt = ShortStr< detail::max_size_strbuf_float + 1 >;

  //Encode in shortest form which preserves value (unless fmtstr is set):
  ShortStrDbl dbl2shortstr( double value, const char * fmtstr = nullptr );

  //Utilities for printing double's to ostreams, either via a format string
  //(e.g. "%.3g"), or without. The latter prints the shortest/nicest string
  //which can be re-encoded as a double without loss of precision (in the sense
  //that if reencoded as a double it will get the original value):
  //
  //This means you can do "std::cout << fmt(myvalue)" (shortest loss-less) or
  //"std::cout << fmt(myvalue,"%.2f")" (specific formatting), without having to
  //bother with iomanip, stream globals (e.g. std::setprecision), or type unsafe
  //printf. The alias fmtg(myvalue) means the same as fmt(myvalue,"%g").
  struct NCRYSTAL_API detail_FmtDbl { double val; const char* fmtstr; };

  NCRYSTAL_API detail_FmtDbl fmt( double );
  NCRYSTAL_API detail_FmtDbl fmt( double, const char * fmtstr );
  NCRYSTAL_API detail_FmtDbl fmtg( double );//same as fmt(..,"%g")
  NCRYSTAL_API std::ostream& operator<<( std::ostream& , const detail_FmtDbl& );

  //Similar, but will use the detectSimpleRationalNumbers to detect and print
  //some simple rational numbers as fractions ("1/3" rather than "0.3333...").
  struct NCRYSTAL_API detail_FmtDblFrac { double val; const char* fmtstr; };
  std::ostream& operator<<( std::ostream& , const detail_FmtDblFrac& );

  NCRYSTAL_API detail_FmtDblFrac fmt_frac( double );
  NCRYSTAL_API detail_FmtDblFrac fmt_frac( double, const char * fmtstr );
  NCRYSTAL_API detail_FmtDblFrac fmtg_frac( double );//same as fmt_frac(..,"%g")

  //Find (A,B) so that value ~= A/B. This is only possible for a small subset
  //of typical values (non-negative integers below 1e9) and some simple fractions with
  //values in (0,1). When not possible, A=B=0.  The main use-case is to allow
  //for pretty-printing double values which were originally input as fractions
  //(e.g. we can print 1/3 instead of 0.33333333):
  std::pair<unsigned,unsigned> detectSimpleRationalNumbers(double value);
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  template <unsigned NMAX>
  inline constexpr const char * ShortStr<NMAX>::data() const noexcept { return &m_data[0]; }

  template <unsigned NMAX>
  inline constexpr typename ShortStr<NMAX>::size_type ShortStr<NMAX>::size() const noexcept { return m_size; }

  template <unsigned NMAX>
  inline constexpr bool ShortStr<NMAX>::empty() const noexcept { return m_size==0; }

  template <unsigned NMAX>
  inline std::string ShortStr<NMAX>::to_string() const noexcept { return std::string( data(), m_size ); }

  template <unsigned NMAX>
  inline ShortStr<NMAX>::ShortStr( const char * cstr ) : ShortStr( cstr, std::strlen(cstr) ) {}

  template <unsigned NMAX>
  inline ShortStr<NMAX>::ShortStr( const char * cstr, unsigned str_len )
    : m_size(str_len)
  {
    nc_assert( cstr );
    if ( str_len+1 > m_data.size() ) {
      std::ostringstream ss;
      ss << "String too long for ShortStr<"<<NMAX<<">: \"";
      ss.write(cstr,str_len);
      ss << '"';
      NCRYSTAL_THROW(BadInput,ss.str());
    }
    std::memcpy(&m_data[0],cstr,str_len);
    m_data[str_len] = '\0';
  }

  template <unsigned NMAX>
  inline std::ostream& operator<<(std::ostream& os, const ShortStr<NMAX>& str)
  {
    os.write(str.data(),str.size());
    return os;
  }

  inline detail_FmtDbl fmt( double val ) { return { val, nullptr }; }
  inline detail_FmtDbl fmt( double val, const char * fmtstr ) { return { val, fmtstr }; }
  inline detail_FmtDbl fmtg( double val ) { return { val, "%g" }; }
  inline std::ostream& operator<<( std::ostream& os , const detail_FmtDbl& fd ) { return os << dbl2shortstr( fd.val, fd.fmtstr ); }

  inline detail_FmtDblFrac fmt_frac( double val ) { return { val, nullptr }; }
  inline detail_FmtDblFrac fmt_frac( double val, const char * fmtstr ) { return { val, fmtstr }; }
  inline detail_FmtDblFrac fmtg_frac( double val ) { return { val, "%g" }; }
}


#endif

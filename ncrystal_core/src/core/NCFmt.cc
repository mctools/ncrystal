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

#include "NCrystal/core/NCFmt.hh"
namespace NC = NCrystal;

#include <streambuf>
#include <istream>

namespace NCRYSTAL_NAMESPACE {
  namespace detail {
    // Inspired by:
    // https://stackoverflow.com/questions/13059091/creating-an-input-stream-from-constant-memory/13059195#13059195
    struct nc_membuf : std::streambuf {
      //NB: base is memory buffer of given size, allocation/deallocation of
      //which should be handled by calling code.
      nc_membuf(char const* base, size_t size)
      {
        char* p = const_cast<char*>(base);
        this->setg(p, p, p + size);
      }
    };
    struct nc_imemstream : virtual nc_membuf, std::istream {
      nc_imemstream(char const* base, size_t size)
        : nc_membuf(base, size),
          std::istream(static_cast<std::streambuf*>(this))
      {
      }
    };

    namespace {
      static constexpr std::size_t denorm_str_maxsize = 64;

      Optional<double> raw_str2dbl_impl( const char * s_data,
                                         std::size_t s_size,
                                         bool allow_denorm_check = true );

      Optional<double> raw_str2dbl_trydenorm( const char * s_data,
                                              std::size_t s_size )
      {
        //NB: If std conversion fails and contains "e-3" it might be a denormal
        //    number (issue seen on some mac/clang CI failure). So we could
        //    replace the 3 with a 2, convert, and then multiply with 1e-100.
        //    This might not catch ALL denormal values
        //    (e.g. 0.000000000000001e-299), but it should hopefully catch those
        //    we ourselves have produced by our fmt() function, ensuring a
        //    conversion round-trip.
        nc_assert( s_size <= denorm_str_maxsize );
        const char * it = s_data;
        const char * itE = it + s_size;
        Optional<double> res;
        for ( ;it!=itE;++it ) {
          if ( *it != 'e' )
            continue;
          ++it;
          if ( it == itE || *it!='-'
               || std::next(it)==itE || *std::next(it)!='3' )
            return res;
          char tmpbuf[denorm_str_maxsize];
          std::memcpy(tmpbuf,s_data,s_size);
          auto idx = std::distance(s_data,it);
          nc_assert( tmpbuf[ idx ] == '-' );
          nc_assert( tmpbuf[ idx+1 ] == '3' );
          if ( s_size != static_cast<std::size_t>(idx+4) )
            return res;//has to be exactly two digits after the '3'
          tmpbuf[ std::distance(s_data,it)+1 ] = '2';
          std::string foo(tmpbuf,s_size);
          auto val100 = raw_str2dbl_impl( tmpbuf, s_size, false );
          if ( val100.has_value() && std::fabs( val100.value() ) < 1e-200 )
            res = val100.value() * 1e-100;
          break;
        }
        return res;
      }

#ifndef NDEBUG
      static int dummy = []()
      {
        auto test = []( std::string s, double vref )
        {
          auto v = raw_str2dbl_trydenorm( s.data(), s.size() );
          nc_assert_always( v.has_value() );
          nc_assert_always( v.value() == vref );
        };
        test("4.94066e-324",4.94066e-324);
        test("5e-320",5e-320);
        test("1e-308",1e-308);
        test("5e-324",5e-324);
        test("1.123456789123456789e-315",1.123456789123456789e-315);
        test("1e-301",1e-301);
        return 1;
      }();
#endif

      Optional<double> raw_str2dbl_impl( const char * s_data,
                                         std::size_t s_size,
                                         bool allow_denorm_check ) {
        //Using streams so we can specify the locale (TODO in c++17 we can
        //possibly use std::from_chars instead!). Using custom stream buffers to
        //reduce need for allocations:
        Optional<double> result;
        detail::nc_imemstream ss(s_data,s_size);
        ss.std::istream::imbue(std::locale::classic());
        double val;
        ss >> val;
        if ( ss.fail() ) {
          if ( allow_denorm_check && s_size <= denorm_str_maxsize )
            result = raw_str2dbl_trydenorm( s_data, s_size );
        } else if ( ss.eof() ) {
          result = val;
        }
        return result;
      }
    }

    Optional<double> raw_str2dbl( const char * s_data, std::size_t s_size ) {
      return raw_str2dbl_impl(s_data,s_size);
    }
    Optional<std::int64_t> raw_str2int64( const char * s_data,
                                          std::size_t s_size ) {
      //Using streams so we can specify the locale (TODO in c++17 we can possibly
      //use std::from_chars instead!). Using custom stream buffers to reduce need
      //for allocations:
      detail::nc_imemstream ss(s_data,s_size);
      ss.std::istream::imbue(std::locale::classic());
      std::int64_t val;
      ss >> val;
      if ( !ss.fail() && ss.eof() )
        return val;
      else
        return NullOpt;
    }

    //This one is not actually added to NCFMT.hh, since it is just needed in
    //NCString.cc, where it is fwd declared.
    Optional<std::uint64_t> raw_str2uint64( const char * s_data,
                                            std::size_t s_size ) {
      detail::nc_imemstream ss(s_data,s_size);
      ss.std::istream::imbue(std::locale::classic());
      std::uint64_t val;
      ss >> val;
      if ( !ss.fail() && ss.eof() ) {
        if ( val && s_data[0]=='-' )
          return NullOpt;//prevent negative values via underflow (but -0 is ok)
        return val;
      } else {
        return NullOpt;
      }
    }
  }
}

NC::ShortStrDbl NC::dbl2shortstr( double value, const char * fmtstr )
{
  //Consistent handling of special values (also means we won't have to worry
  //about nan!=nan issues below):
  if ( std::isnan( value ) )
    return ShortStrDbl::make("nan");
  if ( std::isinf( value ) )
    return value > 0.0 ? ShortStrDbl::make("inf") : ShortStrDbl::make("-inf");
  if ( value == 0.0 && !fmtstr )
    return ShortStrDbl::make("0");

  //printf's "%g" is really hard to replicate with ostreams, and the following
  //is anyway more efficient in terms of malloc's:
  //
  //TODO: When compiler supports it, and we move to C++17, we should likely
  //use std::to_char instead of these locale-dependent methods!

  //Use the delayed_init_t constructor of ShortStrDbl, to avoid a buffer copy
  //(this is complicated, but worth it in this case):
  char * buf;
  ShortStrDbl::size_type* buf_strlen_ptr;
  ShortStrDbl result( ShortStrDbl::delayed_init_t{}, &buf, &buf_strlen_ptr );
  ShortStrDbl::size_type& buf_strlen = *buf_strlen_ptr;

  auto dosnprintf = [&buf,&buf_strlen,value]( const char * fmtstr_ )
  {
    auto nwritten_strlen = std::snprintf(buf,ShortStrDbl::bufsize,fmtstr_,value);
    nc_assert( nwritten_strlen > 0 && std::int64_t(nwritten_strlen+1) <= std::int64_t(ShortStrDbl::bufsize) );
    nc_assert( buf[nwritten_strlen] == '\0' );
    buf_strlen = nwritten_strlen;
    nc_assert( (int)std::strlen(buf) == (int)buf_strlen );
    //Protect against non-C locale creating commas instead of dots for decimal
    //separators:
    auto pcomma = std::memchr(buf,',',buf_strlen);
    if ( pcomma )
      *(char*)pcomma = '.';
  };

  if ( fmtstr ) {
    //honour request for specific fmt string:
    dosnprintf(fmtstr);
  } else {
    //No specific format string requested, aim for lossless and "nice".

    constexpr const char* fmt_patterns[] = {
      "%.1g",  "%.2g",  "%.3g",  "%.4g",  "%.5g",
      "%.6g",  "%.7g",  "%.8g",  "%.9g", "%.10g",
      "%.11g", "%.12g", "%.13g", "%.14g", "%.15g",
      "%.16g", "%.17g", "%.18g", "%.19g", "%.20g",
      "%.21g", "%.22g", "%.23g", "%.24g", "%.25g",
      "%.26g", "%.27g", "%.28g", "%.29g", "%.30g",
    };

    unsigned nstrip = 0;
    double testval = value;
    while ( std::fabs(testval) < std::numeric_limits<double>::min() ) {
      ++nstrip;
      testval *= 10;
    }
    if ( nstrip > 0 )
      --nstrip;

    constexpr int idxfmt_std = std::numeric_limits<double>::max_digits10-1;
    const int idxfmt = std::max<int>(1,idxfmt_std-nstrip);
    const int idxfmtm2 = std::max<int>(0,(idxfmt_std-2)-nstrip);

    const char * fmt_full = fmt_patterns[idxfmt];
    const char * fmt_full_m2 = fmt_patterns[idxfmtm2];

    //First try with a bit less than maxdigits10, which can prevent some
    //originally nice looking values being presented in a less nice looking
    //fashion (e.g. 0.1 -> 0.09999999999999999).
    dosnprintf(fmt_full_m2);
    auto backconv = detail::raw_str2dbl( buf, buf_strlen );
    if ( !backconv.has_value() || backconv.value() != value )
      dosnprintf(fmt_full);
#ifndef NDEBUG
    {
      //loss less check:
      auto backconv2 = detail::raw_str2dbl( buf, buf_strlen );
      nc_assert( backconv2.has_value() && backconv2.value() == value );
    }
#endif
  }
  return result;
}

std::pair<unsigned,unsigned> NC::detectSimpleRationalNumbers(double value)
{
  if (value<=0.0)
    return { 0, ( value == 0.0 ? 1 : 0 ) };
  if (value>=1.0) {
    if (value==1.0)
      return { 1, 1 };
    double intpart;
    if ( std::modf(value, &intpart) != 0.0
         || intpart >= (double)std::numeric_limits<unsigned>::max() )
      return {0,0};
    return { static_cast<unsigned>(intpart), unsigned(1) };
  }
  //ok, value is in (0,1).
  static std::map<uint64_t,std::pair<unsigned,unsigned>> s_cache = []()
  {
    std::map<uint64_t,std::pair<unsigned,unsigned>> result;
    for (unsigned b = 2; b <= 64; ++b) {
      for (unsigned a = 1; a < b; ++a) {
        uint64_t key = static_cast<uint64_t>((double(a)/b)*1e18);
        //NB: key fits in uint64_t since value of a/b is in (0,1)
        auto it = result.find(key);
        if (it==result.end())
          result[key] = { a,b };//only insert if not there (because if a1/b1 =
                                //a2/b2, we prefer the smallest b-value).
      }
    }
    return result;
  }();
  uint64_t key = static_cast<uint64_t>(value*1e18);
  auto it = s_cache.find(key);
  if ( it == s_cache.end() )
    return {0,0};
  return it->second;
}

namespace NCRYSTAL_NAMESPACE {
  std::ostream& operator<<( std::ostream& os , const detail_FmtDblFrac& fd ) {
    auto ab = detectSimpleRationalNumbers(fd.val);
    if ( ab.second == 1 )
      os << ab.first;//integer
    else if ( ab.second == 0 )
      os << dbl2shortstr( fd.val, fd.fmtstr );
    else
      os << ab.first<<"/"<<ab.second;//fraction detected
    return os;
  }
}

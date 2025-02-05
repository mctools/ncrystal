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

#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/core/NCFmt.hh"
#include <istream>
#include <iomanip>
namespace NC = NCrystal;

std::string NC::displayCharSafeQuoted( char ch_raw, char quote_char )
{
  signed char ch = *reinterpret_cast<signed char*>(&ch_raw);
  std::ostringstream os;
  os << quote_char;
  if ( ch == quote_char )
    os << '\\';//e.g. """ -> "\""
  if ( ch >= ' ' && ch < '\x7f' ) {//' ' = 32, 0x7f=127
    os << ch;
  } else {
    os << "\\x";
    int ch_val(ch);
    if ( ch_val < 10 )
      os << '0';
    os << ch_val;
  }
  os << quote_char;
  return os.str();
}

void NC::trim( std::string& input )
{
  char * cB = &input[0];
  const char * cE = cB + input.size();
  const char * c = cB;
  while (c!=cE&&(*c==' '||*c=='\t'||*c=='\n'||*c=='\r'))
    ++c;
  if (c==cE) {
    input.clear();
    return;
  }
  const char * clast = cE-1;
  while (clast>c&&(*clast==' '||*clast=='\t'||*clast=='\n'||*clast=='\r'))
    --clast;
  nc_assert(clast>=c);
  std::size_t len = (clast + 1) - c;
  if (c!=cB) {
    for (std::size_t i = 0; i < len; ++i)
      *(cB + i) = *(c + i);
  }
  if (len<input.size())
    input.resize(len);
}

NC::VectS NC::split2( const std::string& input, std::size_t maxsplit,  char sep )
{
  VectS output;
  split(output,input,maxsplit,sep);
  return output;
}

void NC::split(NC::VectS& output, const std::string& input, std::size_t maxsplit, char sep )
{
  bool keep_empty = (sep!=0);
  //Split input on whitespace (" \t\n") with minimal memory allocations and
  //copying, discarding whitespace which is either leading and trailing or
  //repeating.
  nc_assert(bool(sep!='\0')==bool(sep));
  output.clear();

  if (input.empty()) {
    if (sep)
      output.emplace_back();
    return;
  }

  if (sep&&input[0]==sep)
    output.emplace_back();
  const char * c = &input[0];
  const char * cE = c + input.size();
  const char * partbegin = 0;
  while (true) {
    if ( maxsplit && output.size() == maxsplit ) {
      output.emplace_back(c);
      return;
    }
    if (c==cE||(sep?*c==sep:(*c==' '||*c=='\t'||*c=='\n'||*c=='\r'))) {
      if (partbegin) {
        if (keep_empty||c>partbegin) {
          output.emplace_back(partbegin,c-partbegin);
        }
        partbegin= keep_empty ? c+1 : 0;
      }
      if (c==cE)
        return;
    } else if (!partbegin) {
      partbegin = c;
    }
    ++c;
  }
}

bool NC::startswith(const std::string& str, const std::string& substr)
{
  return str.size()>=substr.size() && str.compare(0, substr.size(), substr) == 0;
}

bool NC::endswith(const std::string& str, const std::string& substr)
{
  return str.size()>=substr.size() && str.compare(str.size()-substr.size(), substr.size(), substr) == 0;
}

bool NC::contains_any(const std::string& haystack, const std::string& needles)
{
  const char * cB = &needles[0];
  const char * cE = cB + needles.size();
  for (const char * c = cB; c!=cE; ++c)
    if (haystack.find(*c) != std::string::npos)
      return true;
  return false;
}

double NC::str2dbl( StrView sv, const char * errmsg)
{
  double result;
  if ( !safe_str2dbl(sv, result ) )
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not a double")<<": \""<<sv<<"\"");
  return result;
}

int NC::str2int( StrView sv, const char * errmsg)
{
  int result;
  if ( !safe_str2int(sv, result ) )
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not an integer")<<": \""<<sv<<"\"");
  return result;
}

int32_t NC::str2int32( StrView sv, const char * errmsg)
{
  int32_t result;
  if ( !safe_str2int(sv, result ) )
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not an integer")<<": \""<<sv<<"\"");
  return result;
}

int64_t NC::str2int64( StrView sv, const char * errmsg)
{
  int64_t result;
  if ( !safe_str2int(sv, result ) )
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not an integer")<<": \""<<sv<<"\"");
  return result;
}

bool NC::safe_str2dbl( StrView s, double& result )
{
  const auto n_size = s.size();
  if ( n_size==0 || isWhiteSpace( s.front() ) || isWhiteSpace( s.back() ) )
    return false;

  {
    auto try_conv = detail::raw_str2dbl( s.data(), n_size );
    if ( try_conv.has_value() ) {
      result = try_conv.value();
      return true;
    }
  }

  //Did not work, check for "+-inf"/"+-INF"/"nan"/"NaN","NAN" for cross-platform
  //consistency:
  if ( n_size == 3 ) {
    if ( isOneOf(s,"inf","INF") ) {
      result = kInfinity;
      return true;
    }
    if ( isOneOf(s,"nan","NAN","NaN") ) {
      static_assert(std::numeric_limits<double>::has_quiet_NaN,"");
      result = std::numeric_limits<double>::quiet_NaN();
      return true;
    }
  }
  if ( n_size == 4 ) {
    if ( isOneOf(s,"+inf","+INF") ) {
      result = kInfinity;
      return true;
    }
    if ( isOneOf(s,"-inf","-INF") ) {
      result = -kInfinity;
      return true;
    }
  }

  //Syntax error:
  return false;
}

bool NC::safe_str2int( StrView s, int32_t& result )
{
  //For robustness/simplicity simply parse as 64 bit value and check range.
  int64_t val64;
  static constexpr int64_t range_min = static_cast<int64_t>( std::numeric_limits<int32_t>::lowest() );
  static constexpr int64_t range_max = static_cast<int64_t>( std::numeric_limits<int32_t>::max() );
  if ( !safe_str2int(s,val64) || val64 > range_max || val64 < range_min )
    return false;
  result = static_cast<int32_t>(val64);
  return true;
}

bool NC::safe_str2int( StrView s, int64_t& result )
{
  const auto n_size = s.size();
  if ( n_size==0 || isWhiteSpace( s.front() ) || isWhiteSpace( s.back() ) )
    return false;
  auto try_conv = detail::raw_str2int64( s.data(), n_size );
  if ( try_conv.has_value() ) {
    result = try_conv.value();
    return true;
  } else {
    return false;
  }
}

bool NC::contains_only(const std::string& haystack, const std::string& needles)
{
  const char * c = &haystack[0];
  const char * cE = c + haystack.size();
  for (;c!=cE;++c) {
    if (!contains(needles,*c))
      return false;
  }
  return true;
}

void NC::strreplace(std::string& str, const std::string& oldtxt, const std::string& newtxt)
{
  if( oldtxt.empty() )
    return;
  std::size_t offset = 0;
  while( ( offset = str.find(oldtxt, offset) ) != std::string::npos ) {
    str.replace( offset, oldtxt.length(), newtxt );
    offset += newtxt.size();
  }
}

std::string NC::joinstr( const Span<const StrView>& parts, StrView sep)
{
  const std::size_t n = parts.size();
  if ( n < 2 )
    return n ? parts.at(0).to_string() : std::string();

  std::size_t newsize( (n-1)*sep.size() );
  for (const auto& p : parts)
    newsize += p.size();
  std::string tmp;
  tmp.reserve(newsize);
  tmp.append(parts.front().data(),parts.front().size());
  for (std::size_t i = 1; i < n; ++i) {
    tmp.append( sep.data(), sep.size() );
    tmp.append( parts[i].data(), parts[i].size() );
  }
  nc_assert(tmp.size()==newsize);
  return tmp;
}

unsigned NC::countTrailingDigits( const std::string& ss )
{
  auto nn = ss.size();
  nc_assert_always(static_cast<uint64_t>(nn)<static_cast<uint64_t>(std::numeric_limits<int>::max()));
  int n = static_cast<int>(nn);

  int nTrailingDigits(0);
  while ( nTrailingDigits < n ) {
    int i = n-(nTrailingDigits+1);
    nc_assert(i>=0);
    auto c = ss.at( i );
    if ( 'c' < '0' || c > '9' )
      break;
    ++nTrailingDigits;
  }
  return static_cast<unsigned>(nTrailingDigits);
}

NC::PairSS NC::decomposeStrWithTrailingDigits( const std::string& ss )
{
  unsigned nTrailingDigits = countTrailingDigits(ss);
  if (nTrailingDigits==0)
    return {ss,std::string()};
  auto nn = ss.size() - nTrailingDigits;
  return { ss.substr(0,nn), ss.substr(nn) };
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    struct GetEnvResult {
      std::string varname;
      const char * val;
    };
    GetEnvResult raw_getenv( std::string& v )
    {
      nc_assert(!v.empty());
      nc_assert(!startswith(v,"NCRYSTAL_"));//common mistake
      GetEnvResult res;
#if defined(NCRYSTAL_NAMESPACE_PROTECTION) && defined( NCRYSTAL_NAMESPACED_ENVVARS )
      res.varname.reserve( 64 );
      res.varname = "NCRYSTAL";
      res.varname += upperCase(ncrystal_xstr(NCRYSTAL_NAMESPACE_PROTECTION));
      res.varname += '_';
#else
      res.varname = "NCRYSTAL_";
#endif
      res.varname += v;
      res.val = std::getenv( res.varname.c_str() );
      return res;
    }
  }
}

//Common access to environment variables (unset and empty vars both return
//empty strings / 0.0 / 0, depending on which function is used):
std::string NC::ncgetenv(std::string v, std::string defval)
{
  auto res = raw_getenv(v);
  return res.val ? std::string(res.val) : defval;
}

double NC::ncgetenv_dbl(std::string v, double defval )
{
  auto res = raw_getenv(v);
  if (!res.val)
    return defval;
  double result;
  if ( !safe_str2dbl(res.val, result ) )
    NCRYSTAL_THROW2(BadInput,"Invalid value of environment variable "
                    <<res.varname
                    <<" (expected a floating point number but got \""
                    <<res.val<<"\").");
  return result;
}

int NC::ncgetenv_int(std::string v, int defval )
{
  auto res = raw_getenv(v);
  if (!res.val)
    return defval;
  int result;
  if ( !safe_str2int(res.val, result ) )
    NCRYSTAL_THROW2(BadInput,"Invalid value of environment variable "<<res.varname
                    <<" (expected an integral number but got \""<<res.val<<"\").");
  return result;
}

std::int64_t NC::ncgetenv_int64(std::string v, std::int64_t defval )
{
  auto res = raw_getenv(v);
  if (!res.val)
    return defval;
  std::int64_t result;
  if ( !safe_str2int(res.val, result ) )
    NCRYSTAL_THROW2(BadInput,"Invalid value of environment variable "<<res.varname
                    <<" (expected an integral number but got \""<<res.val<<"\").");
  return result;
}


bool NC::ncgetenv_bool(std::string v)
{
  auto res = raw_getenv(v);
  if (!res.val)
    return false;
  std::string evs(res.val);
  if (evs.size()==1) {
    if (evs[0]=='0')
      return false;
    if (evs[0]=='1')
      return true;
  }
  NCRYSTAL_THROW2(BadInput,"Invalid value of environment variable "<<res.varname
                  <<" (expected a Boolean value, \"0\" or \"1\", but got \""<<evs<<"\").");
}

std::string NC::bytes2hexstr(const std::vector<uint8_t>& v) {
  const char hexchars[] = "0123456789abcdef";
  std::ostringstream ss;
  for (auto e : v)
    ss << hexchars[ e >> 4 ] << hexchars[ e & 0x0F ];
  return ss.str();
}

std::vector<uint8_t> NC::hexstr2bytes(const std::string& v) {
  std::vector<uint8_t> res;
  auto hex2val = [](unsigned c) {
    unsigned val(99);
    constexpr unsigned val_a = static_cast<unsigned>('a');
    constexpr unsigned val_A = static_cast<unsigned>('A');
    constexpr unsigned val_0 = static_cast<unsigned>('0');
    static_assert( val_0 < val_A && val_A < val_a,"" );
    if ( c >= val_a )
      val = 10 + ( c - val_a );
    else if ( c >= val_A)
      val = 10 + ( c - val_A );
    else if ( c >= val_0)
      val = ( c - val_0 );
    if ( val > 15 )
      NCRYSTAL_THROW2(BadInput, "Invalid character encountered in hex string: "<<c<<" (numeric value)");
    return val;
  };

  auto it = v.begin();
  auto itE = v.end();
  const auto vlen = v.size();
  if ( vlen%2 != 0 ) {
    //0 pad, e.g. interpret (0x)"A" as (0x)"0A".
    res.reserve((vlen+1)/2);
    res.push_back(hex2val(*it++));
  } else {
    res.reserve( vlen/2 );
  }
  while ( it != itE ) {
    auto c = *it++;
    res.push_back ( hex2val(*it++) + 16* hex2val(c) );
  }
  return res;
}

void NC::streamJSON( std::ostream& os, double val )
{
  if ( std::isnan( val ) )
    NCRYSTAL_THROW(CalcError,"Can not represent not-a-number (NaN) values in JSON format!");

  if ( std::isinf(val) ) {
    //infinity is not supported in json, but at least the python json modules
    //parses an out of bounds value as such:
    os << ( val > 0 ? "1.0e99999" : "-1.0e99999" );
    return;
  }
  if ( val == 0.0 ) {
    os << "0.0";//avoid -0, 0.0, etc.
    return;
  }
  //Encode as string, but make sure that we present floating point numbers in a
  //type-preserving way ("5.0", not "5") so we won't end up mapping C++ floating
  //point types to integer types once the json string is decoded in e.g. Python.
  auto sstr = dbl2shortstr(val);
  if ( sstr.to_view().toInt().has_value() )
    os << fmt( val, "%.1f" );
  else
    os << sstr;
}

void NC::streamJSON( std::ostream& os, StrView sv )
{
  //According to answer at:
  //https://stackoverflow.com/questions/3020094/how-should-i-escape-strings-in-json:
  //  "JSON is pretty liberal: The only characters you must escape are \, ", and
  //   control codes (anything less than U+0020)."
  //In addition we use the special forms for "\n" "\t" "\r" for more
  //readable/short results.

  os << '"';
  for ( auto ch : sv ) {
    //Code below written to be valid and hopefully non-warning producing for
    //both signed and unsigned char.
    switch ( ch ) {
    //Special cases first:
    case '\x00' : os << '"'; return;//null-char in data, immediately end!
    case '"'    : os << "\\\""; continue;
    case '\\'   : os << "\\\\"; continue;
    case '\n'   : os << "\\n"; continue;
    case '\r'   : os << "\\r"; continue;
    case '\t'   : os << "\\t"; continue;
    default:
      if ( ch < '\x20' && ch >= '\x01' ) {
        //Control character. Must encode hex value using \u00XX format. We use
        //snprintf rather than iomanip to avoid altering the state of the
        //stream, and to hopefully keep malloc's low:
        char buf[8];
        std::snprintf(buf, sizeof(buf)-1, "%04x", (unsigned char)(ch));
        os << "\\u"<<buf;
      } else {
        //Not control char or special case => simply pass through.
        os << ch;
      }
    }
  }
  os << '"';
}

NC::Optional<std::string> NC::findForbiddenChar( const StrView& teststr,
                                                 const StrView& forbidden_chars,
                                                 ExtraForbidOpt extra )
{
  Optional<char> badchar;
  if ( extra == ExtraForbidOpt::RequireSimpleASCII ) {
    for ( auto ch : teststr ) {
      if (!isSimpleASCII(ch,AllowTabs::Yes,AllowNewLine::Yes)) {
        badchar = ch;
        break;
      }
    }
  }
  if ( !badchar.has_value() && forbidden_chars.has_value() ) {
    auto pos = teststr.find_first_of(forbidden_chars);
    if ( pos != std::string::npos )
      badchar = teststr.at(pos);
  }
  if ( !badchar.has_value() )
    return NullOpt;
  return displayCharSafeQuoted(badchar.value());
}

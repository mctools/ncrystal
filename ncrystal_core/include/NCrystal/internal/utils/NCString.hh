#ifndef NCrystal_String_hh
#define NCrystal_String_hh

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

//String-related utilities

#include "NCrystal/internal/utils/NCSpan.hh"

namespace NCRYSTAL_NAMESPACE {

  class StrView;

  //All bytes must be in range 32..126 (plus optionally new-lines and tabs).
  enum class AllowTabs { Yes, No };
  enum class AllowNewLine { Yes, No };
  constexpr bool isSimpleASCII(char ch, AllowTabs = AllowTabs::Yes, AllowNewLine = AllowNewLine::Yes) noexcept;
  ncconstexpr17 bool isSimpleASCII(const char* str_begin, std::size_t len, AllowTabs = AllowTabs::Yes, AllowNewLine = AllowNewLine::Yes) noexcept;
  bool isSimpleASCII(const std::string&, AllowTabs = AllowTabs::Yes, AllowNewLine = AllowNewLine::Yes) noexcept;
  bool isSimpleASCII(const char*, AllowTabs = AllowTabs::Yes, AllowNewLine = AllowNewLine::Yes) noexcept;

  //Strip excess whitespace (" \t\r\n") from both ends of string:
  void trim( std::string& );
  std::string trim2( std::string );
  VectS trimEntries( VectS&& );//trim all entries

  //Check for whitespace (" \t\r\n"):
  constexpr bool isWhiteSpace( const char ) noexcept;

  //Split input string on separator (default sep=0 means splitting on general
  //whitespace - " \t\r\n"). Empty parts are only kept when sep!=0 (similar to
  //pythons str.split()). Finally, maxsplit can be used to limit the number of
  //splittings performed:
  VectS split2( const std::string& input,
                std::size_t maxsplit = 0,
                char sep = 0);

  //Backwards compatible version. Results are placed in output vector, which is
  //first cleared of existing contents. Eventually we will remove this version
  //and rename split2->split:
  void split(VectS& output,
             const std::string& input,
             std::size_t maxsplit = 0,
             char sep = 0 );

  //Substrings at edges:
  bool startswith(const std::string& str, const std::string& substr);
  bool endswith(const std::string& str, const std::string& substr);
  bool startswith(const std::string& str, char );
  bool endswith(const std::string& str, char );

  //Check if given char or substring (needle) is present in a string (haystack):
  bool contains(const std::string& haystack, char needle );
  bool contains(const std::string& haystack, const std::string& needle);
  //Check if any of the chars in "needles" is present in the string (haystack):
  bool contains_any(const std::string& haystack, const std::string& needles);
  //Check if "haystack" consists entirely of chars from string needles:
  bool contains_only(const std::string& haystack, const std::string& needles);

  //Change string casing:
  std::string lowerCase( std::string );
  std::string upperCase( std::string );

  //String justification:
  std::string str_rjust( std::string, std::size_t n );
  std::string str_ljust( std::string, std::size_t n );

  //Check if letter (defined as containing only a-zA-Z):
  constexpr bool isAlpha( const char ) noexcept;

  //Check if alphanumeric (defined as containing only a-zA-Z0-9):
  constexpr bool isAlphaNumeric( const char ) noexcept;
  bool isAlphaNumeric( const std::string& ) noexcept;

  //Convert strings to numbers. In case of problems, a BadInput exception will
  //be thrown (provide err to modify the message in that exception):
  NCRYSTAL_API double str2dbl( StrView, const char * errmsg = 0);//marked NCRYSTAL_API since used in custom physics example
  int str2int( StrView, const char * errmsg = 0);
  int32_t str2int32( StrView, const char * errmsg = 0);
  int64_t str2int64( StrView, const char * errmsg = 0);

  //Versions which don't throw:
  bool safe_str2dbl( StrView, double& result );
  bool safe_str2int( StrView, int32_t& result );
  bool safe_str2int( StrView, int64_t& result );

  //How many digits does string end with (e.g. 1 in "H1", 0 in "H1a", 3 in "Bla123".
  unsigned countTrailingDigits( const std::string& ss );

  //"Bla123" => ("Bla","123"). Special cases: "Bla" -> ("Bla","") "Bla012" -> ("Bla","012")
  PairSS decomposeStrWithTrailingDigits( const std::string& ss );

  //Replace all occurances of oldtxt in str with newtxt:
  void strreplace(std::string& str, const std::string& oldtxt, const std::string& newtxt);

  //Join ["a","bb","123"] -> "a bb 123":
  std::string joinstr(const VectS& parts );//join by single space
  std::string joinstr(const Span<const StrView>& parts );//join by single space
  std::string joinstr(const Span<const StrView>& parts, StrView sep );
  std::string joinstr(const VectS& parts, StrView sep );

  //Convert values to/from hex strings:
  std::string bytes2hexstr(const std::vector<uint8_t>& v);
  std::vector<uint8_t> hexstr2bytes(const std::string& v);

  //Common access to environment variables - will always be prefixed with
  //NCRYSTAL_. Unset variables means that the default values will be
  //returned. The _dbl/_int versions throws BadInput exceptions in case of
  //problems:
  //NB: Assuming no-one calls setenv/putenv/unsetenv, concurrent calls to these
  //functions are safe (since C++11).
  std::string ncgetenv(std::string, std::string defval = std::string() );
  double ncgetenv_dbl(std::string, double defval = 0.0);
  int ncgetenv_int(std::string, int defval = 0 );
  std::int64_t ncgetenv_int64(std::string, int64_t defval = 0 );
  bool ncgetenv_bool(std::string);//if set to 1 -> true, 0/unset -> false (otherwise exception).

  //Find forbidden characters. Either from a given list and/or by looking for non-ASCII characters.
  //If any are found, return a string representation suitable for printing in error messages.
  enum class ExtraForbidOpt{ RequireSimpleASCII, None };
  Optional<std::string> findForbiddenChar( const StrView& teststr,
                                           const StrView& forbidden_chars,
                                           ExtraForbidOpt extra = ExtraForbidOpt::None );
  //Display char like:>>"a"<<, >>"$"<<, >>"\""<<, >>"\x09"<<, ...
  std::string displayCharSafeQuoted( char ch, char quote_char = '"' );

  ///////////////////////////////////////////////////////
  // A few utilities only intended for constexpr
  // code (inefficient for runtime usage):

  constexpr std::size_t constexpr_strlen( const char* ) noexcept;
  constexpr bool constexpr_cstrequal( const char * c1, const char * c2 ) noexcept;

  //Non-null terminated string:
  constexpr int64_t constexpr_strcmp( const char* c1, std::size_t l1,
                                      const char* c2, std::size_t l2 ) noexcept;

  ///////////////////////////////////////////////////////
  // JSON encoding data:

  //Strings (ascii or utf-8) as properly escaped JSON strings:
  void streamJSON( std::ostream&, StrView );
  void streamJSON( std::ostream&, const char * );
  void streamJSON( std::ostream&, const std::string& );
  //Floating point as JSON numbers (always keeping double precision and FP-type
  //identity, e.g. "5.0" not "5"):
  void streamJSON( std::ostream&, float );
  void streamJSON( std::ostream&, double );
  //Integers as JSON integer numbers (but forbid 8-bit ints and chars to avoid type ambiguity):
  template<class T, typename std::enable_if<std::is_integral<T>::value>::type* = nullptr>
  void streamJSON( std::ostream&, T );
  void streamJSON( std::ostream&, char ) = delete;
  void streamJSON( std::ostream&, uint8_t ) = delete;
  void streamJSON( std::ostream&, int8_t ) = delete;
  //Bools become true / false in json:
  void streamJSON( std::ostream&, bool );
  struct json_null_t{};
  void streamJSON( std::ostream&, json_null_t );//to write out JSON null
  //Containers as JSON arrays:
  template<class TContainer, typename T = typename TContainer::value_type>
  void streamJSON( std::ostream&, const TContainer& );
  template<class T1, class T2>
  void streamJSON( std::ostream&, const std::pair<T1,T2>& );
  //Entries in dictionary "\"key\":<value>".
  enum class JSONDictPos { FIRST, LAST, OTHER };
  template<class T>
  inline void streamJSONDictEntry( std::ostream&,
                                   const char * key,
                                   const T& value,
                                   JSONDictPos pos = JSONDictPos::OTHER);
}

#include "NCrystal/internal/utils/NCStrView.hh"

////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  inline constexpr std::size_t constexpr_strlen( const char* cstr ) noexcept {
    //C++11 compatible => needs single recursive statement.
    return *cstr ? 1 + constexpr_strlen(cstr + 1) : 0;
  }

  inline constexpr bool constexpr_cstrequal( const char * c1, const char * c2 ) noexcept {
    return ( *c1 != *c2
             ? false
             : ( *c1 ? constexpr_cstrequal(c1+1,c2+1) : true ) );
  }

  inline constexpr int64_t constexpr_strcmp( const char* c1, std::size_t l1, const char* c2, std::size_t l2 ) noexcept
  {
    //C++11 compatible => needs single recursive statement.
    return ( l1==0
             ? ( l2==0 ? 0 : -1 )
             : ( l2==0
                 ? 1
                 : ( c1[0]==c2[0] ? constexpr_strcmp(c1+1,l1-1,c2+1,l2-1) : c1[0]-c2[0] )
                 )
             );
  }

  inline bool isSimpleASCII(const char * cstr, AllowTabs allow_tab, AllowNewLine allow_newline) noexcept
  {
    for ( auto c = cstr; *c; ++c ) {
      if ( !isSimpleASCII(*c, allow_tab, allow_newline) )
        return false;
    }
    return true;
  }

  inline ncconstexpr17 bool isSimpleASCII(const char* str_begin, std::size_t len, AllowTabs allow_tab, AllowNewLine allow_newline) noexcept
  {
    auto c = str_begin;
    auto cE = c+len;
    for ( ; c != cE; ++c )
      if ( !isSimpleASCII(*c,allow_tab,allow_newline) )
        return false;
    return true;
  }

  inline constexpr bool isSimpleASCII(char ch, AllowTabs allow_tab, AllowNewLine allow_newline) noexcept
  {
    //If not in 32..126, it must be a control char or char with high bit set
    //(e.g. extended ascii or utf-8 multibyte char). Of these, we only
    //(optionally) allow select whitespace.
    static_assert('\t'<32&&'\n'<32&&'\r'<32,"");
    return ( ch > '\x1f' && ch < '\x7f' )//0x1f=31, 0x7f=127
      || ((allow_tab==AllowTabs::Yes)&&ch=='\t')
      || ((allow_newline==AllowNewLine::Yes)&&(ch=='\n'||ch=='\r'));
  }

  inline bool isSimpleASCII(const std::string& s, AllowTabs allow_tab, AllowNewLine allow_newline) noexcept
  {
    auto c = s.c_str();
    auto cE = c + s.size();
    for ( ; c != cE; ++c )
      if ( !isSimpleASCII(*c,allow_tab,allow_newline) )
        return false;
    return true;
  }

  inline constexpr bool isWhiteSpace(const char c) noexcept
  {
    return c <= ' ' && c >= '\t' && ( c==' ' || c=='\n' || c=='\t' || c=='\r' );
  }

  inline bool contains(const std::string& haystack, char needle )
  {
    return haystack.find(needle) != std::string::npos;
  }

  inline bool contains(const std::string& haystack, const std::string& needle)
  {
    return haystack.find(needle) != std::string::npos;
  }

  inline std::string trim2( std::string s )
  {
    trim(s);
    return s;
  }

  inline std::string lowerCase( std::string s )
  {
    static_assert('A'+32 == 'a',"");
    for (auto& c : s)
      if ( c >= 'A' && c <= 'Z' )
        c += 32;
    return s;
  }

  inline std::string upperCase( std::string s )
  {
    static_assert('A'+32 == 'a',"");
    for (auto& c : s)
      if ( c >= 'a' && c <= 'z' )
        c -= 32;
    return s;
  }

  inline constexpr bool isAlpha( const char c ) noexcept
  {
    return ( c>='a' && c<='z' ) || ( c>='A' && c<='Z' );
  }

  inline constexpr bool isAlphaNumeric( const char c ) noexcept
  {
    return ( c>='a' && c<='z' ) || ( c>='A' && c<='Z' ) || ( c>='0' && c<='9' );
  }

  inline bool isAlphaNumeric( const std::string& s ) noexcept
  {
    for ( auto c : s )
      if (!isAlphaNumeric(c))
        return false;
    return true;
  }

  inline bool startswith(const std::string& str, char c ) { return !str.empty() && str.front()==c; }
  inline bool endswith(const std::string& str, char c ) { return !str.empty() && str.back()==c; }

  inline std::string str_rjust( std::string str, std::size_t n )
  {
    const auto ssize = str.size();
    if ( ssize >= n )
      return str;
    std::string res( n - ssize, ' ' );
    res += str;
    return res;
  }

  inline std::string str_ljust( std::string str, std::size_t n )
  {
    const auto ssize = str.size();
    if ( ssize <= n )
      str.append( n - ssize, ' ' );
    return str;
  }

  inline VectS trimEntries( VectS&& v )
  {
    VectS vv(std::move(v));
    for ( auto& e : vv )
      trim(e);
    return vv;
  }

  inline void streamJSON( std::ostream& os, json_null_t ) { os << "null"; }
  inline void streamJSON( std::ostream& os, bool b ) { os << ( b ? "true" : "false" ); }
  inline void streamJSON( std::ostream& os, const std::string&  s) { streamJSON(os,s.c_str()); }
  inline void streamJSON( std::ostream& os, const char * cstr ) { streamJSON(os,StrView(cstr)); }
  inline void streamJSON( std::ostream& os, float v ) { streamJSON(os,double(v)); }
  template<class T, typename std::enable_if<std::is_integral<T>::value>::type*>
  inline void streamJSON( std::ostream& os, T val )
  {
    os << val;
  }
  template<class T1, class T2>
  inline void streamJSON( std::ostream& os, const std::pair<T1,T2>& val )
  {
    os << '[';
    streamJSON(os,val.first);
    os << ',';
    streamJSON(os,val.second);
    os << ']';
  }

  template<class TContainer, typename T>
  inline void streamJSON( std::ostream& os, const TContainer& arr )
  {
    os << '[';
    bool first(true);
    for ( const auto& e : arr ) {
      if ( first ) {
        first = false;
      } else {
        os << ',';
      }
      streamJSON(os,e);
    }
    os << ']';
  }
  template<class T>
  inline void streamJSONDictEntry( std::ostream& os,
                                   const char * key ,
                                   const T& value,
                                   JSONDictPos pos)
  {
    os << (pos == JSONDictPos::FIRST?'{':',');
    streamJSON(os,key);
    os <<':';
    streamJSON(os,value);
    if (pos == JSONDictPos::LAST)
      os <<'}';
  }

  inline std::string joinstr(const VectS& parts )
  {
    return joinstr(parts,StrView::make(" "));
  }

  inline std::string joinstr(const Span<const StrView>& parts )
  {
    return joinstr(parts,StrView::make(" "));
  }

  inline std::string joinstr(const VectS& parts, StrView sep )
  {
    SmallVector<StrView,8> v;
    for ( auto& e : parts )
      v.emplace_back( e );
    return joinstr( Span<const StrView>(v), sep );
  }

}

#endif

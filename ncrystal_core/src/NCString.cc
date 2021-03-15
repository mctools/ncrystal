////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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

#include "NCrystal/internal/NCString.hh"
#include "NCrystal/internal/NCMath.hh"
#include <cstring>
#include <istream>
#include <iomanip>
#include <cstdlib>
namespace NC = NCrystal;

bool NC::isSimpleASCII(const std::string& input, bool allow_tab, bool allow_newline)
{
  nc_assert('\t'<32&&'\n'<32&&'\r'<32);
  const char * c = &input[0];
  const char * cE = c + input.size();
  for (;c!=cE;++c) {
    if (*c<32) {
      if (allow_tab&&*c=='\t')
        continue;
      if (allow_newline&&(*c=='\n'||*c=='\r'))
        continue;
      return false;
    }
    if ( *c>126 )
      return false;
  }
  return true;
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

bool NC::contains(const std::string& haystack, char needle )
{
  return haystack.find(needle) != std::string::npos;
}

bool NC::contains(const std::string& haystack, const std::string& needle)
{
  return haystack.find(needle) != std::string::npos;
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

namespace NCrystal {
  namespace {

  }
}

double NC::str2dbl(const std::string& s,const char * errmsg)
{
  double result;
  if ( !safe_str2dbl(s, result ) )
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not a double")<<": \""<<s<<"\"");
  return result;
}

int NC::str2int(const std::string& s,const char * errmsg)
{
  int result;
  if ( !safe_str2int(s, result ) )
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not an integer")<<": \""<<s<<"\"");
  return result;
}

bool NC::safe_str2dbl(const std::string& s, double& result )
{
  bool ok(true);
  double val;
  std::stringstream ss(s);
  ss >> val;
  while (!ss.fail()&&!ss.eof()) {
    char c;
    ss >> c;
    if (c!=' '&&c!='\t'&&c!='\n') {
      ok=false;
      break;
    }
  }
  if (!ok||ss.fail()||!ss.eof()) {
    if (!s.empty()) {
      //Always support "inf" and "INF" for cross-platform consistency (possibly prefixed with spaces).
      char last(s.at(s.size()-1));
      if (last=='f'||last=='F') {
        std::string tmp(s);
        trim(tmp);
        if (tmp=="inf"||tmp=="INF") {
          result = kInfinity;
          return true;
        }
        //TODO: Should we add support for "nan"/"infinity"?
      }
    }
    return false;
  }
  result = val;
  return true;
}

bool NC::safe_str2int(const std::string& s, int& result )
{
  int val;
  std::stringstream ss(s);
  ss >> val;
  while (!ss.fail()&&!ss.eof()) {
    char c;
    ss >> c;
    if (c!=' '&&c!='\t'&&c!='\n')
      return false;
  }
  if (ss.fail()||!ss.eof())
    return false;
  result = val;
  return true;
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

std::string NC::joinstr(const NC::VectS& parts, std::string separator)
{
  const std::size_t n = parts.size();
  if ( n < 2 )
    return n ? parts.at(0) : std::string();

  std::size_t newsize( (n-1)*separator.size() );
  for (const auto& p : parts)
    newsize += p.size();
  std::string tmp;
  tmp.reserve(newsize);

  tmp += parts.at(0);
  for (std::size_t i = 1; i < n; ++i) {
    tmp += separator;
    tmp += parts[i];
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

void NC::prettyPrintValue(std::ostream& os, double value, unsigned prec )
{
  auto ab = detectSimpleRationalNumbers(value);
  if ( ab.second == 1 ) {
    //print integer
    os << ab.first;
  } else if ( ab.second == 0 )  {
    //print value
    if ( prec == 0 ) {
      os << value;
    } else {
      //Must set precision. Make sure we don't change this globally:
      std::ostringstream ss;
      ss << std::setprecision(prec) << value;
      os << ss.str();
    }
  } else {
    //Fraction detected!
    os << ab.first<<"/"<<ab.second;
  }
}

std::string NC::prettyPrintValue2Str(double value, unsigned prec )
{
  std::ostringstream ss;
  prettyPrintValue(ss,value,prec);
  return ss.str();
}

//Common access to environment variables (unset and empty vars both return
//empty strings / 0.0 / 0, depending on which function is used):
std::string NC::ncgetenv(std::string v, std::string defval)
{
  nc_assert(!v.empty());
  const char * ev = std::getenv(("NCRYSTAL_"_s+v).c_str());
  return ev ? std::string(ev) : defval;
}

double NC::ncgetenv_dbl(std::string v, double defval )
{
  nc_assert(!v.empty());
  std::string varname("NCRYSTAL_");
  varname += v;
  const char * ev = std::getenv(varname.c_str());
  if (!ev)
    return defval;
  double result;
  if ( !safe_str2dbl(ev, result ) )
    NCRYSTAL_THROW2(BadInput,"Invalid value of environment variable "<<varname
                    <<" (expected a floating point number but got \""<<ev<<"\").");
  return result;
}

int NC::ncgetenv_int(std::string v, int defval )
{
  nc_assert(!v.empty());
  std::string varname("NCRYSTAL_");
  varname += v;
  const char * ev = std::getenv(varname.c_str());
  if (!ev)
    return defval;
  int result;
  if ( !safe_str2int(ev, result ) )
    NCRYSTAL_THROW2(BadInput,"Invalid value of environment variable "<<varname
                    <<" (expected an integral number but got \""<<ev<<"\").");
  return result;
}

bool NC::ncgetenv_bool(std::string v)
{
  nc_assert(!v.empty());
  nc_assert(!startswith(v,"NCRYSTAL_"));
  std::string varname("NCRYSTAL_");
  varname += v;
  const char * ev = std::getenv(varname.c_str());
  if (!ev)
    return false;
  std::string evs(ev);
  if (evs.size()==1) {
    if (evs[0]=='0')
      return false;
    if (evs[0]=='1')
      return true;
  }
  NCRYSTAL_THROW2(BadInput,"Invalid value of environment variable "<<varname
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

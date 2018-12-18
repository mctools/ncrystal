////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2018 NCrystal developers                                   //
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

#include "NCString.hh"
#include "NCrystal/NCDefs.hh"
#include <cstring>
#include <limits>
#include <istream>

bool NCrystal::isSimpleASCII(const std::string& input, bool allow_tab, bool allow_newline)
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

void NCrystal::trim( std::string& input )
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

void NCrystal::split(std::vector<std::string>& output, const std::string& input, std::size_t maxsplit, char sep )
{
  bool keep_empty = (sep!=0);
  //Split input on whitespace (" \t\n") with minimal memory allocations and
  //copying, discarding whitespace which is either leading and trailing or
  //repeating.
  nc_assert(bool(sep!='\0')==bool(sep));
  output.clear();

  if (input.empty()) {
    if (sep)
#if __cplusplus >= 201103L
      output.emplace_back();
#else
      output.push_back(std::string());
#endif
    return;
  }

  if (sep&&input[0]==sep)
#if __cplusplus >= 201103L
    output.emplace_back();
#else
    output.push_back(std::string());
#endif
  const char * c = &input[0];
  const char * cE = c + input.size();
  const char * partbegin = 0;
  while (true) {
    if ( maxsplit && output.size() == maxsplit ) {
#if __cplusplus >= 201103L
      output.emplace_back(c);
#else
      output.push_back(c);
#endif
      return;
    }
    if (c==cE||(sep?*c==sep:(*c==' '||*c=='\t'||*c=='\n'||*c=='\r'))) {
      if (partbegin) {
        if (keep_empty||c>partbegin) {
#if __cplusplus >= 201103L
          output.emplace_back(partbegin,c-partbegin);
#else
          output.push_back(std::string(partbegin,c-partbegin));
#endif
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

std::string NCrystal::basename(const std::string& filename)
{
  std::size_t p = filename.rfind('/');
  return p+1>filename.size() ? filename : filename.substr(p+1);
}

std::string NCrystal::getfileext(const std::string& filename)
{
  std::string bn = basename(filename);
  std::size_t p = bn.rfind('.');
  return p == std::string::npos ? std::string() : bn.substr(p+1);
}

bool NCrystal::startswith(const std::string& str, const std::string& substr)
{
  return str.size()>=substr.size() && strncmp(str.c_str(), substr.c_str(), substr.size() ) == 0;
}

bool NCrystal::contains(const std::string& haystack, char needle )
{
  return haystack.find(needle) != std::string::npos;
}

bool NCrystal::contains(const std::string& haystack, const std::string& needle)
{
  return haystack.find(needle) != std::string::npos;
}

bool NCrystal::contains_any(const std::string& haystack, const std::string& needles)
{
  const char * cB = &needles[0];
  const char * cE = cB + needles.size();
  for (const char * c = cB; c!=cE; ++c)
    if (haystack.find(*c) != std::string::npos)
      return true;
  return false;
}

double NCrystal::str2dbl(const std::string& s,const char * errmsg)
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
        if (tmp=="inf"||tmp=="INF")
          return std::numeric_limits<double>::infinity();
        //TODO for NC2: Should we add support for "nan"/"infinity"?
      }
    }
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not a double")<<": \""<<s<<"\"");
  }
  return val;
}

int NCrystal::str2int(const std::string& s,const char * errmsg)
{
  bool ok(true);
  int val;
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
  if (!ok||ss.fail()||!ss.eof())
    NCRYSTAL_THROW2(BadInput,(errmsg?errmsg:"Invalid number in string is not an integer")<<": \""<<s<<"\"");
  return val;
}

std::ifstream& NCrystal::ignoreCharNTimes(std::ifstream& file, unsigned num, const char& c)
{
  //TODO for NC2: we probably don't really want this function, since it is very
  //accepting of errors in input rather than defensively sanitising (if
  //removing, also remove #include <istream>)
  for(unsigned i=0; i < num ; ++i){
    file.ignore(std::numeric_limits<std::streamsize>::max(),c);
  }
  return file;
}

bool NCrystal::contains_only(const std::string& haystack, const std::string& needles)
{
  const char * c = &haystack[0];
  const char * cE = c + haystack.size();
  for (;c!=cE;++c) {
    if (!contains(needles,*c))
      return false;
  }
  return true;
}

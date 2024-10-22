////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

#include "NCrystal/internal/NCFileUtils.hh"
#include "NCrystal/internal/NCString.hh"
#include "NCFileUtilsWin.hh"
#include <fstream>

namespace NC = NCrystal;

#if defined(__unix__) || (defined (__APPLE__) && defined (__MACH__))
#include <stdlib.h>
#include <limits.h>
std::string NC::tryRealPath( const std::string& fn ) {
#ifndef PATH_MAX
#  define PATH_MAX 1024
#endif
  if ( fn.size() > (PATH_MAX-1) )
    return {};
  char buf[PATH_MAX+1];
  char * res = realpath(fn.c_str(), buf);
  if (!res)
    return {};//failure
  return std::string(res);
}
#else
//Windows or whatever... just give up. [FIXME: To NCFileUtilsWin.cc + add NCFileUtilsUnix.cc]
std::string NC::tryRealPath( const std::string& ) { return {}; }
#endif


NC::Optional<std::string> NC::readEntireFileToString( const std::string& path )
{
  //FIXME: Give error on utf16/utf32 BOM, and skip past any utf8 BOM? See https://www.unicode.org/faq/utf_bom.html#BOM

  //Read entire file into a string while protecting against someone mistakenly
  //trying to open a multi-gigabyte file and bringing their machine to a slow
  //halt.
  using size_type = std::streamsize;
  constexpr size_type read_block_size = 4096;
  static const size_type maxread_megabytes = static_cast<size_type>(ncgetenv_int("MAXREAD_MEGABYTES",100));
  static const size_type maxread_bytes = maxread_megabytes*1048576ull;
  size_type maxread_blocks = maxread_bytes/read_block_size + 1;
  static std::ifstream::char_type buffer[read_block_size] = {};
  auto mode = std::ios_base::binary | std::ios_base::in;
#ifdef NCRYSTAL_USE_WINDOWS_FILEUTILS
  std::ifstream fh = WinFileUtils::open_ifstream_from_path( path, mode );
#else
  std::ifstream fh( path, mode );
#endif
  if ( !fh.good() )
    return NullOpt;//interpret as if file does not exist or is not readable.
  std::string out;
  out.reserve(read_block_size*8);
  while ( fh.read( &buffer[0], read_block_size ) ) {
    out.append( buffer, fh.gcount() );
    if ( maxread_blocks-- == 1 )
      NCRYSTAL_THROW2(DataLoadError,"NCrystal: File too large to read (max size allowed is "
                      <<maxread_megabytes
                      <<"MB - increase by setting NCRYSTAL_MAXREAD_MEGABYTES env. var): "<<path);
  }
  if ( fh.gcount() )
    out.append( buffer, fh.gcount() );
  out.shrink_to_fit();
  return Optional<std::string>(std::move(out));
}

bool NC::path_is_absolute( const std::string& p )
{
  if (p.empty())
    return false;
#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(_WIN32) || defined(__CYGWIN__)//fixme: to NCFileUtilsWin.cc
  //FIXME: We need to consider that c:bla.ncmat means "relative to CWD of C:" not an absolute path!! We should unit test this!!
  if ( p.at(0)=='/' )
    return true;//who are we to argue if people are somehow using forward-slashes on a windows platform...
  //Try to catch e.g. "C:/" or "C:\":
  return p.size()>3 && p.at(1)==':' && ( p.at(2)=='/' || p.at(2)=='\\' ) ;
#else
  return p.at(0)=='/';
#endif
}
std::string NC::path_join(const std::string& p1, const std::string& p2)
{
#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(_WIN32) || defined(__CYGWIN__)//fixme: to NCFileUtilsWin.cc
  if ( contains(p1,'/') || contains(p2,'/') )
    return p1+'/'+p2;//seems like user is already using forward slashes, so assume it is ok.
  return p1+'\\'+p2;
#else
  return p1+'/'+p2;
#endif
}

std::string NC::basename(const std::string& filename)
{
  std::size_t p = filename.rfind('/');
#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(_WIN32) || defined(__CYGWIN__)
  std::size_t p2 = filename.rfind('\\');
  if ( p2 != std::string::npos && ( p == std::string::npos || p2 > p ) )
    p = p2;
#endif
  return p+1>filename.size() ? filename : filename.substr(p+1);
}

std::string NC::getfileext(const std::string& filename)
{
  std::string bn = basename(filename);
  std::size_t p = bn.rfind('.');
  return p == std::string::npos ? std::string() : bn.substr(p+1);
}

bool NC::file_exists(const std::string& name)
{
#ifdef NCRYSTAL_USE_WINDOWS_FILEUTILS
  return WinFileUtils::file_exists( name );
#else
  std::ifstream f(name.c_str());
  return f.good();
#endif
}

#ifdef NCRYSTAL_USE_WINDOWS_FILEUTILS
NC::VectS NC::ncglob(const std::string& pattern)
{
  return WinFileUtils::ncglob(pattern);
}
#else
//POSIX globbing:
#include <glob.h>
NC::VectS NC::ncglob(const std::string& pattern) {
  VectS result;
  glob_t pglob;
  int retval = glob(pattern.c_str(),0,0, &pglob);
  if ( retval != 0 && retval != GLOB_NOMATCH )
    NCRYSTAL_THROW2(CalcError,"Error encountered while"
                    " attempting to glob for \""<<pattern<<"\"");
  if ( retval != GLOB_NOMATCH ) {
    for ( decltype(pglob.gl_pathc) i = 0; i < pglob.gl_pathc; ++i ) {
      auto pv = pglob.gl_pathv[i];
      if ( pv ) {
        std::string s(pv);
        if ( !s.empty() )
          result.push_back(s);
      }
    }
    std::sort(result.begin(),result.end());
  }
  globfree(&pglob);
  return result;
}
#endif

#ifdef NCRYSTAL_USE_WINDOWS_FILEUTILS
std::string NC::ncgetcwd()
{
  return WinFileUtils::get_current_working_dir();
}
#else
//POSIX getcwd:
#include <unistd.h>
std::string NC::ncgetcwd() {
  constexpr std::size_t nfixbuf = 4096;
  char buffix[nfixbuf];
  if (getcwd(&buffix[0], nfixbuf))
    return std::string(&buffix[0]);
  if (errno == ERANGE) {
    //crazy system with crazy long path.
    constexpr std::size_t nlarge = 131072;
#if nc_cplusplus >= 201402L
      //Our make_unique for c++11 seems to have problems with arrays
    auto largebuf = std::make_unique<char[]>(nlarge);
#else
    std::unique_ptr<char[]> largebuf(new char[nlarge]());
#endif
    if (getcwd(&largebuf[0], nlarge))
      return std::string(&largebuf[0]);
    if (errno == ERANGE)
      NCRYSTAL_THROW(CalcError,"current working directory is too long");
  }
  NCRYSTAL_THROW(CalcError,"Could not determine current working directory");
}
#endif

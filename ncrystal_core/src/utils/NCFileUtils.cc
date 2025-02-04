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

#include "NCrystal/internal/utils/NCFileUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCFileUtilsWin.hh"
#include <fstream>

namespace NC = NCrystal;

#ifdef NCRYSTAL_USE_WINDOWS_FILEUTILS
std::string NC::tryRealPath( const std::string& path )
{
  return WinFileUtils::get_absolute_path(path);
}
#elif defined(__unix__) || (defined (__APPLE__) && defined (__MACH__))
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
//Give up:
std::string NC::tryRealPath( const std::string& path )
{
  return {};
}
#endif


NC::Optional<std::string> NC::readEntireFileToString( const std::string& path )
{
  //TODO: Give error on utf16/utf32 BOM, and skip past any utf8 BOM? See
  //https://www.unicode.org/faq/utf_bom.html#BOM

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

namespace NCRYSTAL_NAMESPACE {
  namespace {
    bool is_path_sep( char c ) {
      return c == '/' || c == '\\';
    }
    StrView trim_pathseps_from_edges( StrView p )
    {
      StrView::size_type i = 0;
      StrView::size_type n = p.size();
      while ( n && is_path_sep(p[n-1]) )
        --n;
      while ( i < n && is_path_sep(p[i]) )
        ++i;
      return p.substr( i, n );
    }
    struct AnalysedPath {
      bool m_is_absolute;
      VectS m_parents;
      std::string m_filename;
      std::string m_windows_drive_name;
      void append_dirname( std::string& res ) const
      {
        if ( !m_windows_drive_name.empty() ) {
          res += m_windows_drive_name;
          res += ':';
        }
        char sep = pathsep();
        if ( m_is_absolute )
          res += sep;
        for ( auto& i : ncrange(m_parents.size()) ) {
          if ( i )
            res += sep;
          res += m_parents[i];
        }
      }

    public:
      char pathsep() const
      {
#ifdef NCRYSTAL_USE_WINDOWS_FILEUTILS
        return '\\';
#else
        return ( m_windows_drive_name.empty() ? '/' : '\\' );
#endif
      }
      bool is_absolute() const { return m_is_absolute; }
      const std::string& filename() const { return m_filename; }
      std::string dirname() const
      {
        std::string res;
        res.reserve(128);
        append_dirname(res);
        res.shrink_to_fit();
        return res;
      }

      //Optional<AnalysedPath> parent() const;
      std::string encode() const
      {
        std::string res;
        res.reserve(128);
        append_dirname(res);
        if ( m_filename.empty() ) {
          //Return a string indicating the parent directory, avoiding an empty
          //string:
          if ( res.empty() )
            res += ( m_is_absolute ? pathsep() : '.' );
        } else {
          if ( !m_parents.empty() )
            res += pathsep();
          res += m_filename;
        }
        res.shrink_to_fit();
        return res;
      }
      AnalysedPath( const std::string& filename_str )
      {
        auto filename = StrView(filename_str);

        //Determine (and peel off) windows drive letter:
        if ( filename.size() >= 2
             && filename[1] == ':'
             && isAlphaNumeric(filename[0])
             && filename[0] >= 'A' )
        {
          m_windows_drive_name = filename[0];
          filename = filename.substr(2);
        }

        //Determine if absolute, i.e. starts with a pathsep:
        m_is_absolute = !filename.empty() && is_path_sep(filename[0]);
        //^^^^^ TODO: this is not always correct for windows shared drive paths??

        //Discard path seps at both ends:
        filename = trim_pathseps_from_edges( filename );

        //Now, split by pathsep:
        auto parts = filename.split_any<8,
                                        StrView::SplitKeepEmpty::No,
                                        StrView::SplitTrimParts::No>("/\\");
        if (!parts.empty()) {
          m_filename = parts.back().to_string();
          parts.pop_back();
        }
        m_parents.reserve(parts.size());
        for ( auto& p : parts ) {
          if ( p==".." && !m_parents.empty() )
            m_parents.pop_back(); // ".." steps up from the previous
          else
            m_parents.emplace_back( p.to_string() );
        }
        //Finally, as a special case, we discard "." from filenames:
        if ( m_filename == "." ) {
          m_filename.clear();
          if ( m_parents.empty() && !m_is_absolute )
            m_parents.push_back(".");
        }
        m_parents.shrink_to_fit();
      }
    };

  }
}

bool NC::path_is_absolute( const std::string& p )
{
  return AnalysedPath( p ).is_absolute();
}

std::string NC::basename( const std::string& p )
{
  return AnalysedPath( p ).filename();
}

std::string NC::dirname( const std::string& p )
{
  return AnalysedPath( p ).dirname();
}

std::string NC::normalise( const std::string& p )
{
  return AnalysedPath( p ).encode();
}

std::string NC::getfileext( const std::string& p )
{
  std::string bn = basename( p );
  std::size_t idx = bn.rfind('.');
  return idx == std::string::npos ? std::string() : bn.substr(idx+1);
}

std::string NC::path_join(const std::string& p1, const std::string& p2)
{
#ifdef NCRYSTAL_USE_WINDOWS_FILEUTILS
  if ( contains(p1,'/') || contains(p2,'/') )
    return p1+'/'+p2;//seems like user is already using forward slashes, so
                     //assume it is ok.
  return p1+'\\'+p2;
#else
  return p1+'/'+p2;
#endif
}

bool NC::file_exists( const std::string& p )
{
#ifdef NCRYSTAL_USE_WINDOWS_FILEUTILS
  return WinFileUtils::file_exists( p );
#else
  std::ifstream f(p.c_str());
  return f.good();
#endif
}

#ifndef NCRYSTAL_USE_WINDOWS_FILEUTILS
//POSIX globbing:
#  include <glob.h>
namespace NCRYSTAL_NAMESPACE {
  namespace {
    VectS ncglob_posix_impl(const std::string& pattern) {
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
            if ( !s.empty() && s != "." && s!=".." )
              result.push_back(s);
          }
        }
      }
      globfree(&pglob);
      return result;
    }
  }
}
#endif

NC::VectS NC::ncglob( const std::string& pattern )
{
  if ( pattern.empty() )
    NCRYSTAL_THROW(BadInput,"ncglob does not work with empty patterns");
  AnalysedPath path( pattern );
  if ( path.filename().empty() )
    NCRYSTAL_THROW2(BadInput,"ncglob could not decode pattern: \""
                    <<pattern<<"\"");
  auto contains_wildcard = [](const std::string& s )
  {
    return contains(s,'*') || contains(s,'?');
  };
  std::string path_dirname = path.dirname();
  if ( contains_wildcard( path_dirname ) )
    NCRYSTAL_THROW(BadInput,"ncglob only supports wildcards in the last"
                   " file or directory name");
  if ( !contains_wildcard( path.filename() ) ) {
    //Special case, no wildcards:
    std::string fn = path.encode();
    VectS nonglob_res;
    if ( file_exists(fn) )
      nonglob_res.push_back(std::move(fn));
    return nonglob_res;
  }
  //Ok, some actual work:
#ifdef NCRYSTAL_USE_WINDOWS_FILEUTILS
  VectS res = WinFileUtils::ncglob_impl(pattern);
  //Windows globbing only gives the filenames, we have to re-inject the
  //dirnames:
  if ( !path_dirname.empty() ) {
    for ( auto& e : res )
      e = path_join( path_dirname, e );
  }
#else
  VectS res = ncglob_posix_impl(pattern);
#endif
  res.shrink_to_fit();
  std::sort(res.begin(),res.end());
  return res;
}

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

#if defined(__unix__) && !defined(__APPLE__)
namespace NCRYSTAL_NAMESPACE {
  namespace {
    std::string try_get_self_exe_path_from_proc(const char * procpath)
    {
      //Read a link like /proc/self/exe
      char buf[65536+1];//PATH_MAX is unreliable so we use huge buffer
      auto len = ::readlink(procpath, buf, sizeof(buf)-1 );
      if ( len > 0 && std::size_t(len+1) < std::size_t(sizeof(buf)) ) {
        buf[len] = '\0';//readlink does not add terminating null char
        return { buf };
      } else {
        return {};//error
      }
    }
  }
}
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h>
namespace NCRYSTAL_NAMESPACE {
  namespace {
    std::string get_self_exe_path_apple()
    {
      //First find required buffer size:
      uint32_t bufsize = 0;
      char fakebuf[4];
      auto status = _NSGetExecutablePath(fakebuf, &bufsize);
      if ( status != -1 )
        return {};//unexpected result, buffer should NOT be large enough
      std::string res;
      res.resize( bufsize, 0 );
      status = _NSGetExecutablePath(&res[0], &bufsize);
      if ( status == 0 && res.find('\0') > res.size() )
        return res;
      else
        return {};//error
    }
  }
}
#endif

std::string NC::determine_exe_self_path( int argc, char** argv )
{
  //First try some platform specific methods not relying on argv0:
#ifdef NCRYSTAL_USE_WINDOWS_FILEUTILS
  {
    auto res = WinFileUtils::get_self_exe_path_windows();
    if (!res.empty())
      return res;
  }
#endif
#ifdef __APPLE__
  {
    auto res = get_self_exe_path_apple();
    if (!res.empty())
      return res;
  }
#endif
#if defined(__unix__) && !defined(__APPLE__)
  {
    //Should always work on linux:
    auto res = try_get_self_exe_path_from_proc("/proc/self/exe");
    if (!res.empty())
      return res;
  }
  {
    //Might occasionally work on FreeBSD etc.
    auto res = try_get_self_exe_path_from_proc("/proc/curproc/file");
    if (!res.empty())
      return res;
  }
#endif
  //Fall back to argv0:
  if ( !(argc > 0) )
    return {};//not available (rare, but allowed by C standard);
  std::string argv0( argv[0] );
  if ( NC::path_is_absolute( argv0 ) )
    return argv0;
  std::string guess = NC::path_join( ncgetcwd(), argv0 );
  if ( NC::file_exists(guess) )
    return guess;
  //Give up:
  return {};
}

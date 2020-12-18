////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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

#include "NCrystal/NCFile.hh"

#include <fstream>
#include <cstdlib>

namespace NC = NCrystal;

namespace NCrystal {
  namespace {
    static std::unique_ptr<TextInputManager> s_textInputMgr;
    static std::mutex s_textInputMgr_mutex;//call pfacbFireCallbacks() just before locking this

    std::string path_join(const std::string p1, const std::string p2)
    {
#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(_WIN32) || defined(__CYGWIN__)
      return p1+'\\'+p2;
#else
      return p1+'/'+p2;
#endif
    }

    VectS initDataPathList()
    {
      //Path search order
      VectS list;

      //First take the NCRYSTAL_DATADIR. Prefer the environment variable to the
      //compile-time defined path.

      const char * envpath = std::getenv("NCRYSTAL_DATADIR");
      std::string ncrystal_datadir = envpath ? envpath : "";
#ifdef NCRYSTAL_DATADIR
#  define NCRYSTAL_str(s) #s
#  define NCRYSTAL_xstr(s) NCRYSTAL_str(s)
      //If NCRYSTAL_DATADIR defined *and* environment variable NCRYSTAL_DATADIR is
      //NOT set, check the defined path (thus, setting the environment variable
      //means we no longer look in the compile-time provided NCRYSTAL_DATADIR).
      if (ncrystal_datadir.empty())
        ncrystal_datadir = NCRYSTAL_xstr(NCRYSTAL_DATADIR);
#endif
      if (!ncrystal_datadir.empty())
        list.push_back(ncrystal_datadir);

      //Now the data path:
      const char * envpathlist = std::getenv("NCRYSTAL_DATA_PATH");
      std::string pathlist = envpathlist ? envpathlist : "";

      while (!pathlist.empty()) {
        //A bit annoying that we can't simply use split from NCString.hh here...
        auto i = pathlist.find(':');
        std::string pathtocheck;
        if (i==0) {
          //leading ':', skip it:
          pathlist = pathlist.substr(1);
          continue;
        }
        if (i == std::string::npos) {
          //No ':', this is the last entry:
          std::swap(pathlist,pathtocheck);
        } else {
          pathtocheck = pathlist.substr(0,i);
          pathlist = pathlist.substr(i);
        }
        if (!pathtocheck.empty())
          list.push_back(pathtocheck);
      }


      //Remove duplicates, if any:
      std::set<std::string> seen;
      VectS finalresult;
      finalresult.reserve(list.size());
      for (auto e : list) {
        if (!seen.count(e)) {
          finalresult.emplace_back(e);
          seen.insert(e);
        }
      }
      finalresult.shrink_to_fit();
      return finalresult;

    }
    const VectS& getDataPathList() {
      static VectS thePath = initDataPathList();
      return thePath;
    }


    struct PreFileAccessCallbackCache {
      std::mutex mtx;
      std::vector<std::function<void()>> callbacks;
      std::atomic<bool> duringFire{false};
    };
    PreFileAccessCallbackCache& pfacbCache() {
      static PreFileAccessCallbackCache s_cache;
      return s_cache;
    }
    void pfacbFireCallbacks() {
      auto& cache =pfacbCache();
      if (cache.duringFire)
        return;//already firing
      std::lock_guard<std::mutex> guard(cache.mtx);
      cache.duringFire = true;
      while (!cache.callbacks.empty()) {
        cache.callbacks.back()();
        cache.callbacks.pop_back();
      }
      cache.duringFire = false;
    }
  }
}

void NC::addPreFileAccessCallback(std::function<void()> callback)
{
  auto& cache =pfacbCache();
  std::lock_guard<std::mutex> guard(cache.mtx);
  cache.callbacks.push_back(callback);
}

std::string NC::basename(const std::string& filename)
{
  std::size_t p = filename.rfind('/');
  return p+1>filename.size() ? filename : filename.substr(p+1);
}

std::string NC::getfileext(const std::string& filename)
{
  std::string bn = basename(filename);
  std::size_t p = bn.rfind('.');
  return p == std::string::npos ? std::string() : bn.substr(p+1);
}

bool NC::file_exists(const std::string& name) {
  //Only portable way in C++98 is to attempt to open the file.
  std::ifstream f(name.c_str());
  return f.good();
}

std::string NC::find_file(const std::string& filename) {

  if (filename.empty())
    return std::string();

  if (file_exists(filename))
    return filename;

  if (filename.at(0)=='/')
    return std::string();//don't look further if absolute path (this test might fail on windows...).

  for ( const auto& search_path : getDataPathList() )  {
    std::string tmp = path_join(search_path,filename);
    if (file_exists(tmp))
      return tmp;
  }

  //not found.
  return std::string();
}

std::vector<NC::FileListEntry> NC::listAvailableFiles(std::string extension)
{
  std::vector<FileListEntry> result;
  result.reserve(128);
  std::set<std::string> seen;
  std::set<std::string> seend;
  auto addentry = [&seen,&result](const std::string& name, const std::string& source)
  {
    result.push_back(FileListEntry{name,source,!seen.insert(name).second});
  };

  {
    pfacbFireCallbacks();
    std::lock_guard<std::mutex> guard(s_textInputMgr_mutex);
    if (!!s_textInputMgr) {
      for ( auto& e : s_textInputMgr->getList() ) {
        addentry(e.first,e.second);
      }
    }
  }

  std::string globpattern = "*."_s + extension;
  auto addfromdir = [&addentry,&seend,&globpattern](const std::string& d,const std::string& name )
  {
    if ( d.empty() || !seend.insert(d).second )
      return;
    for ( auto& fn : ncglob(path_join(d,globpattern)))
      addentry(basename(fn),name);
    //Special case, look in ncplugin/ subdir if it is there (to support plugin development):
    std::string subdir("ncplugin");
    auto d2 = path_join(d,subdir);
    if (file_exists(d2))
      for ( auto& fn : ncglob(path_join(d2,globpattern)))
        addentry(path_join(subdir,basename(fn)),name);
  };

  addfromdir(ncgetcwd(),"<current-working-directory>");
  for ( auto& dp : getDataPathList() )
    addfromdir(dp,dp);

  return result;
}



namespace NCrystal {
  TextInputStream::TextInputStream(const std::string& description)
    : m_descr(description)
  {
  }

  const std::string& TextInputStream::onDiskResolvedPath() const
  {
    static std::string str_empty;
    return str_empty;
  }


  class FileTextInputStream : public TextInputStream {
  public:

    virtual ~FileTextInputStream(){}

    FileTextInputStream(const std::string& filename)
      : TextInputStream(filename),
        m_file(filename)
    {
      if (!m_file.good())
        NCRYSTAL_THROW2(FileNotFound,"Failure while trying to open file "<<filename);
      tryReadNext();
    }

    const std::string& onDiskResolvedPath() const final
    {
      return description();
    }

    virtual bool getLine(std::string& line) {
      if (!moreLines()) {
        line.clear();
        return false;
      }
      line = std::move(m_nextLine);
      tryReadNext();
      return true;
    }

    virtual bool moreLines() const
    {
      return !m_nextLine.empty() || m_file.is_open();
    }

    virtual const char * streamType() const
    {
      return "file";
    }

  private:
    void tryReadNext() {
      if (m_file.is_open()) {
        if (!std::getline(m_file,m_nextLine)) {
          m_file.close();
          m_nextLine.clear();
        }
      } else {
        m_nextLine.clear();
      }
    }
    std::ifstream m_file;
    std::string m_nextLine;
  };

  class MemBufTextInputStream : public TextInputStream {
  public:

    virtual ~MemBufTextInputStream(){}

    MemBufTextInputStream(const std::string& buffername,
                          const std::string& buffer)
      : TextInputStream(buffername),
        m_stream(buffer),
        m_more(true)
    {
      tryReadNext();
    }

    virtual bool getLine(std::string& line) {
      if (!m_more) {
        line.clear();
        return false;
      }
      line = std::move(m_nextLine);
      tryReadNext();
      return true;
    }

    virtual bool moreLines() const
    {
      return m_more;
    }

    virtual const char * streamType() const
    {
      return "memory-buffer";
    }

  private:
    void tryReadNext() {
      if (m_more && std::getline(m_stream,m_nextLine))
        return;
      m_more = false;
      m_nextLine.clear();
    }
    std::stringstream m_stream;
    std::string m_nextLine;
    bool m_more;
  };

}

std::unique_ptr<NC::TextInputStream> NC::createTextInputStream( const std::string& sourcename )
{

  {
    pfacbFireCallbacks();
    std::lock_guard<std::mutex> guard(s_textInputMgr_mutex);
    if (!!s_textInputMgr) {
      auto stream = s_textInputMgr->createTextInputStream( sourcename );
      if (stream!=nullptr)
        return stream;
      if (!s_textInputMgr->allowFallbackToUsualDefaults())
        NCRYSTAL_THROW2(FileNotFound,"Could not find input corresponding to name: "<<sourcename);
    }
  }

  //Fall back to looking for on-disk sources, where sourcename are the filenames
  //passed to find_file:
  std::string resolved_name = find_file(sourcename);
  if (resolved_name.empty())
    NCRYSTAL_THROW2(FileNotFound,"Could not find input file: "<<sourcename);
  return createTextInputStreamFromFile(resolved_name);
}

std::unique_ptr<NC::TextInputStream> NC::createTextInputStreamFromBuffer( const std::string& buffername,
                                                                      const std::string& buffer )
{
  return std::make_unique<MemBufTextInputStream>(buffername,buffer);
}

std::unique_ptr<NC::TextInputStream> NC::createTextInputStreamFromFile( const std::string& filepath )
{
  return std::make_unique<FileTextInputStream>(filepath);
}

NC::TextInputManager::TextInputManager()
{
}

NC::TextInputManager::~TextInputManager()
{
}

void NC::registerTextInputManager( std::unique_ptr<TextInputManager> mgr )
{
  pfacbFireCallbacks();
  std::lock_guard<std::mutex> guard(s_textInputMgr_mutex);
  s_textInputMgr = std::move(mgr);
}

#if ( defined (_WIN32) || defined (WIN32) ) && !defined (__CYGWIN__)
//Windows globbing -> untested apart from compilation on godbolt.org
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
NC::VectS NC::ncglob(const std::string& pattern) {
  VectS result;
  WIN32_FIND_DATA fdata;
  HANDLE fh = FindFirstFileA(pattern.c_str(), &fdata);
  if (fh == INVALID_HANDLE_VALUE)
    return result;
  while (true) {
    result.push_back(fdata.cFileName);
    if (!FindNextFileA(fh, &fdata))
      break;
  }
  FindClose(fh);
  std::sort(result.begin(),result.end());
  return result;
}
//Windows getcwd:
std::string NC::ncgetcwd() {
    char buff[MAX_PATH];
    GetModuleFileName( NULL, buff, MAX_PATH );
    string::size_type position = string( buff ).find_last_of( "\\/" );
    return string( buff ).substr( 0, position);
}
#else
//POSIX globbing:
#include <glob.h>
NC::VectS NC::ncglob(const std::string& pattern) {
  VectS result;
  glob_t pglob;
  int retval = glob(pattern.c_str(),0,0, &pglob);
  if ( retval != 0 && retval != GLOB_NOMATCH )
    NCRYSTAL_THROW2(CalcError,"Error encountered while attempting to glob for \""<<pattern<<"\"");
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
#if __cplusplus >= 201402L
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

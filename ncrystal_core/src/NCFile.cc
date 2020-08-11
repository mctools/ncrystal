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

  static UniquePtr<TextInputManager> s_textInputMgr;

  std::string path_join(const std::string p1, const std::string p2)
  {
#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(_WIN32) || defined(__CYGWIN__)
    return p1+'\\'+p2;
#else
    return p1+'/'+p2;
#endif
  }
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

  const char * envpath = std::getenv("NCRYSTAL_DATADIR");
  std::string path = envpath ? envpath : "";
  if (!path.empty()) {
    std::string tmp = path_join(path,filename);
    if (file_exists(tmp))
      return tmp;
  }
#ifdef NCRYSTAL_DATADIR
#  define NCRYSTAL_str(s) #s
#  define NCRYSTAL_xstr(s) NCRYSTAL_str(s)
  path = NCRYSTAL_xstr(NCRYSTAL_DATADIR);
  if (!path.empty()) {
    std::string tmp = path_join(path,filename);
    if (file_exists(tmp))
      return tmp;
  }
#endif


  //not found.
  return std::string();
}

namespace NCrystal {
  TextInputStream::TextInputStream(const std::string& description)
    : m_descr(description)
  {
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
    std::fstream m_file;
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

void NC::createTextInputStream( const std::string& sourcename,
                                NC::UniquePtr<NC::TextInputStream>& streamptr )
{
  streamptr.set(0);

  TextInputManager * mgr = s_textInputMgr.obj();
  if (mgr) {
    mgr->createTextInputStream( sourcename, streamptr );
    if (streamptr.obj())
      return;
    if (!mgr->allowFallbackToUsualDefaults())
      NCRYSTAL_THROW2(FileNotFound,"Could not find input corresponding to name: "<<sourcename);
  }

  //For now only support on-disk sources, where sourcename are the filenames:
  std::string resolved_name = find_file(sourcename);
  if (resolved_name.empty())
    NCRYSTAL_THROW2(FileNotFound,"Could not find input file: "<<sourcename);
  createTextInputStreamFromFile(resolved_name,streamptr);
}

void NC::createTextInputStreamFromBuffer( const std::string& buffername,
                                          const std::string& buffer,
                                          NC::UniquePtr<NC::TextInputStream>& streamptr )
{
  streamptr.set(0);
  streamptr.set(static_cast<TextInputStream*>(new MemBufTextInputStream(buffername,buffer)));
}

void NC::createTextInputStreamFromFile( const std::string& filepath,
                                        NC::UniquePtr<NC::TextInputStream>& streamptr )
{
  streamptr.set(0);
  streamptr.set(static_cast<TextInputStream*>(new FileTextInputStream(filepath)));
}

NC::TextInputManager::TextInputManager()
{
}

NC::TextInputManager::~TextInputManager()
{
}

void NC::registerTextInputManager( TextInputManager* mgr )
{
  s_textInputMgr.set(mgr);
}

#ifndef NCrystal_File_hh
#define NCrystal_File_hh

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

#include "NCrystal/NCMem.hh"

namespace NCrystal {

  //Check if file exists and is readable:
  bool file_exists( const std::string& filename );

  //Search for the file. If filename does not exist relative to the current
  //working directory, is not an absolute path to the file, first look
  //relatively to the directory NCRYSTAL_DATADIR and secondly relatively to a
  //directory (if any) given at compilation time by the preprocessor define
  //-DNCRYSTAL_DATADIR=/some/dir. Returns empty string if not found:
  std::string find_file( const std::string& filename );

  //Interface which abstracts text sources, allowing common interface for
  //reading data from on-disk files and from e.g. in-memory databases. The
  //description should be something identifying the text source, and can be used
  //for debugging purposes or error message (examples would include a filename
  //or a database key).
  class TextInputStream {
  public:
    TextInputStream(const std::string& description);
    virtual ~TextInputStream(){}
    //True if getLine() can yield more data:
    virtual bool moreLines() const = 0;
    //if not endOfInput, replace contents of string with next line of text and
    //advance. Return false if no line was provided because input ran out.
    virtual bool getLine(std::string&) = 0;
    //Access source description
    const std::string& description() const { return m_descr; }
    virtual const char * streamType() const = 0;//e.g. "on-disk file", "memory buffer", ...
  private:
    std::string m_descr;
  };

  //Generic creation of text input streams from source-names. Such source names
  //would normally be a file-name (which will be automatically fed though the
  //find_file(..) function above), but it is possible to register a custom
  //TextInputManager (see below) which overrides this default behaviour. Such
  //customisation would features like alternative file searches, usage of
  //in-memory buffers, etc. Any code reading text files in NCrystal should
  //ultimately move to use this function in order to fully benefit from such
  //usage.
  void createTextInputStream( const std::string& sourcename,
                              UniquePtr<TextInputStream>& streamptr );
  //NB: We can't return unique ptrs by value in C++98 (no move semantics), which
  //    is why the function signature looks like it does.

  class TextInputManager {
  public:
    TextInputManager();
    virtual ~TextInputManager();
    //Reimplement this custom file searching in this function (can throw
    //FileNotFound in case of problems, but one will in any case be thrown if it
    //doesn't supply a result and if the fallback to the usual search patterns
    //is disallowed or fails):
    virtual void createTextInputStream( const std::string& sourcename,
                                        UniquePtr<TextInputStream>& streamptr ) = 0;
    //Override and return false to disable attempts to fall back to the usual
    //search for input files (if the method above supplies a null ptr).
    virtual bool allowFallbackToUsualDefaults() { return true; }
  };
  //Call to register a custom manager (assumes ownership - call with nullptr to
  //delete the custom manager before process shutdown):
  void registerTextInputManager(TextInputManager*);

  //Lower-level utilities for creating input streams directly (for convenience -
  //of course users can reimplemented their own TextInputStream class for even
  //more flexibility):
  void createTextInputStreamFromBuffer( const std::string& buffername,
                                        const std::string& buffer,
                                        UniquePtr<TextInputStream>& streamptr );
  void createTextInputStreamFromFile( const std::string& filepath,
                                      UniquePtr<TextInputStream>& streamptr );//NB: will NOT use the find_file(..) function




}

#endif

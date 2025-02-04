#ifndef NCrystal_FileUtils_hh
#define NCrystal_FileUtils_hh

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

// File system and path related utilities. We might one day make a path object,
// or simple wait for the glorious future when std::filesystem support is
// ubiquitous, but until then...

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  //Check if file exists and is readable:
  bool file_exists( const std::string& filename );

  //Simple file globbing, supporting wildcards only in the trailing directory or
  //file name (sorts results before returning):
  VectS ncglob( const std::string&);

  //Current working directory:
  std::string ncgetcwd();

  //Get dirname (returns "" for a path in the root dir, and "." for a filename
  //without a directory part):
  std::string dirname(const std::string& path);
  //Get the filename after the directory part
  std::string basename(const std::string& path);
  //The extension of the filename (if any):
  std::string getfileext(const std::string& path);
  //Normalise by analysing and reencoding:
  std::string normalise(const std::string& path);

  //Determine if path is absolute:
  bool path_is_absolute( const std::string& );

  //Join paths, trying to pick correct path separator (backwards slash on
  //Windows, unless path already contains forward slashes):
  std::string path_join(const std::string&, const std::string&);

  // tryRealPath: Wrap posix "realpath" or windows GetFullPathNameW function. On
  // any sort of error (including being un an unsupported platform), it returns
  // an empty string. Hopefully this function should work most of the time, and
  // fail gracefully in the rest.
  std::string tryRealPath( const std::string& path );

  //Read entire file into a string while protecting against someone mistakenly
  //trying to open a multi-gigabyte file and bringing their machine to a slow
  //halt. Will return NullOpt in case the file does not exists or is
  //unreadable. Calling on too large files will instead result in a
  //DataLoadError.
  Optional<std::string> readEntireFileToString( const std::string& path );

  //Use in applications (from main) needing to determine their own
  //paths. Returns empty string if unable.
  std::string determine_exe_self_path( int argc, char** argv );

}

#endif

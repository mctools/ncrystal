////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

#include "NCrystal/core/NCDefs.hh"
#include <fstream>

#if ( defined (_WIN32) || defined (WIN32) )
#  define NCRYSTAL_USE_WINDOWS_FILEUTILS

namespace NCRYSTAL_NAMESPACE {

  namespace WinFileUtils {
    //UTF-8 <-> UTF-16 conversion (via MultiByteToWideChar/WideCharToMultiByte
    //with CP_UTF8), needed whenever a native wide-string Win32 API must be
    //fed or read from NCrystal's usual UTF-8 strings. Exposed here (rather
    //than kept file-local) so other code in this component (e.g. the
    //GetEnvironmentVariableW-based env var lookup in NCString.cc) can reuse
    //it instead of duplicating the conversion:
    std::wstring winimpl_str2wstr( const std::string& );
    std::string winimpl_wstr2str( const std::wstring& );
    bool file_exists( const std::string& );
    std::ifstream open_ifstream_from_path( const std::string&,
                                           std::ios_base::openmode
                                           = std::ios_base::in);
    VectS ncglob_impl(const std::string&,
                      const std::string&);
    std::string get_current_working_dir();
    std::string get_self_exe_path_windows();
    std::string get_absolute_path(std::string);
  }
}

#endif

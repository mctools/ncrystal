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

#include "NCFileUtilsWin.hh"

#ifdef NCRYSTAL_USE_WINDOWS_FILEUTILS

#include "NCrystal/internal/NCString.hh"

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#ifdef _MSC_VER
//  Visual Studio adds non-standard std::ifstream(std::wstring) constructor:
#  define NCRYSTAL_FEATURE_IFSTREAM_WSTRING_PATH
#  ifndef NCRYSTAL_AVOID_CPP17FILESYSTEM
#    define NCRYSTAL_AVOID_CPP17FILESYSTEM //Not needed if we can use std::ifstream(std::wstring)
#  endif
#endif

#if (nc_cplusplus >= 201703L) && !defined(NCRYSTAL_AVOID_CPP17FILESYSTEM)//FIXME: Use feature testing macro as well!
#  include <filesystem>
#endif

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace WinFileUtils {

    namespace {
      std::wstring winimpl_str2wstr( const std::string& src )
      {
        const char * in_data = &src[0];
        int in_size = static_cast<int>(src.size());//fixme range check
        std::wstring res;
        if ( !in_size )
          return res;
        int out_size = MultiByteToWideChar( CP_UTF8,
                                            0,//Must be 0 for utf8
                                            in_data, in_size,
                                            nullptr, //dest buffer (nullptr for dry-run)
                                            0//dest buffer size (0 means
                                             //"dry-run" + return needed size).
                                            );
        const char * errmsg = "Failed to convert UTF-8 string to UTF-16";
        if (!out_size)
          NCRYSTAL_THROW(BadInput,errmsg);
        res.resize( out_size );
        wchar_t * out_data = &res[0];
        //Same, but with out_data/out_size provided:
        out_size = MultiByteToWideChar( CP_UTF8, 0,
                                        in_data, in_size,
                                        out_data, out_size );
        if ( out_size != res.size() )
          NCRYSTAL_THROW(BadInput,errmsg);
        return res;
      }

      std::string winimpl_wstr2str( const std::wstring& src )
      {
        const wchar_t * in_data = &src[0];
        int in_size = static_cast<int>(src.size());//fixme range check
        std::string res;
        if ( !in_size )
          return res;
        int out_size = WideCharToMultiByte( CP_UTF8,
                                            0,//Must be 0 for utf8
                                            in_data, in_size,
                                            nullptr, //dest buffer (nullptr for dry-run)
                                            0,//dest buffer size (0 means "dry-run"
                                            //returning, needed size)
                                            nullptr,//Must be null for utf8
                                            nullptr//Must be null for utf8
                                            );
        const char * errmsg = "Failed to convert UTF-16 string to UTF-8";
        if (!out_size)
          NCRYSTAL_THROW(BadInput,errmsg);
        res.resize( out_size );
        char * out_data = &res[0];
        //Same, but with out_data/out_size provided:
        out_size = WideCharToMultiByte( CP_UTF8, 0,
                                        in_data, in_size,
                                        out_data, out_size,
                                        nullptr, nullptr);
        if ( out_size != res.size() )
          NCRYSTAL_THROW(BadInput,errmsg);
        return res;
      }
    }

    bool file_exists( const std::string& path )
    {
      if ( isSimpleASCII( path ) ) {
        std::ifstream f(path.c_str());
        return f.good();
      } else {
        //Need UTF-16 API:
        auto wpath = winimpl_str2wstr( path );
        return (_waccess(wpath.c_str(), 0) == 0);
      }
    }

    std::ifstream open_ifstream_from_path( const std::string& path,
                                           std::ios_base::openmode mode )
    {
      if ( isSimpleASCII( path ) )
        return std::ifstream(path, mode);

      //Need UTF-16 API:
      auto wpath = winimpl_str2wstr( path );

#ifdef NCRYSTAL_FEATURE_IFSTREAM_WSTRING_PATH
      return std::ifstream( wpath, mode );
#elif (nc_cplusplus >= 201703L) && !defined(NCRYSTAL_AVOID_CPP17FILESYSTEM)
      //Rely on presence of C++17 filesystem.
      auto fspath = std::filesystem::path{ path }
      return std::ifstream( fspath, mode );
#else
      NCRYSTAL_THROW(BadInput,"This installation of NCrystal does not support"
                     " non-ASCII filesystem paths on Windows");
      return std::ifstream();
#endif
    }

    std::string get_current_working_dir()
    {
      std::wstring wpath;
      wpath.resize(MAX_PATH);
      auto nsize = GetCurrentDirectoryW(wpath.size(), &wpath[0]);
      nc_assert_always(nsize<=wpath.size());
      wpath.resize(nsize);
      return winimpl_wstr2str( std::move(wpath) );
    }

    VectS ncglob( const std::string& pattern_utf8 )
    {
      auto wpattern = winimpl_str2wstr( pattern_utf8 );
      VectS result;
      WIN32_FIND_DATAW fdata;
      HANDLE fh = FindFirstFileW(wpattern.c_str(), &fdata);
      if (fh == INVALID_HANDLE_VALUE)
        return result;
      while (true) {
        std::wstring hitw( fdata.cFileName );
        std::string hit_utf8 = winimpl_wstr2str( hitw );
        result.push_back(hit_utf8);
        if (!FindNextFileW(fh, &fdata))
          break;
      }
      FindClose(fh);
      result.shrink_to_fit();
      std::sort(result.begin(),result.end());
      return result;
    }

  }
}
#endif

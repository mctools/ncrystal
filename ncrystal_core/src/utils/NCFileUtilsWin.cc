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

#include "NCFileUtilsWin.hh"

#ifdef NCRYSTAL_USE_WINDOWS_FILEUTILS

#include "NCrystal/internal/utils/NCString.hh"

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#ifdef _MSC_VER
//  Visual Studio adds non-standard std::ifstream(std::wstring) constructor:
#  define NCRYSTAL_FEATURE_IFSTREAM_WSTRING_PATH
#endif

//Do not use std::filesystem if we have an alternative:
#ifndef NCRYSTAL_AVOID_CPP17FILESYSTEM
#  ifdef NCRYSTAL_FEATURE_IFSTREAM_WSTRING_PATH
#    define NCRYSTAL_AVOID_CPP17FILESYSTEM
#  endif
#endif

//Do not use std::filesystem if not present:
#ifndef NCRYSTAL_AVOID_CPP17FILESYSTEM
// NB: The __cpp_lib_filesystem macro might only exist in C++20.
#  if nc_cplusplus < 201703L || !defined(__cpp_lib_filesystem)
#    define NCRYSTAL_AVOID_CPP17FILESYSTEM
#  endif
#endif

#ifndef NCRYSTAL_AVOID_CPP17FILESYSTEM
#  include <filesystem>
#endif
#include <climits>

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace WinFileUtils {

    namespace {
      std::wstring winimpl_str2wstr( const std::string& src )
      {
        const char * in_data = &src[0];
        nc_assert_always( src.size() < INT_MAX );
        int in_size = static_cast<int>(src.size());
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
        nc_assert_always( (std::size_t) src.size() < (std::size_t)INT_MAX );
        int in_size = static_cast<int>(src.size());
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
        if ( out_size != static_cast<int>(res.size()) )
          NCRYSTAL_THROW(BadInput,errmsg);
        return res;
      }
    }

    bool file_exists( const std::string& path )
    {
#if 0
      //We could use the dedicated windows API:
      auto wpath = winimpl_str2wstr( path );
      return (_waccess(wpath.c_str(), 0) == 0);
#endif
      //But it might be easier to support using the open_ifstream_from_path
      //function since we already have that implemented:
      std::ifstream f = open_ifstream_from_path( path,
                                                 std::ios_base::in );
      return f.good();
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
#elif !defined(NCRYSTAL_AVOID_CPP17FILESYSTEM)
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

    std::wstring path2wpath( const std::string& path )
    {
      if ( contains( path, '/' ) ) {
        std::string pfix = path;
        for ( auto& c : pfix )
          if ( c == '/' )
            c = '\\';
        return winimpl_str2wstr( pfix );
      } else {
        return winimpl_str2wstr( path );
      }
    }

    std::string get_absolute_path(std::string path)
    {
      if ( path.empty() )
        return {};
      auto wpath = path2wpath( path );

      auto len_with_null_term = GetFullPathNameW( &wpath[0],
                                                  0, nullptr, nullptr );
      if ( len_with_null_term <= 1 )
        return {};
      std::wstring wres;
      wres.resize(len_with_null_term-1,0);
      auto len = GetFullPathNameW( &wpath[0],
                                   len_with_null_term, &wres[0],
                                   nullptr );
      if ( len != wres.size() )
        return {};//failed

      return winimpl_wstr2str( wres );
    }

    VectS ncglob_impl( const std::string& pattern_utf8 )
    {
      std::wstring wpattern;
      if ( contains( pattern_utf8, '/' ) ) {
        std::string pfix = pattern_utf8;
        for ( auto& c : pfix )
          if ( c == '/' )
            c = '\\';
        wpattern = winimpl_str2wstr( pfix );
      } else {
        wpattern = winimpl_str2wstr( pattern_utf8 );
      }
      VectS result;
      WIN32_FIND_DATAW fdata;
      constexpr auto dwAdditionalFlags = ( FIND_FIRST_EX_CASE_SENSITIVE
                                           | FIND_FIRST_EX_LARGE_FETCH );
      HANDLE fh = FindFirstFileExW( wpattern.c_str(),
                                    FindExInfoBasic,
                                    &fdata,
                                    FindExSearchNameMatch,
                                    nullptr,
                                    dwAdditionalFlags );
      if (fh == INVALID_HANDLE_VALUE) {
        auto last_error = GetLastError();
        SetLastError(0);//clear
        if ( last_error == ERROR_FILE_NOT_FOUND ) {
          //Not actually an error condition, just no hits.
        } else if ( last_error == ERROR_INVALID_NAME ) {
          NCRYSTAL_THROW2(BadInput,"Invalid glob pattern: \""
                          <<pattern_utf8<<"\"");
        } else {
          NCRYSTAL_THROW2(CalcError,"Unexpected error (code "<<last_error
                          <<") while globbing via windows"
                          " FindFirstFileW function");
        }
        return result;
      }
      while (true) {
        std::wstring hitw( fdata.cFileName );
        std::string hit_utf8 = winimpl_wstr2str( hitw );
        if ( !hit_utf8.empty() && hit_utf8 != "." && hit_utf8!=".." )
          result.push_back(hit_utf8);
        if (!FindNextFileW(fh, &fdata))
          break;
      }
      FindClose(fh);
      return result;
    }

    std::string get_self_exe_path_windows()
    {
      std::wstring wpath;
      wpath.resize(MAX_PATH);
      auto nsize = GetModuleFileNameW(nullptr, &wpath[0], MAX_PATH);
      nc_assert_always(nsize<=wpath.size());
      wpath.resize(nsize);
      return winimpl_wstr2str( std::move(wpath) );
    }
  }
}
#endif

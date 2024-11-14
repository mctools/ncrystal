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

#include "NCrystal/internal/NCCFileUtils.hh"

// Cross platform string and file system utilities for C.

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////                                                                       ////
////                     C-compatible string utilities                     ////
////                                                                       ////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


#ifdef MCFILEUTILS_CPPNAMESPACE
#  include <cstdlib>
#  include <cstring>
#  include <cassert>
namespace MCFILEUTILS_CPPNAMESPACE {
#else
#  include <stdlib.h>
#  include <string.h>
#  include <assert.h>
#endif

#ifdef STDNS
#  undef STDNS
#endif
#ifdef MCFILEUTILS_CPPNAMESPACE
#  define STDNS std::
#else
#  define STDNS
#endif

  mcu8str mcu8str_create_empty()
  {
    static char dummy[4] = { 0, 0, 0, 0 };
    mcu8str s;
    s.c_str = &dummy[0];
    s.size = 0;
    s.buflen = 0;
    s.owns_memory = 0;
    return s;
  }

  void mcu8str_dealloc( mcu8str* str )
  {
    char * buf_to_free = ( str->owns_memory
                           ? str->c_str
                           : (char*)(0) );
    static char dummy[4] = { 0, 0, 0, 0 };
    str->c_str = &dummy[0];
    str->size = 0;
    str->buflen = 0;
    str->owns_memory = 0;
    if ( buf_to_free )
      STDNS free( buf_to_free );
  }

  mcu8str mcu8str_create( size_t prealloc_size )
  {
    if ( prealloc_size == 0 )
      return mcu8str_create_empty();
    mcu8str s;
    s.c_str = (char*) STDNS malloc( prealloc_size+1 );
    if ( !s.c_str ) {
      //Fixme: exception? (if c++). Otherwise printout?
      return mcu8str_create_empty();
    }
    s.c_str[0] = '\0';
    s.size = 0;
    s.buflen = prealloc_size + 1;
    s.owns_memory = 1;
    return s;
  }

  void mcu8str_append( mcu8str* str, const mcu8str * otherstr )
  {
    if ( otherstr->size == 0 )
      return;
    size_t newsize = str->size + otherstr->size;
    if ( newsize + 1 > str->buflen )
      mcu8str_reserve( str, newsize );
    //memcpy ok since we handle the null char afterwards, there will be no
    //overlapping even if str==otherstr:
    STDNS memcpy( str->c_str + str->size,
                  otherstr->c_str,
                  otherstr->size );
    str->c_str[newsize] = '\0';
    str->size = newsize;
  }

  void mcu8str_append_cstr( mcu8str* str, const char * c_str )
  {
    const size_t o_size = STDNS strlen( c_str );
    if ( o_size == 0 )
      return;
    const size_t newsize = str->size + o_size;
    if ( newsize + 1 > str->buflen )
      mcu8str_reserve( str, newsize );
    //memcpy ok since we handle the null char afterwards, there will be no
    //overlapping even if c_str = str.c_str:
    STDNS memcpy( str->c_str + str->size, c_str, o_size );
    str->c_str[newsize] = '\0';
    str->size = newsize;
  }

  mcu8str mcu8str_create_from_cstr( const char * c_str )
  {
    if ( ! *c_str )
      return mcu8str_create_empty();
    const size_t o_size = STDNS strlen( c_str );
    assert( o_size > 0 );
    mcu8str str = mcu8str_create( o_size );
    //str.c_str is newly allocated, so there can be no overlap:
    STDNS memcpy( str.c_str, c_str, o_size + 1 );
    str.size = o_size;
    return str;
  }

  const mcu8str mcu8str_view_cstr( const char * c_str )
  {
    mcu8str s;
    s.c_str = (char*)c_str;//cast away constness, but we will add it again in
                           //the return type.
    s.size = STDNS strlen( c_str );
    s.buflen = s.size + 1;
    s.owns_memory = 0;
    return s;
  }

  void mcu8str_swap( mcu8str* str1, mcu8str* str2 )
  {
    mcu8str tmp;
    tmp.c_str = str1->c_str;
    tmp.size = str1->size;
    tmp.buflen = str1->buflen;
    tmp.owns_memory = str1->owns_memory;
    str1->c_str = str2->c_str;
    str1->size = str2->size;
    str1->buflen = str2->buflen;
    str1->owns_memory = str2->owns_memory;
    str2->c_str = tmp.c_str;
    str2->size = tmp.size;
    str2->buflen = tmp.buflen;
    str2->owns_memory = tmp.owns_memory;
  }

  void mcu8str_update_size( mcu8str* str )
  {
    if ( str->buflen < 2 ) {
      //This handles not only empty strings, but also the special case of
      //buflen=0, size=0 which occurs after _create_empty or _deallocate:
      str->size = 0;
      assert( str->c_str[0] == '\0' );
      return;
    }
    str->size = STDNS strlen( str->c_str );
    assert( str->size < str->buflen );
  }

  void mcu8str_assign( mcu8str* dest, const mcu8str* src )
  {
    if ( src->size +1 <= dest->buflen ) {
      //Already fits. Using memmove in case of self assignment:
      STDNS memmove( dest->c_str, src->c_str, src->size + 1 );
      dest->size = src->size;
      return;
    }
    mcu8str newstr = mcu8str_create( src->size );
    mcu8str_assign( &newstr, src );
    mcu8str_swap( dest, &newstr );
    mcu8str_dealloc( &newstr );
  }

  mcu8str mcu8str_copy( const mcu8str* str )
  {
    if ( str->size == 0 )
      return mcu8str_create_empty();
    mcu8str newstr = mcu8str_create( str->size );
    mcu8str_assign( &newstr, str );
    return newstr;
  }

  void mcu8str_reserve( mcu8str* str, size_t nsize )
  {
    if ( nsize +1 <= str->buflen )
      return;//already has enough
    mcu8str newstr = mcu8str_create( nsize );
    //newstr.c_str can not overlap since it was just malloc'ed:
    STDNS memcpy( newstr.c_str, str->c_str, str->size + 1 );
    newstr.size = str->size;
    mcu8str_swap( str, &newstr );
    mcu8str_dealloc( &newstr );
  }
  void mcu8str_clear( mcu8str* str )
  {
    if ( str->size ) {
      str->c_str[0] = '\0';
      str->size = 0;
    }
  }

  void mcu8str_ensure_dynamic_buffer( mcu8str* str )
  {
    if ( str->owns_memory )
      return;
    mcu8str str2 = mcu8str_copy( str );
    mcu8str_swap( &str2, str );
    assert( !str2.owns_memory );
    assert( str->owns_memory );
    //No need to call mcu8str_dealloc( &str2 ) since it does not own the memory
    //now.
  }

  //Allow usage of static buffers:
  mcu8str mcu8str_create_from_staticbuffer( char * buf, size_t buflen )
  {
    assert(buflen > 0);
    mcu8str s;
    s.owns_memory = 0;
    s.c_str = buf;
    s.c_str[0] = '\0';
    s.size = 0;
    s.buflen = buflen;
    return s;
  }

  int mcu8str_is_ascii( const mcu8str* str )
  {
    const char * it = str->c_str;
    const char * itE = it + str->size;
    for ( ; it!=itE; ++it )
      if ( (unsigned char)(*it) >= 128 )
        return 0;
    return 1;
  }

  int mcu8str_contains( const mcu8str* str, char needle )
  {
    const char * it = str->c_str;
    const char * itE = it + str->size;
    for ( ; it!=itE; ++it )
      if ( *it == needle )
        return 1;
    return 0;
  }

  void mcu8str_replace( mcu8str* str, char from, char to )
  {
    char * it = str->c_str;
    char * itE = it + str->size;
    for ( ; it!=itE; ++it )
      if ( *it == from )
        *it = to;
  }

#ifdef MCFILEUTILS_CPPNAMESPACE
  std::string mcu8str_to_stdstring( const mcu8str* str )
  {
    return ( str->size
             ? std::string( str->c_str, str->size )
             : std::string() );
  }
#endif
#ifdef MCFILEUTILS_CPPNAMESPACE
}
#endif

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////                                                                       ////
////                         File system utilities                         ////
////                                                                       ////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

#ifdef MCFILEUTILS_CPPNAMESPACE
#  include <cstdlib>
#  include <cstring>
#  include <cassert>
#  include <cstdint>
namespace MCFILEUTILS_CPPNAMESPACE {
#else
#  include <stdlib.h>
#  include <string.h>
#  include <assert.h>
#  include <stdint.h>
#endif
#ifdef MC_IS_WINDOWS
#  undef MC_IS_WINDOWS
#endif
#ifdef MC_IS_APPLE
#  undef MC_IS_APPLE
#endif
#if ( defined (_WIN32) || defined (WIN32) )
#  define MC_IS_WINDOWS
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
  //#  include <shlwapi.h>//for PathIsRelativeW (fixme: remove if needed)
#else
#  include <unistd.h>
#  include <limits.h>
#  include <errno.h>
#  include <sys/stat.h>
#  ifdef __APPLE__
#    define MC_IS_APPLE
#    include <mach-o/dyld.h>
#  endif
#endif

#if 0
  //Fake being on Windows to check basic compilation
#  define MC_IS_WINDOWS
#  define CP_UTF8 1
#  define MAX_PATH 260
#  define DWORD unsigned
#  define INVALID_FILE_ATTRIBUTES 0x123
#  define FILE_ATTRIBUTE_DIRECTORY 0x1000
  int MultiByteToWideChar( int,int,const char*,int,wchar_t*,int);
  int WideCharToMultiByte( int, int,const wchar_t*,int,char*,int,void*,void*);
  int GetCurrentDirectoryW( int, wchar_t* );
  int GetModuleFileNameW( void*, wchar_t*, int );
  DWORD GetFullPathNameW( const wchar_t*, DWORD, wchar_t*, wchar_t** );
  DWORD GetFileAttributesW( const wchar_t* );
  FILE * _wfopen( const wchar_t*,const wchar_t*);
  //int PathIsRelativeW( const wchar_t* );
#endif

#if 0
  //Fake being on macos to check basic compilation
#  define MC_IS_APPLE
  int _NSGetExecutablePath(char*, uint32_t*);
#endif


#ifdef MC_IS_WINDOWS
  namespace {
    //Need wide-char strings for unicode API. We try to be less fancy with these
    //than the mcu8str's, since they are anyway only used on Windows and in this
    //file.
    struct mcwinstr {
      wchar_t * c_str;
      STDNS size_t size;
      STDNS size_t buflen;
    };
    mcwinstr mc_winstr_create( STDNS size_t size );
    mcwinstr mc_winstr_create_empty();
    mcwinstr mc_u8str_to_winstr( const mcu8str* src )
    {
      const int in_size = (int)( src->size );//fixme range check

      //First check for empty string:
      if ( !src->c_str || in_size == 0 )
        return mc_winstr_create_empty();

      const char * in_data = src->c_str;
      //fixme:nc_assert_always( src.size() < INT_MAX );
      int out_size = MultiByteToWideChar( CP_UTF8,
                                          0,//Must be 0 for utf8
                                          in_data, in_size,
                                          nullptr, //dest buffer (nullptr for dry-run)
                                          0//dest buffer size (0 means
                                          //"dry-run" + return needed size).
                                          );
      const char * errmsg = "Failed to convert UTF-8 string to UTF-16";
      (void)errmsg;//fixme
      if (!out_size)
        return mc_winstr_create_empty();//fixme: NCRYSTAL_THROW(BadInput,errmsg);

      mcwinstr res = mc_winstr_create( (STDNS size_t)( out_size ) );

      wchar_t * out_data = res.c_str;
      //Same, but with out_data/out_size provided:
      int out_size2 = MultiByteToWideChar( CP_UTF8, 0,
                                           in_data, in_size,
                                           out_data, out_size );
      if ( out_size != out_size2 || (STDNS size_t)out_size >= res.buflen )
        return mc_winstr_create_empty();//fixme: NCRYSTAL_THROW(BadInput,errmsg);
      res.c_str[out_size] = 0;
      res.size = out_size;
      return res;
    }

    mcwinstr mc_path2wpath( const mcu8str* src )
    {
      //Like mc_u8str_to_winstr but maps any '/' to '\' and upper cases drive
      //letters:
      mcwinstr res;
      char drive_letter = mctools_drive_letter( src );
      if ( drive_letter >= 'a' && drive_letter <= 'z' ) {
        //Let us uppercase it for consistency. For simplicity we do it in the
        //utf8 string at the cost of another malloc/free:
        mcu8str dummy = mcu8str_copy(src);
        dummy.c_str[0] = drive_letter;
        res = mc_u8str_to_winstr( &dummy );
        mcu8str_dealloc(&dummy);
      } else {
        res = mc_u8str_to_winstr( src );
      }
      //Convert to native path:
      wchar_t * it = res.c_str;
      wchar_t * itE = it + res.size;
      for ( ; it!=itE; ++it )
        if ( *it == L'/' )
          *it = L'\\';
      return res;
    }

    mcu8str mc_winstr_to_u8str( const mcwinstr* src )
    {
      const int in_size = (int)( src->size );//fixme range check
      //First check for empty string (also for safeguard in case of deallocated
      //string):
      if ( !src->c_str || in_size == 0 )
        return mcu8str_create_empty();

      const wchar_t * in_data = src->c_str;

      int out_size = WideCharToMultiByte( CP_UTF8,
                                          0,//Must be 0 for utf8
                                          in_data, in_size,
                                          nullptr, //dest buffer (nullptr for dry-run)
                                          0,//dest buffer size (0 means "dry-run"
                                            //returning, needed size)
                                          nullptr,//Must be null for utf8
                                          nullptr//Must be null for utf8
                                          );
      //NB: out_size does not include 0 termination char!
      const char * errmsg = "Failed to convert UTF-16 string to UTF-8";
      if (!out_size)
        return mcu8str_create_from_cstr(errmsg);//fixme: NCRYSTAL_THROW(BadInput,errmsg);
      mcu8str res = mcu8str_create( out_size );
      //res.resize( out_size );
      char * out_data = res.c_str;
      //Same, but with out_data/out_size provided:
      int out_size2 = WideCharToMultiByte( CP_UTF8, 0,
                                           in_data, in_size,
                                           out_data, out_size,
                                           nullptr, nullptr);
      if ( out_size2 != out_size || (STDNS size_t)out_size >= res.buflen )
        return mcu8str_create_from_cstr(errmsg);//fixme: NCRYSTAL_THROW(BadInput,errmsg);
      res.c_str[out_size] = 0;
      res.size = out_size;
      return res;
    }
    mcwinstr mc_winstr_create_empty()
    {
      //empty (keep in allocated buffer for simplicity).
      mcwinstr res_empty = mc_winstr_create(1);
      res_empty.c_str[0] = 0;
      res_empty.size = 0;
      return res_empty;
    }
    mcwinstr mc_winstr_create( STDNS size_t size )
    {
      mcwinstr str;
      str.c_str = (wchar_t*) STDNS malloc( sizeof(wchar_t)*(size + 1) );
      str.c_str[0] = 0;
      str.size = 0;
      str.buflen = size;
      return str;
    }
    void mc_winstr_dealloc( mcwinstr* str )
    {
      if ( !str->c_str )
        return;
      void * bufptr = (void*) str->c_str;
      str->size = 0;
      str->buflen = 0;
      str->c_str = (wchar_t*)(0);
      STDNS free( bufptr );
    }



  }
#endif

  MCTOOLS_FILE_t * mctools_fopen( const mcu8str* path,
                                  const char * mode )
  {
    //Utf8-encoded paths are used directly on platforms except Windows, where
    //only ASCII encoded paths can be used directly:
#ifdef MC_IS_WINDOWS
    const int can_use_directly = ( mcu8str_is_ascii(path) ? 1 : 0 );
    const char nonnative_sep = '/';
#else
    const int can_use_directly = 1;
    const char nonnative_sep = '\\';
#endif
    MCTOOLS_FILE_t * result = (MCTOOLS_FILE_t *)0;
    if ( can_use_directly ) {
      if ( mcu8str_contains( path, nonnative_sep ) ) {
        //Support non-native path seps:
        char buf[4096];
        mcu8str native = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
        mcu8str_assign( &native, path );
        mctools_pathseps_platform( &native );
        result = STDNS fopen( native.c_str, mode );
        mcu8str_dealloc(&native);
      } else {
        result = STDNS fopen( path->c_str, mode );
      }
    }
#ifdef MC_IS_WINDOWS
    else {
      //Ok, non-ASCII path on windows. Reencode as wide strings and use _wfopen:
      mcwinstr wpath = mc_path2wpath( path );//this also maps to native path sep
      char buf[32];
      mcu8str mcu8str_mode = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
      mcu8str_append_cstr(&mcu8str_mode,mode);
      mcwinstr wmode = mc_u8str_to_winstr( &mcu8str_mode );
      result = _wfopen( wpath.c_str, wmode.c_str );//nb: no STDNS here on purpose
      mcu8str_dealloc(&mcu8str_mode);
      mc_winstr_dealloc( &wpath );
      mc_winstr_dealloc( &wmode );
    }
#endif
    return result;
  }

  int mctools_file_exists_and_readable( const mcu8str* path )
  {
    MCTOOLS_FILE_t * fh = mctools_fopen( path, "r" );
    if ( !fh )
      return 0;
    fclose(fh);
    return 1;
  }

  mcu8str mctools_get_current_working_dir()
  {
#ifdef MC_IS_WINDOWS
    mcwinstr wpath = mc_winstr_create( ( MAX_PATH >= 260 ? MAX_PATH : 260 ) );
    auto nsize = GetCurrentDirectoryW(wpath.buflen,wpath.c_str);
    if ( (STDNS size_t)(nsize +1) > wpath.buflen ) {
      //Use larger buffer and try again:
      mc_winstr_dealloc( &wpath );
      wpath = mc_winstr_create( nsize );
      nsize = GetCurrentDirectoryW(wpath.buflen,wpath.c_str);
    }
    if ( nsize == 0 || (STDNS size_t )( nsize +1 ) > wpath.buflen ) {
      //Error (fixme):
      mc_winstr_dealloc( &wpath );//check that we dealloc correctly everywhere
                                  //in this file in case of errors
      return mcu8str_create_empty();
    }
    wpath.c_str[nsize] = 0;
    wpath.size = nsize = 0;
    mcu8str res = mc_winstr_to_u8str( &wpath );
    mc_winstr_dealloc( &wpath );
    return res;
#else
    char buf[4096];//almost always enough in first go (fixme: test with tiny value here)
    mcu8str res = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
    while ( 1 ) {
      if ( getcwd( res.c_str, res.buflen ) ) {
        mcu8str_update_size( &res );
        mcu8str_ensure_dynamic_buffer(&res);//don't return local static buffer
        return res;
      }
      if ( errno == ERANGE ) {
        //Try again with larger buffer: (fixme do not do this forever)
        mcu8str_clear( &res );
        mcu8str_reserve( &res, res.buflen * 2 );
        continue;
      } else {
        //FIXME: ERROR!!
        return mcu8str_create_empty();//fixme
      }
    }
    //fixme: some of these? or perhaps use them in the calling code??? It might
    // be easiest if we return error codes and let the calling code deal with
    // error emission.

    //NCRYSTAL_THROW(CalcError,"current working directory is too
    // long"); NCRYSTAL_THROW(CalcError,"Could not determine current working
    // directory");
#endif
  }

  mcu8str mctools_determine_exe_self_path( int argc, char** argv )
  {
    //First try some platform specific methods not relying on argv0:
#ifdef MC_IS_WINDOWS
    {
      mcwinstr wpath = mc_winstr_create( (STDNS size_t)(MAX_PATH<4096
                                                        ?4096:MAX_PATH) +1 );
      auto nsize = GetModuleFileNameW(nullptr, wpath.c_str, wpath.buflen-1 );
      assert( (STDNS size_t)(nsize) <= wpath.size );
      wpath.c_str[nsize] = 0;//add null terminator
      wpath.size = nsize;
      if ( nsize > 0 ) {
        mcu8str res = mc_winstr_to_u8str( &wpath );
        mc_winstr_dealloc( &wpath );
        return res;
      }
    }
#endif
#ifdef MC_IS_APPLE
    {
      mcu8str path = mcu8str_create( 4096 );//fixme try with tiny value here
      STDNS uint32_t bufsize = path.buflen;
      int status = _NSGetExecutablePath( path.c_str, &bufsize );
      if ( status == -1 ) {
        //buffer too small, try again
        mcu8str_reserve( &path, bufsize+1 );
        bufsize = path.buflen;
        status = _NSGetExecutablePath( path.c_str, &bufsize );
      }
      if ( status == 0 ) {
        //success! Return (but make sure we do not refer to our static buffer)
        mcu8str_update_size( &path );
        if ( path.size > 0 )
          return path;
      }
      mcu8str_dealloc(&path);
    }
#endif
#if defined(__unix__) && !defined(MC_IS_APPLE) && !defined(MC_IS_WINDOWS)
    {
      //Try to invoke readlink on relevant file in /proc. Investigate both
      //"/proc/self/exe" (linux) and "/proc/curproc/file" (might work on some
      //BSDs).
      for ( int ifile = 0; ifile < 2; ++ifile ) {
        const char * filename = ( ifile==0
                                  ? "/proc/self/exe"
                                  : "/proc/curproc/file" );

        char buf[65536+1];//PATH_MAX is unreliable so we use huge buffer for
        //simplicity.
        mcu8str path = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
        //NB: ssize_t is from unistd.h, so never exist in std:: namespace.
        ssize_t len = readlink(filename, path.c_str, path.buflen-1 );
        if ( len > 0 && (STDNS size_t)(len+1) < path.buflen ) {
          path.c_str[len] = '\0';//readlink does not add terminating null char
          path.size = (STDNS size_t)len;
          mcu8str_ensure_dynamic_buffer( &path );
          return path;
        }
        mcu8str_dealloc(&path);
      }
    }
#endif
    //Fall back to argv0 analysis. This will probably only happen outside of the
    //three major platforms (win/macos/linux), so most likely this means we are
    //on some sort of BSD variant. If argv[0] is an absolute path and the file
    //at that path exists, we can take that as a somewhat reliable guess:

    if ( !(argc > 0) )
      return mcu8str_create_empty();//not available (rare, but allowed by C std)

    if ( argv[0][0] == '/' ) {
      //Looks like an absolute path:
      mcu8str path = mcu8str_create_from_cstr( argv[0] );
      if ( mctools_file_exists_and_readable( &path ) )
        return path;
      mcu8str_dealloc( &path );
    }

    //We could do more extensive analysis here, but it might give wrong results
    //if we are unlucky.
    return mcu8str_create_empty();
  }

  mcu8str mctools_real_path( const mcu8str* path )
  {
    if ( path->size == 0 )
      return mcu8str_create_empty();
#ifdef MC_IS_WINDOWS
    mcwinstr wpath = mc_path2wpath( path );
    DWORD len_with_null_term = GetFullPathNameW( wpath.c_str,
                                                 0, nullptr, nullptr );
    if ( len_with_null_term <= 1 ) {
      mc_winstr_dealloc( &wpath );
      return mcu8str_create_empty();
    }
    mcwinstr woutput = mc_winstr_create((STDNS size_t)(len_with_null_term-1));
    DWORD len = GetFullPathNameW( wpath.c_str,
                                  woutput.buflen, woutput.c_str,
                                  nullptr );
    mc_winstr_dealloc( &wpath );
    if ( len+1 != len_with_null_term ) {
      //Failure:
      mc_winstr_dealloc( &woutput );
      return mcu8str_create_empty();
    }
    mcu8str output = mc_winstr_to_u8str( &woutput );
    mc_winstr_dealloc( &woutput );
#else
    char buf[4096];
    mcu8str native = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
    mcu8str_assign( &native, path );
    mctools_pathseps_platform( &native );
    mcu8str output = mcu8str_create( PATH_MAX );
    if ( !realpath( native.c_str, output.c_str ) ) {
      //failure
      mcu8str_dealloc(&output);
      mcu8str_dealloc(&native);
      return mcu8str_create_empty();
    }
    mcu8str_dealloc(&native);

#endif
    mctools_pathseps_platform(&output);
    return output;
  }

  void mctools_pathseps_generic( mcu8str* path )
  {
    mcu8str_replace( path, '\\', '/' );
    //Upper case drive letter for consistency:
    char drive_letter = mctools_drive_letter( path );
    if ( drive_letter )
      path->c_str[0] = drive_letter;

  }

  void mctools_pathseps_platform( mcu8str* path )
  {
#ifdef MC_IS_WINDOWS
    mcu8str_replace( path, '/', '\\' );
#else
    mcu8str_replace( path, '\\', '/' );
#endif
    //Upper case drive letter for consistency:
    char drive_letter = mctools_drive_letter( path );
    if ( drive_letter )
      path->c_str[0] = drive_letter;
  }

  int mctools_path_is_relative( const mcu8str* path )
  {
    return ( mctools_path_is_absolute(path) ? 0 : 1 );
  }

  int mctools_path_is_absolute( const mcu8str* path )
  {
    if ( path->size == 0 )
      return 0;
    const char * path_begin = path->c_str;
    if ( mctools_drive_letter( path ) )
      path_begin += 2;
    if ( *path_begin == '/' || *path_begin == '\\' )
      return 1;
    return 0;
// #ifdef MC_IS_WINDOWS
//     //FIXME: Would it perhaps be better if we implemented this function WITHOUT
//     //platform specific behaviour???
//     mcwinstr wpath = mc_path2wpath( path );
//     int is_abs = ( PathIsRelativeW(wpath.c_str) ? 0 : 1 );
//     mc_winstr_dealloc( &wpath );
//     return is_abs;
// #else
//     return 0;
// #endif
  }

  char mctools_drive_letter( const mcu8str* path )
  {
    if ( path->size < 2 )
      return 0;
    if ( path->c_str[1] != ':' )
      return 0;
    if ( path->c_str[0] >= 'A' && path->c_str[0] <= 'Z' )
      return path->c_str[0];
    if ( path->c_str[0] >= 'a' && path->c_str[0] <= 'z' )
      return (char)(path->c_str[0] - 32);//uppercase
    return 0;
  }

  int mcu8str_equal( const mcu8str* s1, const mcu8str* s2 )
  {
    return ( s1->size == s2->size
             && STDNS memcmp( s1->c_str, s2->c_str, s1->size )==0 );
  }

  int mctools_is_dir( const mcu8str* path )
  {
#ifdef MC_IS_WINDOWS
    mcwinstr wpath = mc_path2wpath( path );
    DWORD fa = GetFileAttributesW( wpath.c_str );
    mc_winstr_dealloc( &wpath );
    if ( fa == INVALID_FILE_ATTRIBUTES )
      return 0;
    return ( fa & FILE_ATTRIBUTE_DIRECTORY) ? 1 : 0;
#else
    char buf[4096];
    mcu8str native = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
    mcu8str_assign( &native, path );
    mctools_pathseps_platform( &native );
    struct stat sinfo;
    int res = ( stat( native.c_str, &sinfo) == 0
                && S_ISDIR(sinfo.st_mode) ) ? 1 : 0;
    mcu8str_dealloc(&native);
    return res;
#endif
  }

  int mctools_is_same_file( const mcu8str* p1, const mcu8str* p2 )
  {
    //First check for trivial string equality:
    if ( mcu8str_equal( p1, p2 ) )
      return mctools_file_exists_and_readable( p1 );

#ifdef MC_IS_WINDOWS
    //For now, just compare the normalised paths (FIXME?):
    mcu8str rp1 = mctools_real_path( p1 );
    mcu8str rp2 = mctools_real_path( p2 );
    if ( mcu8str_equal( &rp1, &rp2 ) ) {
      mcu8str_dealloc( &rp1 );
      mcu8str_dealloc( &rp2 );
      return 1;
    }
    return 0;
#else
    //A hopefully bullet proof way, (st_ino,st_dev) uniquely identifies a file
    //on a POSIX system.
    struct stat sinfo1;
    struct stat sinfo2;
    char buf[4096];
    {
      mcu8str native = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
      mcu8str_assign( &native, p1 );
      mctools_pathseps_platform( &native );
      FILE * fh = fopen(native.c_str,"rb");
      mcu8str_dealloc(&native);
      if ( !fh )
        return 0;
      if ( fstat(fileno(fh), &sinfo1) != 0 ) {
        fclose(fh);
        return 0;
      }
      fclose(fh);
    }
    {
      mcu8str native = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
      mcu8str_assign( &native, p2 );
      mctools_pathseps_platform( &native );
      FILE * fh = fopen(native.c_str,"rb");
      mcu8str_dealloc(&native);
      if ( !fh )
        return 0;
      if ( fstat(fileno(fh), &sinfo2) != 0 ) {
        fclose(fh);
        return 0;
      }
      fclose(fh);
    }
    return ( sinfo1.st_dev == sinfo2.st_dev
             && sinfo1.st_ino == sinfo2.st_ino ) ? 1 : 0;
#endif
  }

  mcu8str mctools_path_join( const mcu8str* p1, const mcu8str* p2 )
  {
    //In general try to mimic os.path.path from Python, although we do not
    //require a drive letter on windows to be considered an absolute path, so we
    //might have slight changes in some weird cases.

#ifdef MC_IS_WINDOWS
    const char native_sep = '\\';
#else
    const char native_sep = '/';
#endif

    if ( p2->size == 0 ) {
      //Simply copy p1, but ensure it ends in a separator:
      char lastchar = p1->c_str[p1->size-1];
      if ( lastchar == '\\' || lastchar == '/' )
        return mcu8str_copy( p1 );
      mcu8str res = mcu8str_create( p1->size + 1 );
      mcu8str_append( &res, p1 );
      res.c_str[res.size] = native_sep;
      ++( res.size );
      res.c_str[res.size] = '\0';
      return res;
    }

    char drive_letter1 = p1->size >= 2 ? mctools_drive_letter(p1) : 0;
    char drive_letter2 = p2->size >= 2 ? mctools_drive_letter(p2) : 0;
    if ( p1->size == 0 || mctools_path_is_absolute(p2) ) {
      //p1 is empty or p2 is an absolute path, so simply discard p1 (like
      //Pythons os.path.join does in this case). However, if there is a drive
      //letter in p1 and not in p2, we use that drive letter:
      if ( drive_letter1 && !drive_letter2 ) {
        mcu8str res = mcu8str_create( p2->size + 2 );
        res.c_str[0] = drive_letter1;
        res.c_str[1] = ':';
        res.c_str[2] = '\0';
        res.size = 2;
        mcu8str_append( &res, p2 );
        return res;
      } else {
        return mcu8str_copy( p2 );
      }
    }

    if ( drive_letter2 && drive_letter1 != drive_letter2 ) {
      //p2 has a drive letter and p1 has a different drive letter => simply
      //return p2.
      return mcu8str_copy( p2 );
    }

    //Neither p1 or p2 are empty here, so we can simply join them. However, if
    //there is a common drive letter or repeated slashes at the end of p1, we do
    //not repeat them needlessley.

    //Selected ranges:
    STDNS size_t iB1 = 0;
    STDNS size_t iE1 = p1->size;
    STDNS size_t iB2 = 0;
    STDNS size_t iE2 = p2->size;

    if ( drive_letter1 && drive_letter1 == drive_letter2 )
      iB2 += 2;//skip common drive letter from p2

    //Strip repeated trailing slashes down to one trailing slash:
    STDNS size_t minlen1 = ( drive_letter1 ? 3 : 1 );
    while ( iE1-iB1 > minlen1
            && ( p1->c_str[iE1-1]=='/' || p1->c_str[iE1-1]=='\\' )
            && ( p1->c_str[iE1-2]=='/' || p1->c_str[iE1-2]=='\\' ) )
      --iE1;

    STDNS size_t needs_slash = ( ( p1->c_str[iE1-1]=='/'
                                   || p1->c_str[iE1-1]=='\\' ) ? 0 : 1 );

    STDNS size_t s1 = iE1-iB1;
    STDNS size_t s2 = iE2-iB2;
    STDNS size_t newsize = s1 + s2 + needs_slash;
    mcu8str res =  mcu8str_create( newsize  );
    STDNS memcpy( res.c_str, p1->c_str, s1 );
    STDNS size_t used = s1;
    if ( needs_slash ) {
      res.c_str[used] = native_sep;
      ++used;
    }
    STDNS memcpy( res.c_str + used, p2->c_str, s2+1 );//+1 to include null char
    res.size = newsize;
    return res;
  }

  mcu8str mctools_fileextension( const mcu8str* path )
  {
    const char * it = mctools_fileextension_view( path );
    return mcu8str_create_from_cstr( it );
  }

  mcu8str mctools_basename( const mcu8str* path )
  {
    const char * itB = path->c_str;
    const char * itE = path->c_str + path->size;
    if ( mctools_drive_letter( path ) )
      itB += 2;//always skip drive letter for basename

    if ( itB == itE )
      return mcu8str_create_empty();

    //Find last path sep:
    const char * it = itE - 1;
    while ( it > itB && *it!='/' && *it!='\\' )
      --it;
    if ( *it=='/' || *it=='\\' )
      ++it;
    STDNS size_t bnsize = (STDNS size_t)(itE - it);
    if ( bnsize == 0 || ( bnsize==1 && *it=='.' ) )
      return mcu8str_create_empty();
    mcu8str res = mcu8str_create( bnsize );
    STDNS memcpy( res.c_str, it, bnsize + 1 );
    return res;
  }

  const char * mctools_basename_view( const mcu8str* path )
  {
    const char * itB = path->c_str;
    const char * itE = path->c_str + path->size;
    if ( mctools_drive_letter( path ) )
      itB += 2;//always skip drive letter for basename
    if ( itB == itE )
      return itE;//empty
    //Find last path sep:
    const char * it = itE - 1;
    while ( it > itB && *it!='/' && *it!='\\' )
      --it;
    if ( *it=='/' || *it=='\\' )
      ++it;
    if ( it+1 == itE && *it == '.' )
      return itE;//empty
    return it;
  }

  const char* mctools_fileextension_view( const mcu8str* path )
  {
    const char * it = mctools_basename_view( path );
    const char * itLastDot = (const char *)0;
    while ( true ) {
      if ( *it == 0 )
        break;//end of string reached
      if ( *it == '.' )
        itLastDot = it;
      ++it;
    }
    if ( itLastDot )
      return itLastDot + 1;
    return it;//empty string
  }

  mcu8str mctools_dirname( const mcu8str* path )
  {
    char drive_letter = mctools_drive_letter( path );
    const char * itB = path->c_str;
    const char * itE = path->c_str + path->size;
    if ( drive_letter )
      itB += 2;//ignore drive letter while analysing

    if ( itB == itE ) {
      //empty:
      if ( drive_letter ) {
        auto res = mcu8str_create_from_cstr("::");
        res.c_str[0] = drive_letter;
        return res;
      } else {
        return mcu8str_create_empty();
      }
    }

    //Find last path sep, skipping consecutive separatorss:
    const char * it0 = itB;
    //const char * itB = ( drive_letter ? it0 + 2 : it0 );
    const char * it = (itE-1);//start at last

    while ( it > itB && *it!='/' && *it!='\\' )
      --it;
    //it is now at the last path sep:
    while ( it > itB && ( *(it-1)=='/' || *(it-1)=='\\' ) )
      --it;

    //If this turned a path starting with a path sep into an empty path, add
    //back a single slash (so that dirname("/")="/" and not ""):
    if ( it == itB && ( *it == '/' || *it == '\\' ) )
      ++it;

    //Apart from a potential drive_letter, the string at [it0,it) is now the
    //result:
    STDNS size_t ressize = (STDNS size_t)(it-it0);

    if ( ressize == 1 && *it0 == '.' ) {
      //special case, "." or "D:."
      if ( drive_letter ) {
        auto res = mcu8str_create_from_cstr("::");
        res.c_str[0] = drive_letter;
        return res;
      } else {
        return mcu8str_create_from_cstr(".");
      }
    }

    if ( ressize == 0 ) {
      if ( path->c_str[0] == '.' && !drive_letter )
        return mcu8str_create_from_cstr(".");
      if ( drive_letter ) {
        auto res = mcu8str_create_from_cstr("::");
        res.c_str[0] = drive_letter;
        return res;
      } else {
        return mcu8str_create_empty();
      }
    }
    if ( drive_letter ) {
      mcu8str res = mcu8str_create( ressize + 2 );
      res.c_str[0] = drive_letter;
      res.c_str[1] = ':';
      STDNS memcpy( res.c_str + 2, it0, ressize );
      res.c_str[ressize+2] = '\0';
      res.size = ressize+2;
      return res;
    } else {
      mcu8str res = mcu8str_create( ressize );
      STDNS memcpy( res.c_str, it0, ressize );
      res.c_str[ressize] = '\0';
      res.size = ressize;
      return res;
    }
  }

#ifdef MCFILEUTILS_CPPNAMESPACE
}
#endif

//Fixme: We should support (and test) windows paths starting with "\\?\"...

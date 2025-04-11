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

// Cross platform string and file system utilities for C.

#if ( defined (_WIN32) || defined (WIN32) )
#  define MC_IS_WINDOWS
#endif
#ifndef MC_IS_WINDOWS
#  ifndef __FreeBSD__
     //These cause problems on freebsd:
#    ifndef _POSIX_C_SOURCE
#      define _POSIX_C_SOURCE 200809L
#    endif
#    ifndef _XOPEN_SOURCE
#      define _XOPEN_SOURCE 500
#    endif
#  endif
#  include <unistd.h>
#  include <limits.h>
#  include <errno.h>
#  include <sys/stat.h>
#endif

//Not pretty, but it works:
#if __cplusplus
#  include "NCrystal/internal/utils/NCCFileUtils.hh"
#else
#  include "NCCFileUtils.h"
#endif

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
#  include <new>//for std::bad_alloc
#  include <cstdint>
namespace MCFILEUTILS_CPPNAMESPACE {
#else
#  include <stdlib.h>
#  include <string.h>
#  include <assert.h>
#  include <stdint.h>
#endif

#define MCTOOLS_STATIC_ASSERT0(COND,MSG) { typedef char mctools_##MSG[(COND)?1:-1]; mctools_##MSG dummy; (void)dummy; }
#define MCTOOLS_STATIC_ASSERT3(expr,x) MCTOOLS_STATIC_ASSERT0(expr,fail_at_line_##x)
#define MCTOOLS_STATIC_ASSERT2(expr,x) MCTOOLS_STATIC_ASSERT3(expr,x)
#define MCTOOLS_STATIC_ASSERT(expr)    MCTOOLS_STATIC_ASSERT2(expr,__LINE__)

#ifdef STDNS
#  undef STDNS
#endif
#ifdef MCFILEUTILS_CPPNAMESPACE
#  define STDNS std::
#else
#  define STDNS
#endif


//MSVC C4702 ("unreachable code") triggers in the function below, but there are
//no branches in the function, and it IS being called?!?:
#ifdef _MSC_VER
#  pragma warning( push )
#  pragma warning( disable : 4702 )
#endif
  mcu8str mcu8str_create_empty(void)
  {
    static char dummy[4] = { 0, 0, 0, 0 };
    mcu8str s;
    s.c_str = &dummy[0];//NB: C4702 triggers on this line!
    s.size = 0;
    s.buflen = 0;
    s.owns_memory = 0;
    return s;
  }
#ifdef _MSC_VER
#  pragma warning( pop )
#endif

  void mcu8str_dealloc( mcu8str* str )
  {
    char * buf_to_free = ( str->owns_memory
                           ? str->c_str
                           : NULL );
    static char dummy[4] = { 0, 0, 0, 0 };
    str->c_str = &dummy[0];
    str->size = 0;
    str->buflen = 0;
    str->owns_memory = 0;
    if ( buf_to_free )
      STDNS free( buf_to_free );
  }

#ifdef MCFILEUTILS_CPPNAMESPACE
  namespace {
#endif
    void mctools_impl_error( const char * msg )
    {
#ifdef __cplusplus
      throw std::runtime_error(msg);
#else
      fprintf(stderr, "%s\n",msg);
      exit(EXIT_FAILURE);
#endif
    }
#ifdef MCFILEUTILS_CPPNAMESPACE
  }
#endif

#define MCTOOLS_STROKFORUINT( x ) ( (x) < (UINT_MAX-1) )
  mcu8str_size_t mctools_strlen( const char * c_str, mcu8str_size_t nmax )
  {
    //safe strlen which guarantees result can fit in an unsigned integer (or
    //less, when 0<nmax<UINT_MAX is provided).
    if ( nmax == 0 )
      nmax = ( UINT_MAX > PTRDIFF_MAX ? PTRDIFF_MAX : UINT_MAX );

    if ( nmax > PTRDIFF_MAX )
      mctools_impl_error("str length out of range");

    const char * nullchr = (const char *) STDNS memchr( c_str, '\0', nmax );
    mcu8str_size_t o_size;
    if ( nullchr )
      o_size = nullchr - c_str;
    else
      o_size = nmax;
    if ( !MCTOOLS_STROKFORUINT(o_size) )
      mctools_impl_error("str length out of range");
    return o_size;
  }

  mcu8str mcu8str_create( mcu8str_size_t prealloc_size )
  {
    MCTOOLS_STATIC_ASSERT( sizeof(int)==sizeof(unsigned) );
    MCTOOLS_STATIC_ASSERT( sizeof(size_t)>=sizeof(unsigned) );
    MCTOOLS_STATIC_ASSERT( sizeof(unsigned)>=4 );
    MCTOOLS_STATIC_ASSERT( sizeof(size_t)>=4 );
    MCTOOLS_STATIC_ASSERT( sizeof(void*)>=4 );

    if ( prealloc_size == 0 )
      return mcu8str_create_empty();
    mcu8str s;
    //Limit prealloc_size, since we use unsigned type to store the size in
    //buflen below:
    s.c_str = (char*) ( MCTOOLS_STROKFORUINT(prealloc_size)
                        ? STDNS malloc( prealloc_size+1 ) : NULL);
    if ( !s.c_str ) {
#ifdef MCFILEUTILS_CPPNAMESPACE
      throw std::bad_alloc();
#else
      fprintf(stderr, "ERROR: Memory allocation failed in mcu8str_create\n");
      exit(EXIT_FAILURE);
#endif
      //gives (correctly but annoyingly) unreachable warning with MSVC
      //return mcu8str_create_empty();
    }
    s.c_str[0] = '\0';
    s.size = 0;
    s.buflen = (unsigned)(prealloc_size + 1);//cast ok, checked above
    s.owns_memory = 1;
    return s;
  }

  void mcu8str_append( mcu8str* str, const mcu8str * otherstr )
  {
    if ( otherstr->size == 0 )
      return;
    mcu8str_size_t newsize = str->size + otherstr->size;
    if ( newsize + 1 > str->buflen )
      mcu8str_reserve( str, newsize );
    //memcpy ok since we handle the null char afterwards, there will be no
    //overlapping even if str==otherstr:
    STDNS memcpy( str->c_str + str->size,
                  otherstr->c_str,
                  otherstr->size );
    str->c_str[newsize] = '\0';
    str->size = (unsigned)newsize;//cast ok, mcu8str_reserve checked.
  }

  void mcu8str_append_cstr( mcu8str* str, const char * c_str )
  {
    const mcu8str_size_t o_size = mctools_strlen( c_str, 0 );
    if ( o_size == 0 )
      return;
    const mcu8str_size_t newsize = str->size + o_size;
    if ( newsize + 1 > str->buflen )
      mcu8str_reserve( str, newsize );
    //memcpy ok since we handle the null char afterwards, there will be no
    //overlapping even if c_str = str.c_str:
    STDNS memcpy( str->c_str + str->size, c_str, o_size );
    str->c_str[newsize] = '\0';
    str->size = (unsigned)newsize;//cast ok, mcu8str_reserve checked.
  }

  mcu8str mcu8str_create_from_cstr( const char * c_str )
  {
    if ( ! *c_str )
      return mcu8str_create_empty();
    const mcu8str_size_t o_size = mctools_strlen( c_str, 0 );
    assert( o_size > 0 );
    mcu8str str = mcu8str_create( o_size );
    //str.c_str is newly allocated, so there can be no overlap:
    STDNS memcpy( str.c_str, c_str, o_size + 1 );
    str.size = (unsigned)o_size;//cast ok, mctools_strlen checked.
    return str;
  }

#ifdef MCFILEUTILS_CPPNAMESPACE
  const
#endif
 mcu8str mcu8str_view_cstr( const char * c_str )
  {
    mcu8str s;
    s.c_str = (char*)c_str;//cast away constness, but we will add it again in
                           //the return type (in c++).
    s.size = (unsigned)mctools_strlen( c_str, 0 );//cast due to mctools_strlen
    s.buflen = s.size + 1;
    s.owns_memory = 0;
    return s;
  }

#ifdef MCFILEUTILS_CPPNAMESPACE
  const
#endif
 mcu8str mcu8str_view_str( const mcu8str * str )
  {
    mcu8str s;
    s.c_str = (char*)str->c_str;//cast away constness, but we will add it again in
                           //the return type (in c++).
    s.size = str->size;
    s.buflen = s.size + 1;//limit as much as possible.
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
    mcu8str_size_t n = mctools_strlen( str->c_str, str->buflen );
    if ( n >= str->buflen ) {
      mctools_impl_error("mcu8str_update_size logic error");
      //return;
    }
    str->size = (unsigned)n;
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

  void mcu8str_reserve( mcu8str* str, mcu8str_size_t nsize )
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

  int mcu8str_equal( const mcu8str* s1, const mcu8str* s2 )
  {
    return ( s1->size == s2->size
             && STDNS memcmp( s1->c_str, s2->c_str, s1->size )==0 );
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
  mcu8str mcu8str_create_from_staticbuffer( char * buf, mcu8str_size_t buflen )
  {
    if ( buflen==0 || !MCTOOLS_STROKFORUINT(buflen-1) )
      mctools_impl_error("static buffer length out of range");
    assert(buflen > 0);
    mcu8str s;
    s.owns_memory = 0;
    s.c_str = buf;
    s.c_str[0] = '\0';
    s.size = 0;
    s.buflen = (unsigned)buflen;
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
#else
#  ifdef __APPLE__
#    define MC_IS_APPLE
#    include <mach-o/dyld.h>
#  endif
#endif

#ifdef MCFILEUTILS_CPPNAMESPACE
#  include <cstdlib>
#  include <cstring>
#  include <cassert>
#else
#  include <stdlib.h>
#  include <string.h>
#  include <assert.h>
#endif

#ifdef __cplusplus
#  ifdef NULL
#    undef NULL
#  endif
#  define NULL nullptr
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
  DWORD GetLongPathNameW( const wchar_t*, wchar_t*, DWORD);
  DWORD GetLastError(void);
  FILE * _wfopen( const wchar_t*,const wchar_t*);
#  define errno_t int
  errno_t _wfopen_s( FILE **, const wchar_t*,const wchar_t*);
  typedef struct { int dummy; } FILETIME;
  typedef struct {
    DWORD    dwFileAttributes;
    FILETIME ftCreationTime;
    FILETIME ftLastAccessTime;
    FILETIME ftLastWriteTime;
    DWORD    dwVolumeSerialNumber;
    DWORD    nFileSizeHigh;
    DWORD    nFileSizeLow;
    DWORD    nNumberOfLinks;
    DWORD    nFileIndexHigh;
    DWORD    nFileIndexLow;
  } BY_HANDLE_FILE_INFORMATION;
  //int PathIsRelativeW( const wchar_t* );
  int CompareFileTime( const FILETIME *,const FILETIME *);
  int _fileno( FILE * );
  typedef void* HANDLE;
  int GetFileInformationByHandle( HANDLE, BY_HANDLE_FILE_INFORMATION* );
  DWORD GetFinalPathNameByHandleW( HANDLE, wchar_t*, DWORD, DWORD);
  void CloseHandle( HANDLE );
  HANDLE CreateFileW( const wchar_t*,DWORD,DWORD,void*,DWORD,DWORD,HANDLE );
  typedef void* LPSECURITY_ATTRIBUTES;
  //intptr_t _get_osfhandle(int);
#  define FILE_SHARE_READ 0x00000001
#  define FILE_SHARE_WRITE 0x00000002
#  define FILE_SHARE_DELETE 0x00000004
#  define OPEN_EXISTING 12345
#  define FILE_ATTRIBUTE_NORMAL 12345
#  define INVALID_HANDLE_VALUE 0
#  define FILE_READ_ATTRIBUTES 12345
#  define FILE_FLAG_BACKUP_SEMANTICS 12345
#  define GENERIC_READ 12345
#  define FILE_ATTRIBUTE_REPARSE_POINT 12345
#  define ERROR_CANT_ACCESS_FILE 12345
#  define FILE_FLAG_OPEN_REPARSE_POINT 12345
#endif

#if 0
  //Fake being on macos to check basic compilation
#  define MC_IS_APPLE
  int _NSGetExecutablePath(char*, uint32_t*);
#endif

#ifdef MCFILEUTILS_CPPNAMESPACE
namespace MCFILEUTILS_CPPNAMESPACE {
namespace {
#endif

#ifdef MC_IS_WINDOWS
#  ifndef DWORD_MAX
#    define DWORD_MAX 0xffffffffUL
#  endif

  int mctools_impl_str_size2int( mcu8str_size_t v )
  {
#  ifdef INT_MAX
    const int lim = INT_MAX - 1;
#  else
    const int lim = 2147483646;
#  endif
    if ( !( v < lim && v + 1 < lim) )
      mctools_impl_error("str length out of range");
    return (int)(v);
  }
#endif

  int mctools_impl_has_winnamespace( const mcu8str * path )
  {
    return ( path->size >= 4 && path->c_str[2] == '?'
             && ( path->c_str[0] == '/' || path->c_str[0] == '\\' )
             && ( path->c_str[1] == '/' || path->c_str[1] == '\\' )
             && ( path->c_str[3] == '/' || path->c_str[3] == '\\' ) );
  }

#ifdef MCFILEUTILS_CPPNAMESPACE
  //We want to return a const view, but can't in c:
  const
#endif
  mcu8str mctools_impl_view_no_winnamespace( const mcu8str * path )
  {
    //Windows paths starting with "\\?\" are a bit special. For our purposes we
    //simply skip those initial characters. See more on:
    //https://learn.microsoft.com/en-us/windows/win32/fileio/naming-a-file
    mcu8str res = mcu8str_view_str( path );
    if ( mctools_impl_has_winnamespace(path) ) {
      //Begins with "<sep><sep>?<sep>", modify the view:
      res.c_str += 4;
        res.size -= 4;
    }
    return res;
  }

#ifdef MCFILEUTILS_CPPNAMESPACE
}
#endif



#ifdef MC_IS_WINDOWS
#  ifdef __cplusplus
  namespace {
#  endif
    //Need wide-char strings for unicode API. We try to be less fancy with these
    //than the mcu8str's, since they are anyway only used on Windows and in this
    //file.
    typedef struct {
      wchar_t * c_str;
      mcu8str_size_t size;
      mcu8str_size_t buflen;
    } mcwinstr;

    void mc_winstr_swap( mcwinstr* p1, mcwinstr* p2 )
    {
      mcwinstr tmp = *p1;
      p1->c_str = p2->c_str;
      p1->size = p2->size;
      p1->buflen = p2->buflen;
      p2->c_str = tmp.c_str;
      p2->size = tmp.size;
      p2->buflen = tmp.buflen;
    }

    mcwinstr mc_winstr_create( mcu8str_size_t );
    mcwinstr mc_winstr_create_empty(void);
    void mc_winstr_dealloc( mcwinstr* );

    mcwinstr mc_u8str_to_winstr( const mcu8str* src )
    {
      const int in_size = mctools_impl_str_size2int( src->size );

      //First check for empty string:
      if ( !(src->c_str[0]) || in_size == 0 )
        return mc_winstr_create_empty();

      const char * in_data = src->c_str;
      int out_size = MultiByteToWideChar( CP_UTF8,
                                          0,//Must be 0 for utf8
                                          in_data, in_size,
                                          NULL, //dest buffer (NULL for dry-run)
                                          0//dest buffer size (0 means
                                          //"dry-run" + return needed size).
                                          );
      const char * errmsg = "Failed to convert UTF-8 string to UTF-16";
      if (!out_size) {
        mctools_impl_error(errmsg);
        //return mc_winstr_create_empty();
      }

      mcwinstr res = mc_winstr_create( (mcu8str_size_t)( out_size ) );

      wchar_t * out_data = res.c_str;
      //Same, but with out_data/out_size provided:
      int out_size2 = MultiByteToWideChar( CP_UTF8, 0,
                                           in_data, in_size,
                                           out_data, out_size );
      if ( out_size != out_size2 || (mcu8str_size_t)out_size >= res.buflen ) {
        mc_winstr_dealloc( &res );
        mctools_impl_error(errmsg);
        //return mc_winstr_create_empty();
      }
      res.c_str[out_size] = 0;
      res.size = out_size;
      return res;
    }

    mcwinstr mc_path2wpath( const mcu8str* src )
    {
      //Like mc_u8str_to_winstr but runs input through mctools_pathseps_platform
      //for consistency. Usually without any unnecessary malloc/free calls.
      char buf[4096];
      mcu8str srccopy = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
      mcu8str_assign( &srccopy, src );
      mctools_pathseps_platform( &srccopy );
      mcwinstr res = mc_u8str_to_winstr( &srccopy );
      mcu8str_dealloc(&srccopy);
      return res;
    }

    mcwinstr mc_winstr_copy( const mcwinstr* path )
    {
      if ( path->size == 0 )
        return mc_winstr_create_empty();
      mcwinstr str = mc_winstr_create( path->size );
      STDNS memcpy( str.c_str, path->c_str, sizeof(wchar_t)*(path->size + 1) );
      str.size = path->size;
      return str;
    }

    mcwinstr mc_impl_expand_wpath_to_longpathname( const mcwinstr* sp )
    {
      if ( sp->size == 0 )
        return mc_winstr_create_empty();
      mcwinstr out = mc_winstr_create( 4096 );
      assert( out.buflen < DWORD_MAX );
      DWORD len = GetLongPathNameW( sp->c_str, out.c_str, (DWORD)out.buflen );
      if ( len == 0 ) {
        //failure (e.g. file does not exist). We return the original string
        //unchanged in this case:
        mc_winstr_dealloc( &out );
        return mc_winstr_copy( sp );
      }
      if ( len < 4096 ) {
        out.c_str[len] = 0;
        out.size = (mcu8str_size_t)len;
        return out;
      }
      mc_winstr_dealloc( &out );
      if ( len >= DWORD_MAX )
        mctools_impl_error("str length out of range for DWORD");
      out = mc_winstr_create( (DWORD)len );
      len = GetLongPathNameW( sp->c_str, out.c_str, (DWORD)out.buflen );
      if ( (mcu8str_size_t) len < out.buflen ) {
        out.c_str[len] = 0;
        out.size = (mcu8str_size_t)len;
        return out;
      } else {
        //should not happen
        mc_winstr_dealloc( &out );
        return mc_winstr_copy( sp );
      }
    }

    mcwinstr mc_winstr_expand_to_fullpath( const mcwinstr* wpath )
    {
      //nb: in case of problems might throw, prints to stderr or return empty
      //string.

      DWORD len_with_null_term = GetFullPathNameW( wpath->c_str,
                                                   (DWORD)0, NULL, NULL );
      if ( len_with_null_term <= 1 ) {
        mctools_impl_error( "Full path expansion with GetFullPathNameW failed" );
        //return mc_winstr_create_empty();
      }
      mcwinstr woutput = mc_winstr_create((mcu8str_size_t)(len_with_null_term-1));
      DWORD len = GetFullPathNameW( wpath->c_str,
                                    (DWORD)woutput.buflen, woutput.c_str,
                                    NULL );
      if ( len+1 != len_with_null_term ) {
        mc_winstr_dealloc( &woutput );
        mctools_impl_error( "Full path expansion with GetFullPathNameW failed" );
        //return mc_winstr_create_empty();
      }
      //OK, update output size:
      woutput.c_str[len] = 0;
      woutput.size = (mcu8str_size_t)len;
      //For some reason the above might have returned a path with short 8.3
      //filenames, so we expand once more (only works if file exists):
      if ( woutput.size > 0 ) {
        mcwinstr woutput2 = mc_impl_expand_wpath_to_longpathname( &woutput );
        if ( woutput2.size > 0 ) {
          //ok, use expanded instead:
          mc_winstr_swap( &woutput, &woutput2 );
        }
        mc_winstr_dealloc(&woutput2);
      }
      return woutput;
    }

    mcwinstr mc_path2longwpath( const mcu8str* src )
    {
      mcwinstr wp = mc_path2wpath( src );
      if ( wp.size == 0 )
        return wp;
      mcwinstr lwp = mc_impl_expand_wpath_to_longpathname( &wp );
      if ( lwp.size == 0 ) {
        //something went wrong, return wp.
        mc_winstr_dealloc(&lwp);
        return wp;
      }
      mc_winstr_dealloc(&wp);
      return lwp;
    }

    mcu8str mc_winstr_to_u8str( const mcwinstr* src )
    {
      const int in_size = mctools_impl_str_size2int( src->size );

      //First check for empty string (also for safeguard in case of deallocated
      //string):
      if ( !src->c_str[0] || in_size == 0 )
        return mcu8str_create_empty();

      const wchar_t * in_data = src->c_str;

      int out_size = WideCharToMultiByte( CP_UTF8,
                                          0,//Must be 0 for utf8
                                          in_data, in_size,
                                          NULL, //dest buffer (NULL for dry-run)
                                          0,//dest buffer size (0 means "dry-run"
                                            //returning, needed size)
                                          NULL,//Must be null for utf8
                                          NULL//Must be null for utf8
                                          );
      //NB: out_size does not include 0 termination char!
      const char * errmsg = "Failed to convert UTF-16 string to UTF-8";

      if (!out_size) {
        mctools_impl_error(errmsg);
        //return mcu8str_create_empty();
      }
      mcu8str res = mcu8str_create( out_size );
      //res.resize( out_size );
      char * out_data = res.c_str;
      //Same, but with out_data/out_size provided:
      int out_size2 = WideCharToMultiByte( CP_UTF8, 0,
                                           in_data, in_size,
                                           out_data, out_size,
                                           NULL, NULL);
      if ( out_size2 != out_size || (mcu8str_size_t)out_size >= res.buflen ) {
        mcu8str_dealloc(&res);
        mctools_impl_error(errmsg);
        //return mcu8str_create_empty();
      }
      res.c_str[out_size] = 0;
      res.size = out_size;
      return res;
    }

    mcu8str mc_winstr_expand_to_fullpath_u8str( const mcwinstr* wpath )
    {
      //nb: in case of problems might throw, prints to stderr or return empty
      //string.
      mcwinstr fp = mc_winstr_expand_to_fullpath( wpath );
      mcu8str res = mc_winstr_to_u8str( &fp );
      mc_winstr_dealloc(&fp);
      if ( mctools_impl_has_winnamespace( &res ) ) {
        mcu8str res2 = mctools_impl_view_no_winnamespace( &res );
        mcu8str_ensure_dynamic_buffer(&res2);
        mcu8str_swap( &res, &res2 );
        mcu8str_dealloc(&res);
      }
      return res;
    }

    mcwinstr mc_winstr_create_empty(void)
    {
      //empty (keep in allocated buffer for simplicity).
      return mc_winstr_create(0);
    }

    mcwinstr mc_winstr_create( mcu8str_size_t size )
    {
      MCTOOLS_STATIC_ASSERT( sizeof(DWORD)==sizeof(uint32_t) );
      MCTOOLS_STATIC_ASSERT( sizeof(unsigned)>=sizeof(uint32_t) );
      MCTOOLS_STATIC_ASSERT( sizeof(size_t)>=sizeof(uint32_t) );
      mcwinstr str;
      str.c_str = (wchar_t*)( MCTOOLS_STROKFORUINT(size)
                              ?  STDNS malloc( sizeof(wchar_t)*(size + 1) )
                              : NULL );
      if ( !str.c_str ) {
#ifdef MCFILEUTILS_CPPNAMESPACE
        throw std::bad_alloc();
#else
        fprintf(stderr, "ERROR: Memory allocation failed in mc_winstr_create\n");
        exit(EXIT_FAILURE);
#endif
        //gives (correctly but annoyingly) unreachable warning with MSVC
        //return mc_winstr_create_empty();
      }
      str.c_str[0] = 0;
      str.size = 0;
      str.buflen = size + 1;
      return str;
    }

    void mc_winstr_dealloc( mcwinstr* str )
    {
      if ( !str->c_str )
        return;
      void * bufptr = (void*) str->c_str;
      str->size = 0;
      str->buflen = 0;
      str->c_str = NULL;
      STDNS free( bufptr );
    }

#  ifdef __cplusplus
  }//end of anonymous namespace
#  endif
#endif

  MCTOOLS_FILE_t * mctools_fopen( const mcu8str* rawpath,
                                  const char * mode )
  {
    mcu8str path = mctools_impl_view_no_winnamespace( rawpath );

    //Utf8-encoded paths are used directly on platforms except Windows, where
    //only ASCII encoded paths can be used directly. We also run paths through
    //mctools_pathseps_platform.
#ifdef MC_IS_WINDOWS
    //Ok, non-ASCII path on windows. Reencode as wide strings and use _wfopen:
    mcwinstr wpath = mc_path2wpath( &path );//this also maps to native path sep
    char buf[32];
    mcu8str mcu8str_mode = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
    mcu8str_append_cstr(&mcu8str_mode,mode);
    mcwinstr wmode = mc_u8str_to_winstr( &mcu8str_mode );
    //    MCTOOLS_FILE_t * result = _wfopen( wpath.c_str, wmode.c_str );//nb: no STDNS here on purpose
    MCTOOLS_FILE_t * result = NULL;
    errno_t errno_result = _wfopen_s( &result, wpath.c_str, wmode.c_str );
    mcu8str_dealloc(&mcu8str_mode);
    mc_winstr_dealloc( &wpath );
    mc_winstr_dealloc( &wmode );
    if ( errno_result != 0 ) {
      if ( result )
        fclose(result);//a bit weird, but most likely we won't get here.
      return NULL;
    }
    return result;
#else
    if ( !mcu8str_contains( &path, '\\' ) ) {
      //Usual case on unix -> no extra work needed
      return STDNS fopen( path.c_str, mode );
    } else {
      //Copy path so we can run it through mctools_pathseps_platform:
      char buf[4096];
      mcu8str pathcopy = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
      mcu8str_assign( &pathcopy, &path );
      mctools_pathseps_platform( &pathcopy );
      MCTOOLS_FILE_t * fh = STDNS fopen( pathcopy.c_str, mode );
      mcu8str_dealloc(&pathcopy);
      return fh;
    }
#endif
  }

  int mctools_is_dir( const mcu8str* rawpath )
  {
    mcu8str path = mctools_impl_view_no_winnamespace( rawpath );
#ifdef MC_IS_WINDOWS
    mcwinstr wpath = mc_path2wpath( &path );
    DWORD fa = GetFileAttributesW( wpath.c_str );
    mc_winstr_dealloc( &wpath );
    if ( fa == INVALID_FILE_ATTRIBUTES )
      return 0;
    return ( fa & FILE_ATTRIBUTE_DIRECTORY) ? 1 : 0;
#else
    char buf[4096];
    mcu8str native = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
    mcu8str_assign( &native, &path );
    mctools_pathseps_platform( &native );
    struct stat sinfo;
    int res = ( stat( native.c_str, &sinfo) == 0
                && S_ISDIR(sinfo.st_mode) ) ? 1 : 0;
    mcu8str_dealloc(&native);
    return res;
#endif
  }

  int mctools_exists( const mcu8str* rawpath )
  {
    mcu8str path = mctools_impl_view_no_winnamespace( rawpath );
#ifdef MC_IS_WINDOWS
    mcwinstr wpath = mc_path2wpath( &path );
    DWORD fa = GetFileAttributesW( wpath.c_str );
    mc_winstr_dealloc( &wpath );
    return fa == INVALID_FILE_ATTRIBUTES ? 0 : 1;
#else
    char buf[4096];
    mcu8str native = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
    mcu8str_assign( &native, &path );
    mctools_pathseps_platform( &native );
    struct stat sinfo;
    int res = ( stat( native.c_str, &sinfo) == 0 ) ? 1 : 0;
    mcu8str_dealloc(&native);
    return res;
#endif
  }

  int mctools_is_file( const mcu8str* rawpath )
  {
    mcu8str path = mctools_impl_view_no_winnamespace( rawpath );
#ifdef MC_IS_WINDOWS
    mcwinstr wpath = mc_path2wpath( &path );
    DWORD fa = GetFileAttributesW( wpath.c_str );
    mc_winstr_dealloc( &wpath );
    if ( fa == INVALID_FILE_ATTRIBUTES )
      return 0;
    if ( fa & FILE_ATTRIBUTE_DIRECTORY)
      return 0;
    return 1;
#else
    char buf[4096];
    mcu8str native = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
    mcu8str_assign( &native, &path );
    mctools_pathseps_platform( &native );
    struct stat sinfo;
    int res = ( stat( native.c_str, &sinfo) == 0
                && !S_ISDIR(sinfo.st_mode) ) ? 1 : 0;
    mcu8str_dealloc(&native);
    return res;
#endif
  }

  mcu8str mctools_get_current_working_dir(void)
  {
#ifdef MC_IS_WINDOWS
    mcwinstr wpath = mc_winstr_create( 4096 );
    DWORD nsize = GetCurrentDirectoryW((DWORD)wpath.buflen, wpath.c_str);
    if ( (mcu8str_size_t)(nsize + 1) > wpath.buflen ) {
      //Use larger buffer and try again:
      mc_winstr_dealloc( &wpath );
      wpath = mc_winstr_create( nsize );
      nsize = GetCurrentDirectoryW((DWORD)wpath.buflen,wpath.c_str);
    }
    if ( nsize == 0 || (mcu8str_size_t )( nsize + 1 ) > wpath.buflen ) {
      mc_winstr_dealloc( &wpath );
      mctools_impl_error("Failed to get current working directory");
      //return mcu8str_create_empty();
    }
    wpath.c_str[nsize] = 0;
    wpath.size = nsize;
    mcwinstr woutput = mc_impl_expand_wpath_to_longpathname( &wpath );
    mc_winstr_dealloc( &wpath );
    mcu8str res = mc_winstr_to_u8str( &woutput );
    mc_winstr_dealloc( &woutput );
    mctools_pathseps_platform( &res );
    return res;
#else
    char buf[4096];//almost always enough in first go
    mcu8str res = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
    while ( 1 ) {
      if ( getcwd( res.c_str, res.buflen ) ) {
        mcu8str_update_size( &res );
        mcu8str_ensure_dynamic_buffer(&res);//don't return local static buffer
        mctools_pathseps_platform( &res );
        return res;
      }
      if ( errno == ERANGE && res.buflen < 2000000 ) {
        //Try again with larger buffer:
        mcu8str_clear( &res );
        mcu8str_reserve( &res, (res.buflen-1) * 2 );
      } else {
        mcu8str_dealloc( &res );
        mctools_impl_error("Failed to get current working directory");
        //return mcu8str_create_empty();
      }
    }
#endif
  }

  mcu8str mctools_determine_exe_self_path( int argc, char** argv )
  {
    //First try some platform specific methods not relying on argv0:
#ifdef MC_IS_WINDOWS
    {
      mcwinstr wpath = mc_winstr_create( (mcu8str_size_t)(MAX_PATH<32768
                                                          ?32768:MAX_PATH) +1 );
      DWORD nsize = GetModuleFileNameW(NULL, wpath.c_str, (DWORD)wpath.buflen );
      if ( nsize == 0 || (mcu8str_size_t)(nsize + 2) >= wpath.buflen ) {//+2 is safety
        //failure:
        mc_winstr_dealloc( &wpath );
      } else {
        //success:
        assert( (mcu8str_size_t)(nsize)+2 <= wpath.buflen );//+2 is safety
        wpath.c_str[nsize] = 0;//add null terminator
        wpath.size = nsize;
        mcu8str res = mc_winstr_to_u8str( &wpath );
        mc_winstr_dealloc( &wpath );
        mctools_pathseps_platform(&res);
        return res;
      }
    }
#endif
#ifdef MC_IS_APPLE
    {
      mcu8str path = mcu8str_create( 4096 );
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
        if ( path.size > 0 ) {
          mctools_pathseps_platform(&path);
          return path;
        }
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
        //NB: ssize_t is from unistd.h, so never exists in std:: namespace.
        ssize_t len = readlink(filename, path.c_str, path.buflen-1 );
        if ( len > 0 && (mcu8str_size_t)(len+1) < path.buflen ) {
          path.c_str[len] = '\0';//readlink does not add terminating null char
          path.size = (mcu8str_size_t)len;
          mcu8str_ensure_dynamic_buffer( &path );
          mctools_pathseps_platform(&path);
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
      mctools_pathseps_platform( &path );
      if ( mctools_is_file( &path ) ) {
        mctools_pathseps_platform(&path);
        return path;
      }
      mcu8str_dealloc( &path );
    }

    //We could do more extensive analysis here, but it might give wrong results
    //if we are unlucky.
    return mcu8str_create_empty();
  }

  mcu8str mctools_real_path( const mcu8str* rawpath )
  {
    //NB: Might return empty in case of errors.
    mcu8str path = mctools_impl_view_no_winnamespace(rawpath);
    if ( path.size == 0 )
      return mcu8str_create_empty();

#ifdef MC_IS_WINDOWS
    //To resolve symlinks we must open a file handle with CreateFileW. If it
    //fails, we revert back to at least returning an absolute path.
    mcwinstr wpath = mc_path2wpath( &path );
    DWORD share = FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE;
    DWORD access = FILE_READ_ATTRIBUTES;
    DWORD flags = FILE_FLAG_BACKUP_SEMANTICS;
    HANDLE fh1 = CreateFileW( wpath.c_str, access, share, NULL,
                              OPEN_EXISTING, flags, NULL );
    if ( fh1 == INVALID_HANDLE_VALUE ) {
      //Failed to open. This can happen for a variety of reasons. As a fallback
      //option, we return the full path (which won't resolve symlinks):
      mcu8str output = mc_winstr_expand_to_fullpath_u8str( &wpath );
      mc_winstr_dealloc( &wpath );
      return output;
    }

    //We have an open handle fh1, use it to find the resolved path:
    mcwinstr resolvedpath = mc_winstr_create(4096);
    DWORD len = GetFinalPathNameByHandleW( fh1, resolvedpath.c_str,
                                           (DWORD)resolvedpath.buflen, 0 );
    if ( (mcu8str_size_t)len >= resolvedpath.buflen ) {
      //Too short buffer, try again:
      mc_winstr_dealloc(&resolvedpath);
      resolvedpath = mc_winstr_create(len);
      len = GetFinalPathNameByHandleW( fh1, resolvedpath.c_str,
                                       (DWORD)resolvedpath.buflen, 0 );
    }
    //Close the open file handle:
    CloseHandle( fh1 );
    if ( len == 0 || (mcu8str_size_t)len >= resolvedpath.buflen ) {
      //Failure, fall-back to unresolved absolute path:
      mc_winstr_dealloc(&resolvedpath);
      mcu8str output = mc_winstr_expand_to_fullpath_u8str( &wpath );
      mc_winstr_dealloc( &wpath );
      return output;
    }
    mc_winstr_dealloc( &wpath );
    resolvedpath.c_str[len] = 0;
    resolvedpath.size = len;
    //All ok, get rid of any 8.3-form path names and return:
    mcwinstr woutput = mc_impl_expand_wpath_to_longpathname( &resolvedpath );
    mc_winstr_dealloc(&resolvedpath);
    mcu8str output = mc_winstr_to_u8str( &woutput );
    mc_winstr_dealloc( &woutput );
    if ( mctools_impl_has_winnamespace( &output ) ) {
      mcu8str tmp = mctools_impl_view_no_winnamespace( &output );
      mcu8str_ensure_dynamic_buffer(&tmp);
      mcu8str_swap( &tmp, &output );
      mcu8str_dealloc(&tmp);
    }
#else
    char buf[4096];
    mcu8str native = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
    mcu8str_assign( &native, &path );
    mctools_pathseps_platform( &native );
    mcu8str output = mcu8str_create( PATH_MAX );
    if ( !realpath( native.c_str, output.c_str ) ) {
      //failure
      mcu8str_dealloc(&output);
      mcu8str_dealloc(&native);
      return mcu8str_create_empty();
    }
    mcu8str_dealloc(&native);
    mcu8str_update_size(&output);
#endif
    mctools_pathseps_platform(&output);
    return output;
  }

  void mctools_pathseps_generic( mcu8str* path )
  {
    mcu8str_replace( path, '\\', '/' );
    //Upper case drive letter for consistency:
    char drive_letter = mctools_drive_letter( path );
    if ( drive_letter ) {
      if ( mctools_impl_has_winnamespace( path ) )
        path->c_str[4] = drive_letter;
      else
        path->c_str[0] = drive_letter;
    }
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
    if ( drive_letter ) {
      if ( mctools_impl_has_winnamespace( path ) )
        path->c_str[4] = drive_letter;
      else
        path->c_str[0] = drive_letter;
    }
  }

  int mctools_path_is_relative( const mcu8str* path )
  {
    return ( mctools_path_is_absolute(path) ? 0 : 1 );
  }

  int mctools_path_is_absolute( const mcu8str* rawpath )
  {
    mcu8str path = mctools_impl_view_no_winnamespace(rawpath);
    if ( path.size == 0 )
      return 0;
    if ( path.size >= 2
         && path.c_str[0] == '~'
         && ( path.c_str[1] == '/' || path.c_str[1] == '\\' ) )
      return 1;//leading "~/"
    const char * path_begin = path.c_str;
    if ( mctools_drive_letter( &path ) )
      path_begin += 2;
    if ( *path_begin == '/' || *path_begin == '\\' )
      return 1;
    return 0;
  }

  char mctools_drive_letter( const mcu8str* rawpath )
  {
    mcu8str path = mctools_impl_view_no_winnamespace(rawpath);
    if ( path.size < 2 )
      return 0;
    if ( path.c_str[1] != ':' )
      return 0;
    if ( path.c_str[0] >= 'A' && path.c_str[0] <= 'Z' )
      return path.c_str[0];
    if ( path.c_str[0] >= 'a' && path.c_str[0] <= 'z' )
      return (char)(path.c_str[0] - 32);//uppercase
    return 0;
  }

  int mctools_is_same_file( const mcu8str* p1raw, const mcu8str* p2raw )
  {
    mcu8str p1 = mctools_impl_view_no_winnamespace(p1raw);
    mcu8str p2 = mctools_impl_view_no_winnamespace(p2raw);

    //First check for trivial string equality:
    if ( mcu8str_equal( &p1, &p2 ) )
      return mctools_is_file( &p1 );

#ifdef MC_IS_WINDOWS
    //Open files and query file information. Note that we try to keep both
    //handles open while querying info, for consistent results.
    mcwinstr wp1 = mc_path2wpath( &p1 );
    mcwinstr wp2 = mc_path2wpath( &p2 );

    DWORD share = FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE;
    DWORD access = FILE_READ_ATTRIBUTES;
    DWORD flags = FILE_FLAG_BACKUP_SEMANTICS;
    HANDLE fh1 = CreateFileW( wp1.c_str,
                              access,
                              share,
                              NULL,
                              OPEN_EXISTING,
                              flags,
                              NULL );
    mc_winstr_dealloc( &wp1 );
    if ( fh1 == INVALID_HANDLE_VALUE ) {
      mc_winstr_dealloc( &wp2 );
      return 0;//failed
    }
    HANDLE fh2 = CreateFileW( wp2.c_str,
                              access,
                              share,
                              NULL,
                              OPEN_EXISTING,
                              flags,
                              NULL );
    mc_winstr_dealloc( &wp2 );
    if ( fh2 == INVALID_HANDLE_VALUE ) {
      CloseHandle( fh1 );
      return 0;//failed
    }

    BY_HANDLE_FILE_INFORMATION info1;
    BY_HANDLE_FILE_INFORMATION info2;
    memset( &info1, 0, sizeof(info1) );
    memset( &info2, 0, sizeof(info2) );

    int got_info = ( GetFileInformationByHandle( fh1, &info1 )
                     && GetFileInformationByHandle( fh2, &info2 ) ) ? 1 : 0;
    CloseHandle(fh1);
    CloseHandle(fh2);
    if ( !got_info )
      return 0;

    //always false for directories:
    if ( info1.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY )
      return 0;
    if ( info2.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY )
      return 0;

    //Now compare index:
    if ( info1.nFileIndexLow != info2.nFileIndexLow )
      return 0;
    if ( info1.nFileIndexHigh != info2.nFileIndexHigh )
      return 0;
    if ( info1.dwVolumeSerialNumber != info2.dwVolumeSerialNumber )
      return 0;
    //Extra checks, added for safety in cases of rare clashes in the above:
    if ( info1.nFileSizeLow != info2.nFileSizeLow )
      return 0;
    if ( info1.nFileSizeHigh != info2.nFileSizeHigh )
      return 0;
    if ( info1.nNumberOfLinks != info2.nNumberOfLinks )
      return 0;
    if ( info1.dwFileAttributes != info2.dwFileAttributes )
      return 0;
    //Finally, look at the creation time:
    if ( CompareFileTime( &info1.ftCreationTime,
                          &info2.ftCreationTime ) != 0 )
      return 0;

    //Must be the same file (even though the two paths might still represent two
    //separate files in the file-system, but hardlinked to the same actual file
    //node).
    return 1;
#else
    //A hopefully bullet proof way, (st_ino,st_dev) uniquely identifies a file
    //on a POSIX system.
    struct stat sinfo1;
    struct stat sinfo2;
    char buf[4096];
    {
      mcu8str native = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
      mcu8str_assign( &native, &p1 );
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
    if ( S_ISDIR(sinfo1.st_mode) )
      return 0;//always false for directories
    {
      mcu8str native = mcu8str_create_from_staticbuffer( buf, sizeof(buf) );
      mcu8str_assign( &native, &p2 );
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
    if ( S_ISDIR(sinfo2.st_mode) )
      return 0;//always false for directories
    return ( sinfo1.st_dev == sinfo2.st_dev
             && sinfo1.st_ino == sinfo2.st_ino ) ? 1 : 0;
#endif
  }

  mcu8str mctools_path_join( const mcu8str* p1raw, const mcu8str* p2raw )
  {
    //In general try to mimic os.path.path from Python, although we do not
    //require a drive letter on windows to be considered an absolute path, so we
    //might have slight changes in some weird cases.

    mcu8str p1 = mctools_impl_view_no_winnamespace(p1raw);
    mcu8str p2 = mctools_impl_view_no_winnamespace(p2raw);

#ifdef MC_IS_WINDOWS
    const char native_sep = '\\';
#else
    const char native_sep = '/';
#endif

    if ( p2.size == 0 ) {
      //Simply copy p1, but ensure it ends in a separator:
      if ( !p1.size )
        return mcu8str_create_empty();//special case: join("","")->""
      char lastchar = p1.c_str[p1.size-1];
      if ( lastchar == '\\' || lastchar == '/' ) {
        mcu8str p1copy = mcu8str_copy( &p1 );
        mctools_pathseps_platform(&p1copy);
        return p1copy;
      }
      mcu8str res = mcu8str_create( p1.size + 1 );
      mcu8str_append( &res, &p1 );
      res.c_str[res.size] = native_sep;
      ++( res.size );
      res.c_str[res.size] = '\0';
      mctools_pathseps_platform(&res);
      return res;
    }

    char drive_letter1 = p1.size >= 2 ? mctools_drive_letter(&p1) : 0;
    char drive_letter2 = p2.size >= 2 ? mctools_drive_letter(&p2) : 0;
    if ( p1.size == 0 || mctools_path_is_absolute(&p2) ) {
      //p1 is empty or p2 is an absolute path, so simply discard p1 (like
      //Pythons os.path.join does in this case). However, if there is a drive
      //letter in p1 and not in p2, we use that drive letter:
      if ( drive_letter1 && !drive_letter2 ) {
        mcu8str res = mcu8str_create( p2.size + 2 );
        res.c_str[0] = drive_letter1;
        res.c_str[1] = ':';
        res.c_str[2] = '\0';
        res.size = 2;
        mcu8str_append( &res, &p2 );
        mctools_pathseps_platform(&res);
        return res;
      } else {
        mcu8str p2copy = mcu8str_copy( &p2 );
        mctools_pathseps_platform(&p2copy);
        return p2copy;
      }
    }

    if ( drive_letter2 && drive_letter1 != drive_letter2 ) {
      //p2 has a drive letter and p1 has a different drive letter => simply
      //return p2.
      mcu8str p2copy = mcu8str_copy( &p2 );
      mctools_pathseps_platform(&p2copy);
      return p2copy;
    }

    //Neither p1 or p2 are empty here, so we can simply join them. However, if
    //there is a common drive letter or repeated slashes at the end of p1, we do
    //not repeat them needlessley.

    //Selected ranges:
    mcu8str_size_t iB1 = 0;
    mcu8str_size_t iE1 = p1.size;
    mcu8str_size_t iB2 = 0;
    mcu8str_size_t iE2 = p2.size;

    if ( drive_letter1 && drive_letter1 == drive_letter2 )
      iB2 += 2;//skip common drive letter from p2

    //Strip repeated trailing slashes down to one trailing slash:
    mcu8str_size_t minlen1 = ( drive_letter1 ? 3 : 1 );
    while ( iE1-iB1 > minlen1
            && ( p1.c_str[iE1-1]=='/' || p1.c_str[iE1-1]=='\\' )
            && ( p1.c_str[iE1-2]=='/' || p1.c_str[iE1-2]=='\\' ) )
      --iE1;

    mcu8str_size_t needs_slash = ( ( p1.c_str[iE1-1]=='/'
                                   || p1.c_str[iE1-1]=='\\' ) ? 0 : 1 );

    mcu8str_size_t s1 = iE1-iB1;
    mcu8str_size_t s2 = iE2-iB2;
    mcu8str_size_t newsize = s1 + s2 + needs_slash;
    mcu8str res =  mcu8str_create( newsize  );
    STDNS memcpy( res.c_str, p1.c_str, s1 );
    mcu8str_size_t used = s1;
    if ( needs_slash ) {
      res.c_str[used] = native_sep;
      ++used;
    }
    STDNS memcpy( res.c_str + used, p2.c_str, s2+1 );//+1 to include null char
    res.size = (unsigned)newsize;//cast ok due to mcu8str_create
    mctools_pathseps_platform(&res);
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

    if ( mctools_impl_has_winnamespace( path ) )
      itB += 4;//skip '\\?\'
    if ( mctools_drive_letter( path ) )
      itB += 2;//skip drive letter

    if ( itB == itE )
      return mcu8str_create_empty();

    //Find last path sep:
    const char * it = itE - 1;
    while ( it > itB && *it!='/' && *it!='\\' )
      --it;
    if ( *it=='/' || *it=='\\' )
      ++it;
    mcu8str_size_t bnsize = (mcu8str_size_t)(itE - it);
    if ( bnsize == 0 || ( bnsize==1 && *it=='.' ) )
      return mcu8str_create_empty();
    mcu8str res = mcu8str_create( bnsize );
    STDNS memcpy( res.c_str, it, bnsize + 1 );
    res.size = (unsigned)bnsize;//cast ok due to mcu8str_create
    return res;
  }

  const char * mctools_basename_view( const mcu8str* path )
  {
    const char * itB = path->c_str;
    const char * itE = path->c_str + path->size;
    if ( mctools_impl_has_winnamespace( path ) )
      itB += 4;//skip '\\?\'
    if ( mctools_drive_letter( path ) )
      itB += 2;//skip drive letter
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
    const char * itLastDot = NULL;
    while ( 1 ) {
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

  mcu8str mctools_expand_path( const mcu8str* praw )
  {
    mcu8str p = mctools_impl_view_no_winnamespace(praw);
    if ( p.size == 0 )
      return mcu8str_create_empty();
    mcu8str res;
    res.owns_memory = 0;
    res.size = 0;
    res.buflen = 0;
#ifdef MC_IS_WINDOWS
    mcwinstr wp = mc_path2longwpath( &p );
    res = mc_winstr_to_u8str(&wp);
    mc_winstr_dealloc(&wp);
#else
    if ( p.size >= 2
         && p.c_str[0] == '~'
         && ( p.c_str[1] == '/' || p.c_str[1] == '\\' ) ) {
      const char * home = getenv("HOME");
      if ( home ) {
        if ( p.size == 2 ) {
          mcu8str rhome = mcu8str_create_from_cstr(home);
          mctools_pathseps_platform(&rhome);
          return rhome;
        }
        const mcu8str_size_t home_size = mctools_strlen( home, 0 );
        const mcu8str_size_t newsize = (mcu8str_size_t)(home_size + p.size - 1);
        if ( ( newsize + 1 != (mcu8str_size_t)(home_size + p.size) )
             || newsize <= home_size || newsize <= p.size )
          mctools_impl_error("string length overflow in mctools_expand_path");
        res = mcu8str_create( newsize );
        mcu8str_append_cstr(&res,home);
        mcu8str_append_cstr(&res, p.c_str + 1 );
      }
    }
    if ( res.size == 0 )
      res = mcu8str_copy( &p );
#endif
    mctools_pathseps_platform(&res);
    return res;
  }

  mcu8str mctools_absolute_path( const mcu8str* pathraw )
  {
    mcu8str path = mctools_impl_view_no_winnamespace(pathraw);
    mcu8str res = mcu8str_create_empty();
    if ( path.size == 0 )
      return res;
#ifdef MC_IS_WINDOWS
    {
      mcwinstr wpath = mc_path2wpath( &path );
      res = mc_winstr_expand_to_fullpath_u8str( &wpath );
      mc_winstr_dealloc(&wpath);
    }
    const int not_done = (res.size == 0?1:0);
#else
    const int not_done = 1;
#endif
    if ( not_done ) {
      if ( mctools_path_is_absolute(&path) ) {
        res = mcu8str_copy( &path );
      } else {
        mcu8str cwd = mctools_get_current_working_dir();
        res = mctools_path_join( &cwd, &path );
        mcu8str_dealloc(&cwd);
      }
    }
    mctools_pathseps_platform(&res);
    return res;
  }

  mcu8str mctools_dirname( const mcu8str* pathraw )
  {
    mcu8str path = mctools_impl_view_no_winnamespace(pathraw);
    char drive_letter = mctools_drive_letter( &path );
    const char * itB = path.c_str;
    const char * itE = path.c_str + path.size;
    if ( drive_letter )
      itB += 2;//ignore drive letter while analysing

    if ( itB == itE ) {
      //empty:
      if ( drive_letter ) {
        mcu8str res = mcu8str_create_from_cstr("::");
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
    mcu8str_size_t ressize = (mcu8str_size_t)(it-it0);

    if ( ressize == 1 && *it0 == '.' ) {
      //special case, "." or "D:."
      if ( drive_letter ) {
        mcu8str res = mcu8str_create_from_cstr("::");
        res.c_str[0] = drive_letter;
        return res;
      } else {
        return mcu8str_create_from_cstr(".");
      }
    }

    if ( ressize == 0 ) {
      if ( path.c_str[0] == '.' && !drive_letter )
        return mcu8str_create_from_cstr(".");
      if ( drive_letter ) {
        mcu8str res = mcu8str_create_from_cstr("::");
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
      res.size = (unsigned)ressize+2;//cast ok, checked by mcu8str_create
      mctools_pathseps_platform(&res);
      return res;
    } else {
      mcu8str res = mcu8str_create( ressize );
      STDNS memcpy( res.c_str, it0, ressize );
      res.c_str[ressize] = '\0';
      res.size = (unsigned)ressize;//cast ok, checked by mcu8str_create
      mctools_pathseps_platform(&res);
      return res;
    }
  }

#ifdef _WIN32
  wchar_t* mctools_path2wpath( const mcu8str* path )
  {
    mcu8str p = mctools_impl_view_no_winnamespace(path);
    if ( p.size == 0 )
      return NULL;
    mcwinstr wp = mc_path2longwpath( &p );
    return wp.c_str;
  }

  mcu8str mctool_wcharstr_to_u8str( const wchar_t* wstr_raw )
  {
    mcwinstr wstr_view;
    wstr_view.c_str = (wchar_t*)wstr_raw;
    wstr_view.size = wcslen(wstr_raw);
    wstr_view.buflen = wstr_view.size + 1;
    const mcwinstr* wstrcp = &wstr_view;
    return mc_winstr_to_u8str( wstrcp );
  }

#endif

#ifdef MCFILEUTILS_CPPNAMESPACE
}
#endif

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

#include "NCrystal/internal/utils/NCDynLoader.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCFileUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include <functional>
#include <stdexcept>
#include <string>

//NCRYSTAL_DISABLE_DYNLOAD can be used to disable dynamic loading (controlled
//via CMake NCRYSTAL_ENABLE_DYNLOAD flag).

#if ( defined (_WIN32) || defined (WIN32) ) && !defined (__CYGWIN__) && !defined (NCRYSTAL_WIN_LOADLIB)
#  define NCRYSTAL_WIN_LOADLIB
#endif

#ifndef NCRYSTAL_WIN_LOADLIB
#  include <dlfcn.h>
#else
#  define NCRYSTAL_WIN_LOADLIB
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#endif

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
#ifndef NCRYSTAL_DISABLE_DYNLOAD
    std::mutex& getMutex() {
      static std::mutex theMutex;
      return theMutex;
    }
#endif
  }

  //Wrap dlsym. Returns { errormsg, address }. Success is indicated by errormsg.empty().
  std::pair<std::string,void *> implLookupSymbol( void* handle,  const std::string& symbol )
  {
#ifdef NCRYSTAL_DISABLE_DYNLOAD
    (void)handle;
    (void)symbol;
    return {std::string("NCrystal dynamic loading is disabled"),nullptr};
#else
    NCRYSTAL_LOCK_GUARD( getMutex() );
    const char* errMsg(nullptr);
    void* result;
#  ifdef NCRYSTAL_WIN_LOADLIB
    FARPROC fproc;
    result = (void *)(intptr_t)GetProcAddress((HINSTANCE)handle, symbol.c_str());
    if (!result) {
      errMsg = "GetProcAddress failed";
      //NB: We could perhaps useGetLastError() to get an error code and extract an
      //actual message with that.;
    }
#  else
    dlerror();//clear previous errors
    result = dlsym( handle, symbol.c_str());
    if ( !result ) {
      //NB: This might not be an error, functions are allowed to have address of
      //0x0! However, it might be an error, and a call to dlerror() will give a
      //nullptr in case there is not actually an error.
      errMsg = dlerror();
    }
#  endif
    if (errMsg) {
      std::string errstr(errMsg);
      if (errstr.empty())
        errstr = "<unknown>";
      return { errstr, nullptr };
    } else {
      return { std::string(), result };
    }
#endif
  }
}

NC::DynLoader::DynLoader( const std::string& filename,
                          ScopeFlag scopeflag,
                          LazyFlag lazyflag )
  : m_handle(nullptr), m_filename(filename)
{
  (void)lazyflag;
  (void)scopeflag;
#ifndef NCRYSTAL_DISABLE_DYNLOAD
  NCRYSTAL_LOCK_GUARD( getMutex() );
#  ifdef NCRYSTAL_WIN_LOADLIB
  m_handle = (void*)LoadLibrary(filename.c_str());
#  else
  dlerror();//clear previous errors
  int flag = ( scopeflag==ScopeFlag::global ? RTLD_GLOBAL : RTLD_LOCAL );
  flag |= ( lazyflag==LazyFlag::lazy ? RTLD_LAZY : RTLD_NOW );
  m_handle = dlopen( filename.c_str(), flag );
  if (!m_handle && !startswith(filename,"/")) {
    //Perhaps it failed because library path was relative. Let us try absolute:
    std::string abspathfilename = ncgetcwd()+"/"+filename;
    if (file_exists(abspathfilename))
      m_handle = dlopen( abspathfilename.c_str(), flag );
  }
#  endif
#endif
  if ( !m_handle ) {
#ifdef NCRYSTAL_DISABLE_DYNLOAD
    const char* errMsg = "NCrystal dynamic loading is disabled";
#else
#  ifdef NCRYSTAL_WIN_LOADLIB
    //    const char* errMsg = GetLastError();
    const char* errMsg = "";
#  else
    const char* errMsg = dlerror();
#  endif
#endif
    NCRYSTAL_THROW2( DataLoadError, "Could not load shared library: " << filename
                     << " (error was: " << ( errMsg ? errMsg : "<unknown>" ) << ")" );
  }
}

NC::DynLoader::~DynLoader()
{
  if ( !m_handle || !m_doClose )
    return;

  const char* errMsg = nullptr;
#ifdef NCRYSTAL_DISABLE_DYNLOAD
  errMsg = "NCrystal dynamic loading is disabled";
#else
  NCRYSTAL_LOCK_GUARD( getMutex() );
#  ifdef NCRYSTAL_WIN_LOADLIB
  if (!FreeLibrary((HINSTANCE)m_handle))
    errMsg = "";//tried errMsg = GetLastError();
#  else
  dlerror();//clear
  if ( dlclose(m_handle)!=0 )
    errMsg = dlerror();
#  endif
#endif
  if (errMsg) {
    //Never throw from destructors, but we hope it is at least safe to call the
    //message handlers.
    NCRYSTAL_WARN("Problems releasing handle to shared library: "
                  << m_filename << " (error was: "
                  << ( errMsg ? errMsg : "<unknown>" ) << ")");
  }
}

void * NC::DynLoader::findSymbolAddr( const std::string& symbol ) const
{
  auto result = implLookupSymbol( m_handle,  symbol );
  if ( !result.first.empty() )
    NCRYSTAL_THROW2( DataLoadError, "Problems looking up symbol \""<<symbol
                     <<"\" in shared library: " << m_filename
                     << " (error was: " << result.first << ")" );
  return result.second;
}

std::pair<bool,void *> NC::DynLoader::tryFindSymbolAddr( const std::string& symbol ) const
{
  auto result = implLookupSymbol( m_handle,  symbol );
  return { result.first.empty() , result.second };
}

void NC::DynLoader::doNotClose()
{
  m_doClose = false;
}

NC::DynLoader::DynLoader( NC::DynLoader&& o )
  : m_handle(o.m_handle),
    m_filename(std::move(o.m_filename)),
    m_doClose(o.m_doClose)
{
  o.m_handle = nullptr;
  o.m_filename.clear();
  o.m_doClose = false;
}

NC::DynLoader& NC::DynLoader::operator=( DynLoader&& o )
{
  m_doClose = o.m_doClose;
  m_handle = o.m_handle;
  m_filename = std::move(o.m_filename);
  m_doClose = o.m_doClose;
  o.m_handle = nullptr;
  o.m_filename.clear();
  o.m_doClose = false;
  return *this;
}

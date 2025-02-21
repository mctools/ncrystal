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

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// As an addition to the usual NCrystal license pasted above, note that THIS  //
// PARTICULAR FILE is also placed into the Public Domain, so that it might    //
// be easily adopted into the code-base of any project wishing to use it:     //
//                                                                            //
// This is free and unencumbered software released into the public domain.    //
//                                                                            //
// Anyone is free to copy, modify, publish, use, compile, sell, or            //
// distribute this software, either in source code form or as a compiled      //
// binary, for any purpose, commercial or non-commercial, and by any          //
// means.                                                                     //
//                                                                            //
// In jurisdictions that recognize copyright laws, the author or authors      //
// of this software dedicate any and all copyright interest in the            //
// software to the public domain. We make this dedication for the benefit     //
// of the public at large and to the detriment of our heirs and               //
// successors. We intend this dedication to be an overt act of                //
// relinquishment in perpetuity of all present and future rights to this      //
// software under copyright law.                                              //
//                                                                            //
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,            //
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF         //
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.     //
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR          //
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,      //
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR      //
// OTHER DEALINGS IN THE SOFTWARE.                                            //
//                                                                            //
// For more information, please refer to <https://unlicense.org/>             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "ncrystal_load.hh"
#include <string>
#include <mutex>
#include <stdexcept>
#include <stdio.h>

#if !defined(NCLOAD_WINDOWS) && ( defined (_WIN32) || defined (WIN32) )
#  define NCLOAD_WINDOWS
#endif
#ifdef NCLOAD_WINDOWS
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#else
#  include <dlfcn.h>
#endif

namespace {

  struct ncrystal_config {
    std::string shlibpath;
    unsigned long intversion = 0;
    std::string symbol_namespace;
  };

  ncrystal_config query_ncrystal_config()
  {
    char buffer[4096];
#ifdef NCLOAD_WINDOWS
    FILE* pipe = _popen("ncrystal-config --show "
                        "intversion shlibpath namespace", "r");
#else
    FILE* pipe = popen("ncrystal-config --show "
                       "intversion shlibpath namespace 2>/dev/null", "r");
#endif
    if (!pipe)
      return {};//failure
    auto readLine = [pipe]( std::string& tgt ) -> bool
    {
      //Read line and discard trailing whitespace (including newline chars).
      char buffer[4096];
      if (fgets(buffer, sizeof(buffer), pipe) == NULL)
        return false;
      tgt = buffer;
      while( !tgt.empty() && std::isspace( tgt.back() ) )
        tgt.pop_back();
      return true;
    };
    auto parseIntVersion = []( const std::string& s )
    {
      char * str_end = nullptr;
      unsigned long v = std::strtoul( s.c_str(), &str_end, 10 );
      return ( v >= 2002000
               && v < 999999999
               && str_end == s.c_str() + s.size() ) ? v : 0;
    };

    ncrystal_config res;
    bool all_ok(true);
    if ( !readLine(res.shlibpath)
         || !(res.intversion = parseIntVersion( res.shlibpath ))
         || !readLine( res.shlibpath )
         || res.shlibpath.empty()
         || !readLine(res.symbol_namespace ) ) {
      res.intversion = 0;//failure
    }

#ifdef NCLOAD_WINDOWS
    auto returnCode = _pclose(pipe);
#else
    auto returnCode = pclose(pipe);
#endif
    if ( returnCode == 0 && res.intversion >= 2002000  )
      return res;
    return {};//failure
  }

  struct NCrystalAPIDB {
    std::mutex mtx;
    std::shared_ptr<const NCrystalAPI> api;
    typedef void* (*FctSignature)(int);
    FctSignature ncrystal_access_dynapi_fct = nullptr;
  };

  void* load_dynapi_raw( unsigned interface_id, NCrystalAPIDB& db )
  {
    if ( !db.ncrystal_access_dynapi_fct ) {
      auto cfg = query_ncrystal_config();
      if ( !cfg.intversion >= 4000003 )
        throw std::runtime_error("Could not locate a functioning and"
                                 " recent enough NCrystal installation.");
#ifdef NCLOAD_WINDOWS
      auto handle = LoadLibrary(cfg.shlibpath.c_str());
#else
      dlerror();//clear previous errors
      void * handle = dlopen( cfg.shlibpath.c_str(), RTLD_LOCAL | RTLD_LAZY );
#endif
      if (!handle)
        throw std::runtime_error("Loading of the NCrystal library failed");

      std::string symbol("ncrystal");
      symbol += cfg.symbol_namespace;
      symbol += "_access_dynamic_api";

#ifdef NCLOAD_WINDOWS
      FARPROC fproc;
      void * addr = (void *)(intptr_t)GetProcAddress(handle,
                                                     symbol.c_str());
      if (!addr)
        throw std::runtime_error("GetProcAddress("
                                 "ncrystal_access_dynamic_api) failed");
#else
      dlerror();//clear previous errors
      void * addr = dlsym( handle, symbol.c_str());
      if ( !addr )
        throw std::runtime_error("dlsym(ncrystal_access_dynamic_api) failed");
#endif
      db.ncrystal_access_dynapi_fct
        = reinterpret_cast<NCrystalAPIDB::FctSignature>(addr);
    }

    return (*db.ncrystal_access_dynapi_fct)( interface_id );
  }

  NCrystalAPIDB& getNCrystalAPIDB()
  {
    static NCrystalAPIDB db;
    return db;
  }
}

std::shared_ptr<const NCrystalAPI> loadNCrystalAPI()
{
  auto& db = getNCrystalAPIDB();
  std::lock_guard<std::mutex> lock(db.mtx);
  if ( ! db.api ) {
    void * raw_api = load_dynapi_raw(NCrystalAPI::interface_id,db);
    if (!raw_api)
      throw std::runtime_error("Failed to access required NCrystal interface.");
    db.api = *reinterpret_cast<std::shared_ptr<const NCrystalAPI>*>(raw_api);
  }
  return db.api;
}

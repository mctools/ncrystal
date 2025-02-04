#ifndef NCrystal_DynLoader_hh
#define NCrystal_DynLoader_hh

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

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  //Dynamic load libraries and symbols.

  class DynLoader : private MoveOnly {
  public:

    //Constructor, getFunction, and findSymbolAddr functions throw DataLoadError
    //exceptions in case of problems.

    //Flags (same meaning as dlopen's RTLD_ flags). Might not be used on Windows:
    enum class ScopeFlag { global, local, default_value = local };
    enum class LazyFlag { lazy, now, default_value = lazy };

    //Constructor, loads library:
    DynLoader( const std::string& filename,
               ScopeFlag = ScopeFlag::default_value,
               LazyFlag = LazyFlag::default_value );

    //Destructor, unloads library unless doNotClose was called (or object was moved-from):
    ~DynLoader();

    //Lookup function by symbol and reinterpret_cast it to a given signature:
    template<class Signature>
    std::pair<bool,std::function<Signature>> tryGetFunction( const std::string& fctsymbol) const
    {
      auto a = tryFindSymbolAddr(fctsymbol);
      return { a.first, reinterpret_cast<Signature*>(a.second) };
    }

    template<class Signature>
    std::function<Signature> getFunction( const std::string& fctsymbol) const
    {
      //The use of dummy below, in the very particular manner it is written, is
      //a way to workaround a quirk in the C/POSIX specification. With gcc 4.8.5
      //and -pedantic, it would generate a misleading warning ("ISO C++ forbids
      //casting between pointer-to-function and pointer-to-object") otherwise.
      //The workaround is taken from the dlsym manpage
      //(https://linux.die.net/man/3/dlsym).
      void (*dummy)();
      *(void **) (&dummy) = findSymbolAddr(fctsymbol);
      return reinterpret_cast<Signature*>(dummy);
    }

    //Just lookup the address of symbols:
    void * findSymbolAddr( const std::string& ) const;
    std::pair<bool,void *> tryFindSymbolAddr( const std::string& ) const;

    DynLoader( DynLoader&& );
    DynLoader& operator=( DynLoader&& );

    //Call to avoid dlclose call in destructor:
    void doNotClose();

  private:
    void* m_handle;
    std::string m_filename;
    bool m_doClose = true;
  };

}

#endif

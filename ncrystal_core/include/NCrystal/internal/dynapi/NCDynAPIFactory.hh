#ifndef NCrystal_DynAPIFactory_hh
#define NCrystal_DynAPIFactory_hh

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

// Factory function ncrystal_access_dynamic_api function which can be used to
// instantiate a given dynamic api. This function has C-bindings since it is
// meant to be accessible via dlopen or GetProcAddress. Note also that it is
// namespaced when NCrystal is built with namespaced symbols (e.g. if the
// namespace is "dev", then the actual symbol to dlopen is
// "ncrystaldev_access_dynamic_api"). The intention is for client applications
// to use the ncrystal-config command to get both the path to the NCrystal
// shared library and the NCrystal symbol namespace (they might also opt to get
// the NCrystal version if they wish for whatever reason to impose runtime
// restrictions on the version of NCrystal).

#include "NCrystal/core/NCException.hh"

#ifdef ncrystal_access_dynamic_api
#  undef ncrystal_access_dynamic_api
#endif
#define ncrystal_access_dynamic_api NCRYSTAL_APPLY_C_NAMESPACE(access_dynamic_api)

extern "C" {
  NCRYSTAL_API void * ncrystal_access_dynamic_api( int interface_id );
}

//Convenience wrapper to be able to use the API directly from C++ code via
//normal build-time linking:
namespace NCRYSTAL_NAMESPACE {

  template<class TDynAPI>
  inline const TDynAPI * createDynAPI( bool allow_fail = false )
  {
    void * o = ncrystal_access_dynamic_api( TDynAPI::interface_id );
    if (o)
      return reinterpret_cast<const TDynAPI*>(const_cast<const void*>(o));
    if ( !allow_fail )
      NCRYSTAL_THROW2( BadInput,
                       "Unable to provide DynAPI with interface_id "
                       << TDynAPI::interface_id );
    return nullptr;

  }
}

#endif

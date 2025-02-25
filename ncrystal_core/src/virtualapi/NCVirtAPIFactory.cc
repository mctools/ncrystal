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

#include "NCrystal/virtualapi/NCVirtAPIFactory.hh"
#include "NCVirtAPI_Type1_v1_impl.hh"

void * ncrystal_access_virtual_api( unsigned interface_id )
{
  using t1v1 = NCrystal::VirtAPI::Type1_v1_Impl;
  if ( interface_id == t1v1::interface_id ) {
    static std::shared_ptr<const t1v1> sp_t1v1 = std::make_shared<t1v1>();
    return (void*)(&sp_t1v1);
  }

  return nullptr;
}

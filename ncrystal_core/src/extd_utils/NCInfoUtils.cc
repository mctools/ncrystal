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

#include "NCrystal/internal/extd_utils/NCInfoUtils.hh"

namespace NC = NCrystal;

const NC::DynamicInfo*
NC::InfoUtils::findDynInfo( const Info& info,
                            Optional<std::string> atomDisplayLabel )
{
  const DynamicInfo* res = nullptr;
  for ( auto& di : info.getDynamicInfoList() ) {
    if ( atomDisplayLabel.has_value()
         && ( info.displayLabel(di->atom().index)
              != atomDisplayLabel.value() ) ) {
      continue;
    }
    if (res) {
      NCRYSTAL_THROW(BadInput,"Could not uniquely determine DynInfo object."
                     " Perhaps providing an explicit atom label will help.");
    } else {
      nc_assert_always( di.get() != nullptr );
      res = di.get();
    }
  }
  if (!res)
    NCRYSTAL_THROW(BadInput,"Could not find requested DynInfo object.");
  return res;
}


#ifndef NCrystal_ParseNCMAT_hh
#define NCrystal_ParseNCMAT_hh

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

#include "NCrystal/internal/ncmat/NCNCMATData.hh"
#include "NCrystal/text/NCTextData.hh"

namespace NCRYSTAL_NAMESPACE {

  //Parse NCMAT data. Will throw BadInput exceptions in case of problems.
  //
  //It will do some rudimentary syntax checking (including presence/absence of
  //data sections and will call NCMATData::validate), but not a full validation
  //of data (a more complete validation is typically carried out afterwards by
  //the NCMAT Loader code).
  //
  //If doFinalValidation is false, the parser won't call NCMatData::validate()
  //before returning (although some other validations will still take place
  //during parsing). Calls with doFinalValidation=false should only happen if
  //the calling code intends to call NCMatData::validate() on the returned data.

  NCRYSTAL_API NCMATData parseNCMATData( const TextData&,
                                         bool doFinalValidation = true );

}

#endif

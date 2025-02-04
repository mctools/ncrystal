#ifndef NCrystal_LoadNCMAT_hh
#define NCrystal_LoadNCMAT_hh

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

#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/text/NCTextData.hh"
#include "NCrystal/factories/NCMatCfg.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace FactImpl {
    class InfoRequest;
  }
  class NCMATData;

  // Read .ncmat file and return a corresponding NCrystal::Info object from it
  // (it will have a reference count of 0 when returned).
  //
  // Parameters must be set via a NCMATCfgVars struct.  Parameters "temp",
  // "dcutoff", "dcutoffup", and "atomdb" have the same meaning as the
  // corresponding parameters described in NCMatCfg.hh (although atomdb must
  // here be already be split into "lines" and "words").
  //
  // Setting "temp" to -1.0 will result in a temperature of 293.15K unless
  // something in the input indicates another value (i.e. if a scatterkernel is
  // valid at 200K, then temp=-1 and temp=200 will both result in 200K, while
  // any other value results in an error).

  struct NCRYSTAL_API NCMATCfgVars : private MoveOnly {
    Temperature temp = Temperature{-1.0};//kelvin
    double dcutoff = 0.0;//angstrom
    double dcutoffup = kInfinity;//angstrom
    std::vector<VectS> atomdb;
    DataSourceName dataSourceName;

    //Optional (purely for passing on top-level info parameters to non-ncmat
    //sub-phases):
    const FactImpl::InfoRequest * originalInfoRequest = nullptr;
  };

  //The core feature is to load from parsed NCMAT data structure (which will be
  //consumed):
  NCRYSTAL_API Info loadNCMAT( NCMATData&& ncmat_data,
                               NCMATCfgVars&& cfgvars = NCMATCfgVars() );

  //For conveniece, can also parse TextData and load result:
  NCRYSTAL_API Info loadNCMAT( const TextData&,
                               NCMATCfgVars&& cfgvars = NCMATCfgVars() );

  //For additional convenience can also load TextData from filename, parse it,
  //and then load it:
  NCRYSTAL_API Info loadNCMAT( const char * ncmat_file,
                               NCMATCfgVars&& cfgvars = NCMATCfgVars() );

  NCRYSTAL_API Info loadNCMAT( const std::string& ncmat_file,
                               NCMATCfgVars&& cfgvars = NCMATCfgVars() );

  //Finally, it is of course possible to load directly from InfoRequest's, as
  //they, contain both configuration and TextData:
  NCRYSTAL_API Info loadNCMAT( const FactImpl::InfoRequest& );

  //@CUSTOM_xxx sections in NCMAT input will produce warnings by default, unless
  //the env var NCRYSTAL_NCMAT_NOWARNFORCUSTOM was set when NCrystal was
  //loaded. In any case, the behaviour can be queried and modified with:

  NCRYSTAL_API bool getNCMATWarnOnCustomSections();
  NCRYSTAL_API void setNCMATWarnOnCustomSections(bool);

}

#endif

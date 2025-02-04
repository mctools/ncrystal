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

#include "NCrystal/text/NCTextData.hh"
#include "NCrystal/factories/NCFactRequests.hh"
#include "NCrystal/internal/utils/NCAtomUtils.hh"

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Utilities for parsing the .laz/.lau files associated with the McStas       //
// components PowderN and Single_crystal, and turning them into Info objects. //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace Lazy {

    struct HKLFsq {
      double fsquared;
      int h,k,l;
    };

    struct ParsedLazyData : private MoveOnly {
      StructureInfo structInfo;//NB: unit cell volume not initialised yet
      DecodedChemForm chemform;
      Optional<Temperature> temp;
      DebyeTemperature debye_temp;
      VectS raw_header;
      std::vector<HKLFsq> hkllist;
    };

    struct LazyCfgVars : private MoveOnly {
      Optional<Temperature> temp;
      double dcutoff = 0.0;
      double dcutoffup = kInfinity;
      std::vector<VectS> atomdb;
      DataSourceName dataSourceName;
    };

    ParsedLazyData parseLazyTextData( const TextData&, const double& dcutoff=0 );
    InfoPtr buildInfo( const LazyCfgVars&, const ParsedLazyData& );
    InfoPtr buildInfoFromLazyData( const FactImpl::InfoRequest& );

    //Utility fct which could be useful for other purposes (i.e. if we want to
    //support hard-wired HKL lists in NCMAT files):
    std::vector<HKLFsq> validateAndNormaliseHKLFsqList( int spacegroup,
                                                        const std::vector<HKLFsq>& );
  }

}

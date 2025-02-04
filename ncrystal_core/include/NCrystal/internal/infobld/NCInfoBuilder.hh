#ifndef NCrystal_InfoBuilder_hh
#define NCrystal_InfoBuilder_hh

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

////////////////////////////////////////////////////////////////////////////
// Builder functions and associated mutable structures used to construct  //
// Info objects. Once constructed, Info objects are immutable.            //
////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace InfoBuilder {

    struct NCRYSTAL_API UnitCell final {
      StructureInfo structinfo;
      Optional<AtomInfoList> atomlist;
    };

    struct NCRYSTAL_API SinglePhaseBuilder final : MoveOnly {

      DataSourceName dataSourceName;

      //Most materials have at least unitcell or dynamics information (it is OK
      //to leave unitcell.structinfo.volume = 0.0 as it can be calculated):
      Optional<UnitCell> unitcell;
      Optional<DynamicInfoList> dynamics;

      //Composition (do not provide if unitcell.atomlist or dynamics are available)
      Optional<Info::Composition> composition;

      //Most materials also have a temperature:
      Optional<Temperature> temperature;

      //A density value must be provided if unitcell is absent (otherwise they
      //will be derived from the unit cell volume and structure). If one of the
      //two densities is provided, the other will be derived from it based on
      //calculated average atomic mass:
      Optional<Density> density;
      Optional<NumberDensity> numberDensity;

      struct NCRYSTAL_API HKLPlanes final : MoveOnly {
        PairDD dspacingRange;//Dspacing interval searched for planes (usually a
                              //bit bigger than that deduced from the list of
                              //planes since it is unlikely that planes exists
                              //with dspacing values exactly at the cutoff
                              //values).

        //Either provide the HKL list directly or via a generator function
        //allowing for delayed initialisation. Such a generator function will
        //receive ptrs to the unit cell info if one is specified - otherwise
        //just nullptrs. It will also get a requested dspacing range (usually
        //the same as the dcutoff_range above, but it might be a more narrow
        //range for inspections of the largest dspacings with less loading
        //time). Finally, the generator function should be thread-safe and might
        //be invoked more than once.
        using HKLListGenFct = std::function<HKLList(const StructureInfo*,
                                                    const AtomInfoList*,
                                                    PairDD dspacingRangeRequest)>;
        Variant<HKLList,HKLListGenFct> source;
      };
      Optional<HKLPlanes> hklPlanes;

      //Optional non-Bragg xsect curve (for esoteric usage only):
      std::function<CrossSect(NeutronEnergy)> bkgdxsectprovider;

      //Various other fields:
      Optional<Info::CustomData> customData;
      Info::StateOfMatter stateOfMatter = Info::StateOfMatter::Unknown;
    };

    struct NCRYSTAL_API MultiPhaseBuilder final : MoveOnly {
      Info::PhaseList phases;
    };

    //Build from provided information:
    NCRYSTAL_API Info buildInfo( SinglePhaseBuilder&& );
    NCRYSTAL_API Info buildInfo( MultiPhaseBuilder&& );

    //Build by modifying the density of an existing Info object (this is cheap
    //as heavy underlying data structures will be shared):
    NCRYSTAL_API InfoPtr buildInfoPtr( InfoPtr, Density );
    NCRYSTAL_API InfoPtr buildInfoPtr( InfoPtr, NumberDensity );
    NCRYSTAL_API InfoPtr buildInfoPtrWithScaledDensity( InfoPtr, double scaleFactor );

    //Record any non-Info parameters of cfg data on Info object:
    NCRYSTAL_API InfoPtr recordCfgDataOnInfoObject( InfoPtr, const Cfg::CfgData& );

    //Convenience:
    NCRYSTAL_API InfoPtr buildInfoPtr( SinglePhaseBuilder&& );
    NCRYSTAL_API InfoPtr buildInfoPtr( MultiPhaseBuilder&& );

    //Various utility functions follow here.

    //Construct composition from a chemical formula like "Al2O3", "H2O", or "Be":
    NCRYSTAL_API Info::Composition buildCompositionFromChemForm( const std::string& );

  }
}

#endif

#ifndef NCrystal_PowderBraggUtils_hh
#define NCrystal_PowderBraggUtils_hh

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

#include "NCrystal/internal/utils/NCExtraTypes.hh"
#include "NCrystal/interfaces/NCInfo.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace PowderBraggUtils {

    //Common utilities used for modelling Bragg diffraction in isotropic
    //materials. This is intended for usage by both ideal powder models or
    //models with extinction.

    using CellData = PowderBraggInput::CellData;
    using PIMergedData = PowderBraggInput::MergedData;
    using PIData = PowderBraggInput::Data;

    CellData prepareCellData( double volume, unsigned n_atoms );
    CellData prepareCellData( const StructureInfo& );

    ncnodiscard17 PIMergedData prepareMergedData( const Info& );
    ncnodiscard17 PIMergedData prepareMergedData( const StructureInfo&,
                                                  PIMergedData::PlaneList&& );
    ncnodiscard17 PIMergedData prepareMergedData( const CellData&,
                                                  PIMergedData::PlaneList&& );

    ncnodiscard17 PIData prepareData( const Info& );
    ncnodiscard17 PIData prepareData( const StructureInfo&,
                                      PIData::PlaneList&& );
    ncnodiscard17 PIData prepareData( const CellData&,
                                      PIData::PlaneList&& );

    ncnodiscard17 PIMergedData prepareMergedData( const PIData& );

    //Sort the plane list so larger d-spacing comes first.
    void sortData( PIMergedData& );
    void sortData( PIData& );
    void sortData( PIMergedData::PlaneList& );
    void sortData( PIData::PlaneList& );

    //NCBadInput if problems:
    void checkData( const PIMergedData& );
    void checkData( const PIData& );
    void checkData( const CellData& );

  }
}
#endif

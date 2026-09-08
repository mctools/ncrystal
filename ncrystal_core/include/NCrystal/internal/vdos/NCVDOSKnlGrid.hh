#ifndef NCrystal_VDOSKnlGrid_hh
#define NCrystal_VDOSKnlGrid_hh

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

#include "NCrystal/internal/vdos/NCVDOSGn.hh"

namespace NCRYSTAL_NAMESPACE {

  /////////////////////////////////////////////////////////////////////////
  // Utilities for determining suitable beta and alpha grids based on Gn //
  // functions.                                                          //
  /////////////////////////////////////////////////////////////////////////

  namespace VDOS {

    std::pair<unsigned,unsigned> gridDimFromLux( unsigned vdoslux );

    std::pair<VectD,VectD>
    determineAlphaBetaGridFromGn( const GnExpansion&,
                                  unsigned nalpha, unsigned nbeta );//fixme: nalpha before nbeta!

    //Intermediate function, representing a weighted combination of Gn
    //function. Returns (betavals,combgnvals), and is exposed here for
    //diagnostics:
    std::pair<VectD,VectD> getCombinedGnFct( const VDOSGn& Gn );

    //Intermediate function which determines alpha=beta grid points necessary to
    //model neutrons in the E->0 limit, where the available phasespace
    //approaches a narrowing region around the line alpha=beta. For reference it
    //also returns the relative contribution to S at E=0 at these points.
    std::pair<VectD,VectD> setupE0ABGrid( const GnExpansion&, unsigned npts );

  }
}

#endif

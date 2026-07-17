#ifndef NCrystal_SABCellInteg_hh
#define NCrystal_SABCellInteg_hh

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

#include "NCrystal/internal/sab/NCSABSurveyor.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"//fixme disentangle old infrastructure?
#include "NCrystal/internal/sab/NCSABCfg.hh"
#include "NCrystal/internal/utils/NCStrView.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {

    ////////////////////////////////////////
    // Utilities for SAB Cell integration //
    ////////////////////////////////////////

    struct CellData {
      double a1, a2; // alpha extent of cell, a2>a1>=0
      double b1, b2; // beta extent of cell, b2 > b1
      double S[4];   // {S(a1,b1), S(a2,b1), S(a1,b2), S(a2,b2).
      double logS[4];// {log(S[i])} (0 where S==0).
    };
    std::ostream& operator<<( std::ostream&, const CellData& );

    //Fixme: decide if better with a class or namespace for the next:
    struct StdLogLinCellIntegrator final : private NonInstantiable {

      ~StdLogLinCellIntegrator() = delete;

      using Region = SABCellSurvey::Region;
      using RegionList = SABCellSurvey::RegionList;

      using IntegrationScheme = SABCfg::IntegrationScheme;//fixme for migration

      ////////////////////////////////////////////////////////////////////////
      //Integrate within kinematic bounds.
      //
      //Important note: This is intended for usage in the case where we are
      //already certain that the phase space boundary CROSSES the cell.

      static void integrateWithinKB( const CellData&,
                                     double E_div_kT,
                                     IntegrationScheme,
                                     StableSumKahan& tgt );

      ////////////////////////////////////////////////////////////////////////
      //Integrate a full cell (this is a lot faster):

      static void integrateFullCell( const CellData& c, StableSumKahan& tgt )
      {
        double f = 0.5 * (c.b2-c.b1);
        tgt.add( f * integrateAlphaInterval_fast(c.a1,c.S[0],c.a2 , c.S[1],
                                                 c.logS[0], c.logS[1]) );
        tgt.add( f * integrateAlphaInterval_fast(c.a1,c.S[2], c.a2 , c.S[3],
                                                 c.logS[2], c.logS[3]) );
      }


    };

  }
}

#endif

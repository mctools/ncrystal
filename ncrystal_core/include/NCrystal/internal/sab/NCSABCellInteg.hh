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

      //Supported numerical integration schemes:
      enum class IntegrationScheme : std::uint_fast32_t {
        //Basic fixed-order schemes:
        Trapez2   = 0x020002,
        Trapez3   = 0x020003,
        Trapez5   = 0x020005,
        Trapez9   = 0x020009,
        Trapez17  = 0x020011,
        Trapez33  = 0x020021,
        Simpson3  = 0x040003,
        Simpson5  = 0x040005,
        Simpson9  = 0x040009,
        Simpson17 = 0x040011,
        Simpson33 = 0x040021,
        Romberg5  = 0x010005,
        Romberg9  = 0x010009,
        Romberg17 = 0x010011,
        Romberg33 = 0x010021,
        //Flex schemes, based on potentially adaptive Romberg integration:
        Flex5  = 0x100005,//mostly Romberg5, occasionally adaptive Romberg17+
                          //with target prec=1e-3
        Flex9  = 0x100009,//mostly Romberg9, occasionally adaptive Romberg17+
                          //with target prec=1e-4
        Flex17 = 0x100011,//mostly Romberg17, occasionally adaptive Romberg17+
                          //with target prec=1e-6
        Flex33 = 0x100021,//Adaptive Romberg33+, target prec=1e-9
        Flex65 = 0x100041,//Adaptive Romberg65+, target prec=1e-12
        MaxPrec = Flex65,
        Default = Flex17//fixme revisit
      };
      //Utilities for encoding to/from strings
      static const char * integSchemeToStr( IntegrationScheme );
      static IntegrationScheme str2IntegScheme( StrView );
      static const char * allIntegSchemesAsStr();//';' separated list

      ////////////////////////////////////////////////////////////////////////
      //Integrate within kinematic bounds.
      //
      //Important note: This is intended for usage in the case where we are
      //already certain that the phase space boundary CROSSES the cell.

      static void integrateWithinKB( const CellData&,
                                     double E_div_kT,
                                     IntegrationScheme,
                                     StableSum& tgt );

      ////////////////////////////////////////////////////////////////////////
      //Integrate a full cell (this is a lot faster):

      static void integrateFullCell( const CellData& c, StableSum& tgt )
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

#ifndef NCrystal_ExtnEval_hh
#define NCrystal_ExtnEval_hh

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

#include "NCrystal/internal/utils/NCExtraTypes.hh"

#if 0//fixme
namespace NCRYSTAL_NAMESPACE {

  ///////////////////////////////////////////////////////////////////////////
  //                                                                       //
  // SIMD-friendly infrastructure for extinction factor model evaluators.  //
  //                                                                       //
  ///////////////////////////////////////////////////////////////////////////

  namespace Extn {

    //Fixme: some of these simd baskets to common utils?
    //       ... and maybe we can use it for SC normal search speedup??
    //

    template<unsigned NCOUNT = 8>
    struct DblBasket {
      //Data structure optimized for SIMD, containing a fixed number of
      //doubles. The number should be low enough that any excess padding is
      //minimised, and high enough that it can easily fill SIMD
      //registers. AVX-512bit registers have space for 8 doubles, so it makes
      //sense to have at least 8. For now we aim for those 8, since we anyway
      //usually have non-vectorisable parts of the calculations as well.
      //
      //We maximally align the data to increase the likelihood of efficient
      //vectorisation.
      //
      //Note that in generic builds, like PyPI wheels, we normally can not have
      //AVX512 available, but might just have 128 bit simd registers. However,
      //coding like this will still make things faster + at some point even PyPI
      //wheels will catch up.
      static constexpr unsigned count = NCOUNT;
      alignas( alignof(double)*count ) double data[count];

      double& operator[]( unsigned idx ) ncnoexceptndebug {
#ifndef NDEBUG
        if ( i >= count )
          NCRYSTAL_THROW(CalcError,"DblBasket index out of range");
#endif
        return data[idx];
      }

      const double& operator[]( unsigned idx ) const ncnoexceptndebug {
#ifndef NDEBUG
        if ( i >= count )
          NCRYSTAL_THROW(CalcError,"DblBasket index out of range");
#endif
        return data[idx];
      }

    };

    template<class TDblBasket = DblBasket>
    struct FDMBasket {
      static constexpr unsigned count = TDblBasket::count;
      TDblBasket fsq;
      TDblBasket inv2dsp;// = 1/(2*dspacing)
      TDblBasket fdm;// = fsq * dsp * mult
    };

    std::vector<FDMBasket> vectorizeFDM( const PowderBraggInput::Data& );


  }

}

////////////////////////////
// Inline implementations //
////////////////////////////

//fixme: ...

#endif
#endif


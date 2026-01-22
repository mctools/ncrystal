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

#if 0//fixme

#include "NCrystal/internal/extn_utils/NCExtnEval.hh"

namespace NC = NCrystal;
namespace NCE = NCrystal::Extn;

std::vector<NCE::FDMBasket> NCE::vectorizeFDM( const PowderBraggInput::Data& data )
{
  std::vector<FDMBasket> res;
  unsigned i = FDMBasket::count;
  for ( auto& e : data.planes ) {
    if ( i >= FDMBasket::count ) {
      res.emplace_back();
      i = 0;
    }
    auto& b = res.back();
    b.fsq[i] = e.fsq;
    b.inv2dsp[i] = 1.0 / (2.0 * e.dsp);
    b.fdm[i] = e.fsq * e.dsp * e.mult;
    ++i;
  }
  //Fill remainder of last basket with dummy entries having Fsq=0:
  while ( i < FDMBasket::count ) {
    auto& b = res.back();
    b.fsq[i] = 0.0;
    b.inv2dsp[i] = 1.0;
    b.fdm[i] = 0.0;
    ++i;
  }
  return res;
}
#endif

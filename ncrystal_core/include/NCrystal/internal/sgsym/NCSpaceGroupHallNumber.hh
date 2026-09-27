#ifndef NCrystal_SpaceGroupHallNumber_hh
#define NCrystal_SpaceGroupHallNumber_hh

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

//NB: Only includes from the core component are allowed here, since this
//class is intended to eventually be moved to NCTypes.hh.
#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  //Identifies one of the 530 space group settings listed in ITVB 2001 Table
  //A1.4.2.7, by its index (1..530) in that table. This index is known as the
  //"Hall number" (e.g. in spglib). The class sgsym/NCSpaceGroup.hh provides
  //the corresponding space group information.

  class SpaceGroupHallNumber final
    : public EncapsulatedValue<SpaceGroupHallNumber,std::uint16_t> {
  public:
    using EncapsulatedValue::EncapsulatedValue;
    static constexpr const char * unit() noexcept { return ""; }
    static constexpr std::uint16_t max_value = 530;
    void validate() const;//throws BadInput if not in 1..530
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  inline void SpaceGroupHallNumber::validate() const
  {
    if ( ! ( m_value >= 1 && m_value <= max_value ) )
      NCRYSTAL_THROW2(BadInput,"Invalid space group Hall number: "<<m_value
                      <<" (must be in range 1.."<<max_value<<")");
  }
}

#endif

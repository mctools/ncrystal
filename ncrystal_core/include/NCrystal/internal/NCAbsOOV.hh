#ifndef NCrystal_AbsOOV_hh
#define NCrystal_AbsOOV_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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

#include "NCrystal/NCAbsorption.hh"


namespace NCrystal {

  class Info;

  class AbsOOV : public Absorption {
  public:

    // Provide absorption cross section based on simple 1/velocity
    // (=OneOverVelocity=OOV) scaling. This is non-oriented.

    //Constructor:
    AbsOOV(const Info*);
    virtual ~AbsOOV();

    bool isOriented() const { return false; };
    double crossSection(double ekin, const double (&neutron_direction)[3] ) const;
    double crossSectionNonOriented( double ekin ) const;

  private:
    double m_c;
  };
}

#endif

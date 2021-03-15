#ifndef NCrystal_SCOrientation_hh
#define NCrystal_SCOrientation_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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

#include "NCrystal/NCVariant.hh"
#include "NCrystal/NCTypes.hh"

namespace NCrystal {

  class NCRYSTAL_API SCOrientation {
  public:

    //Class which defines the orientation of a given single-crystal in the
    //laboratory frame, and which is needed to construct an SCBragg
    //instance.

    static constexpr const double default_tolerance = 1e-4;

    //Construct unspecified orientation, which must subsequently be completed
    //with calls to both setPrimaryDirection and setSecondaryDirection:
    SCOrientation() = default;

    //Provide vectors in the crystal frame and specify their desired direction
    //in the laboratory frame. This is well-formed as long as the opening angle
    //between primary and secondary vectors is similar in the two frames within
    //the stated tolerance in radians, because only the component of each
    //secondary direction ortogonal to the corresponding primary direction will
    //be used. The vectors do not need to be normalised:
    void setPrimaryDirection( const CrystalAxis&, const LabAxis& );
    void setSecondaryDirection( const CrystalAxis&, const LabAxis&,
                                double tolerance = default_tolerance );

    //Often, however, it is more convenient and useful to specify the directions
    //in the crystal frame by simply indicating a point in hkl space, and
    //letting the corresponding direction be the normal of the hkl-plane in
    //question, defined as:
    //
    // a_reci*h + b_reci*k + c_reci*l
    //
    //with (a_reci, b_reci, c_reci) being the reciprocal lattice indices. To
    //define directions in this fashion, use the following methods:

    void setPrimaryDirection( const HKLPoint&, const LabAxis& );
    void setSecondaryDirection( const HKLPoint&, const LabAxis&,
                                double tolerance = default_tolerance );

    //Complete when both primary and secondary directions have been set:
    bool isComplete() const;

    void clear();
  private:
    Variant<CrystalAxis,HKLPoint> m_crystal[2];
    Optional<LabAxis> m_lab[2];
    double m_tolerance = default_tolerance;
    struct Impl;
    friend struct Impl;
  public:
    //Access contents (idir=0: primary, idir=1: secondary):
    const Variant<CrystalAxis,HKLPoint>& getCrysDir(unsigned idir) const { nc_assert(idir<2); return m_crystal[idir]; }
    const Optional<LabAxis>& getLabDir(unsigned idir) const { nc_assert(idir<2); return m_lab[idir]; }
    double getTolerance() const { return m_tolerance; }
  };
}

#endif

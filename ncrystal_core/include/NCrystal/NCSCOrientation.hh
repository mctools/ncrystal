#ifndef NCrystal_SCOrientation_hh
#define NCrystal_SCOrientation_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

namespace NCrystal {

  class SCBragg;
  class MatCfg;

  class SCOrientation {
  public:

    //Class which defines the orientation of a given single-crystal in the
    //laboratory frame, and which is needed to construct an SCBragg
    //instance.

    //Construct unspecified orientation, which must subsequently be completed
    //with calls to both setPrimaryDirection and setSecondaryDirection:
    SCOrientation();

    //Provide vectors in the crystal frame and specify their desired direction
    //in the laboratory frame. This is well-formed as long as the opening angle
    //between primary and secondary vectors is similar in the two frames within
    //the stated tolerance in radians, because only the component of each
    //secondary direction ortogonal to the corresponding primary direction will
    //be used. The vectors do not need to be normalised:
    void setPrimaryDirection( const double (&crystal_direction)[3],
                              const double (&direction)[3] );
    void setSecondaryDirection( const double (&crystal_direction)[3],
                                const double (&direction)[3],
                                double tolerance = 1e-4 );

    //Often, however, it is more convenient and useful to specify the directions
    //in the crystal frame by simply indicating a point in hkl space, and
    //letting the corresponding direction be the normal of the hkl-plane in
    //question, defined as:
    //
    // a_reci*h + b_reci*k + c_reci*l
    //
    //with (a_reci, b_reci, c_reci) being the reciprocal lattice indices. To
    //define directions in this fashion, use the following methods:

    void setPrimaryDirection( double h, double k, double l,
                              const double (&direction)[3] );
    void setSecondaryDirection( double h, double k, double l,
                                const double (&direction)[3],
                                double tolerance = 1e-4 );

    //Complete when both primary and secondary directions have been set:
    bool isComplete() const;

    //Copy/assignment is allowed:
    SCOrientation(const SCOrientation&);
    SCOrientation& operator=(const SCOrientation&);

    ~SCOrientation();
  private:
    friend class SCBragg;
    friend class MatCfg;
    double m_crystal[2][3];
    double m_lab[2][3];
    double m_tolerance;
    bool m_crystal_is_hkl[2];
    void checkInput() const;
  };
}

#endif

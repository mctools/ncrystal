#ifndef NCrystal_SCGeoComputation_hh
#define NCrystal_SCGeoComputation_hh

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

#include "NCReflections.hh"

namespace NCrystal {

  class RotMatrix;

  class SCGeoComputation {
  public:

    //Constructor needs lattice spacings a,b,c and angles alpha, beta, gamma in degrees:
    SCGeoComputation( double a, double b, double c,
                      double alpha, double beta, double gamma );
    ~SCGeoComputation();

    Vector getReciDir(const Vector& hkl) const;
    Vector getReciDir(double h, double k, double l) const;

    //Calculates and caches transform:
    RotMatrix * calcTransform( Vector lab1, Vector lab2,
                               Vector crystal1, Vector crystal2 );

    //calculate the reciprocal vector in the crystal that has been aligned with lab using calcTransform
    //must only be called after calcTransform
    Vector getReciVecInRotCry(const Vector& reci_vec) const;

  protected:
    RotMatrix * m_cry2lab;
    RotMatrix * m_reci_lattice;
  };

}
#endif

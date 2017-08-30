#ifndef NCrystal_RotMatrix_hh
#define NCrystal_RotMatrix_hh

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

#include "NCMatrix.hh"
#include "NCVector.hh"

namespace NCrystal {

  class RotMatrix : public Matrix {
  public:

    //Specialised Matrix which always has dimensions of 3x3 and which can multiply
    //NCVector, as well as provide column access as NCVectors.

    virtual ~RotMatrix();

    RotMatrix();//null 3x3 matrix
    RotMatrix( const RotMatrix & );
    RotMatrix( const Matrix & );//must be dimension 3x3
    RotMatrix( const double * );//must be length 9
    RotMatrix( const std::vector<double>& );//must be length 9

    Vector operator*(const Vector&) const;
    const Vector& colX() const;
    const Vector& colY() const;
    const Vector& colZ() const;
  };

}

////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::Vector NCrystal::RotMatrix::operator*(const NCrystal::Vector& v) const
{
  nc_assert( m_colcount==3 && m_rowcount==3 );
  const double * d = &m_data[0];
  return Vector( d[0]*v.x() + d[1]*v.y() + d[2]*v.z(),
                 d[3]*v.x() + d[4]*v.y() + d[5]*v.z(),
                 d[6]*v.x() + d[7]*v.y() + d[8]*v.z()  );
}

inline const NCrystal::Vector& NCrystal::RotMatrix::colX() const
{
  return *reinterpret_cast<const NCrystal::Vector*>((*this)[0]);
}

inline const NCrystal::Vector& NCrystal::RotMatrix::colY() const
{
  return *reinterpret_cast<const NCrystal::Vector*>((*this)[1]);
}

inline const NCrystal::Vector& NCrystal::RotMatrix::colZ() const
{
  return *reinterpret_cast<const NCrystal::Vector*>((*this)[2]);
}

#endif

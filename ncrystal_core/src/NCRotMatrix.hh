#ifndef NCrystal_RotMatrix_hh
#define NCrystal_RotMatrix_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2018 NCrystal developers                                   //
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

    //Specialised Matrix which always has dimensions of 3x3 and which can
    //multiply NCVector, as well as provide column access as NCVectors. Note
    //that in the strict mathematical sense, the contents are not always
    //guaranteed to be a true rotation matrix.

    virtual ~RotMatrix();

    RotMatrix();//null 3x3 matrix

    //Construct rotation which transforms vectors v and u into v_trf and u_trf
    //respectively. As a sanity check, it is required (within the given
    //tolerance) that angle(v,u)==angle(v_trf,u_trf) and that v and u are not
    //parallel.
    RotMatrix(const Vector v, const Vector& v_trf,
              const Vector u, const Vector& u_trf, double tolerance = 1.0e-6 );

    double determinant() const;

    //Copy semantics:
    RotMatrix( const RotMatrix & );
    RotMatrix( const Matrix & );//must be dimension 3x3
    RotMatrix& operator=(const RotMatrix&);
    RotMatrix& operator=(const Matrix&);//must be dimension 3x3

#if __cplusplus >= 201103L
    //Move semantics:
    RotMatrix(RotMatrix && o);
    RotMatrix& operator=(RotMatrix&& o);
    RotMatrix(Matrix && o);//must be dimension 3x3
    RotMatrix& operator=(Matrix&& o);//must be dimension 3x3
#endif
    RotMatrix( const double * );//must be length 9
    RotMatrix( const std::vector<double>& );//must be length 9

    Vector operator*(const Vector&) const;
    const Vector& colX() const;
    const Vector& colY() const;
    const Vector& colZ() const;
  };

  class PhiRot {
  public:
    //Simpler class representing a simple rotation around a given axis, caching
    //just cosing and sine of the rotation angle, and rotating via Rodrigues'
    //rotation formula.
    PhiRot(double phi);//will calculate cosphi+sinphi, assuming phi is in [-pi,pi].
    PhiRot(double cosphi, double sinphi);//provide precalculated cos/sin values
    ~PhiRot();
    //Rotate v and angle of phi around an axis. The axis vector should be normalised, and for
    //efficiency, the function takes precalculated cross and dot products between the two vectors.
    Vector rotateVectorAroundAxis( const Vector& v, const Vector& axis,
                                   const Vector& axis_cross_v, double axis_dot_v,
                                   bool reverse_direction = false) const;
    //Version not using precalculated cross/dot products:
    Vector rotateVectorAroundAxis( const Vector& v, const Vector& axis,
                                   bool reverse_direction = false) const;
    double getCosPhi() const { return m_cosphi; }
    double getSinPhi() const { return m_sinphi; }
  private:
    double m_cosphi;
    double m_sinphi;
  };

  //With |a|=|b|=1, cosab=dot(a,b), sinab=sqrt(1-cosab^2), Rotate vector v from
  //the frame in which a=(sina,0,cosa), b=(0,0,1) into the frame in which they
  //have the provided values. If the two vectors are parallel (|sinab|~=0), the
  //rotation will not be fully determined, in which case the indeterminate
  //direction will be randomised if rand!=0, or if rand==0 an exception will be
  //thrown.
  class RandomBase;
  void rotateToFrame( double sinab, double cosab,
                      const Vector& a, const Vector& b,
                      Vector&v, RandomBase *  rand = 0);

}

////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCrystal {
  inline Vector RotMatrix::operator*(const Vector& v) const
  {
    nc_assert( m_colcount==3 && m_rowcount==3 );
    const double * d = &m_data[0];
    return Vector( d[0]*v.x() + d[1]*v.y() + d[2]*v.z(),
                   d[3]*v.x() + d[4]*v.y() + d[5]*v.z(),
                   d[6]*v.x() + d[7]*v.y() + d[8]*v.z()  );
  }

  inline const Vector& RotMatrix::colX() const
  {
    return *reinterpret_cast<const Vector*>((*this)[0]);
  }

  inline const Vector& RotMatrix::colY() const
  {
    return *reinterpret_cast<const Vector*>((*this)[1]);
  }

  inline const Vector& RotMatrix::colZ() const
  {
    return *reinterpret_cast<const Vector*>((*this)[2]);
  }

  inline RotMatrix::RotMatrix()
    : Matrix()
  {
    m_rowcount=3;
    m_colcount=3;
    m_data.resize(9,0.0);
  }

  inline RotMatrix::RotMatrix(const double *data)
    : Matrix(3,3,data)
  {
  }

  inline RotMatrix::RotMatrix(const std::vector<double>& data)
    : Matrix(3,3,data)
  {
    if (data.size()!=9)
      NCRYSTAL_THROW(BadInput,"Input vector to NCRotMatrix must have size 3*3=9");
  }

  inline RotMatrix::RotMatrix(const Matrix &m)
    : Matrix(m)
  {
    if (m.nRows()!=3||m.nCols()!=3)
      NCRYSTAL_THROW(BadInput,"Can only convert 3x3 Matrix to RotMatrix");
  }

  inline RotMatrix::RotMatrix(const RotMatrix &m)
    : Matrix(m)
  {
  }

    // RotMatrix& operator=(const RotMatrix&);
    // RotMatrix& operator=(const Matrix&);//must be dimension 3x3

  inline RotMatrix& RotMatrix::operator=(const Matrix& m)
  {
    if (m.nRows()!=3||m.nCols()!=3)
      NCRYSTAL_THROW(BadInput,"Can only convert 3x3 Matrix to RotMatrix");
    Matrix::operator=(m);
    return *this;
  }

  inline RotMatrix& RotMatrix::operator=(const RotMatrix& m)
  {
    Matrix::operator=(m);
    return *this;
  }
#if __cplusplus >= 201103L
  //Move semantics:
  inline RotMatrix::RotMatrix(RotMatrix && m)
    : Matrix(std::move(m))
  {
  }
  inline RotMatrix& RotMatrix::operator=(RotMatrix&& m)
  {
    Matrix::operator=(std::move(m));
    return *this;
  }
  inline RotMatrix::RotMatrix(Matrix && m)
    : Matrix(std::move(m))
  {
    if (m.nRows()!=3||m.nCols()!=3)
      NCRYSTAL_THROW(BadInput,"Can only convert 3x3 Matrix to RotMatrix");
  }
  inline RotMatrix& RotMatrix::operator=(Matrix&& m)
  {
    if (m.nRows()!=3||m.nCols()!=3)
      NCRYSTAL_THROW(BadInput,"Can only convert 3x3 Matrix to RotMatrix");
    Matrix::operator=(std::move(m));
    return *this;
  }
#endif

  inline PhiRot::PhiRot(double thephival) { sincos_mpipi(thephival,m_cosphi,m_sinphi); }
  inline PhiRot::PhiRot(double cp, double sp) : m_cosphi(cp), m_sinphi(sp) { nc_assert(ncabs(cp*cp+sp*sp-1.0)<1e-6); }
  inline PhiRot::~PhiRot(){}
  inline Vector PhiRot::rotateVectorAroundAxis( const Vector& v,
                                                const Vector& axis,
                                                const Vector& axis_cross_v,
                                                double axis_dot_v,
                                                bool reverse_direction ) const
  {
    nc_assert(axis.isUnitVector());
    Vector result = v;
    result *= m_cosphi;
    result += axis_cross_v * ( m_sinphi * (reverse_direction?-1.0:1.0) );
    result += axis * ( axis_dot_v*(1.0-m_cosphi) );
    return result;
  }

  inline Vector PhiRot::rotateVectorAroundAxis( const Vector& v,
                                                const Vector& axis,
                                                bool reverse_direction ) const
  {
    return rotateVectorAroundAxis(v,axis,axis.cross(v),axis.dot(v),reverse_direction);
  }

}

#endif

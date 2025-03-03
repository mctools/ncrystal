////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
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

#include "NCrystal/internal/utils/NCRotMatrix.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"

NCrystal::RotMatrix::RotMatrix(const NCrystal::Vector v, const NCrystal::Vector& v_trf,
                               const NCrystal::Vector u, const NCrystal::Vector& u_trf,
                               double tolerance )
  : Matrix()
{
  if (v.isParallel(u,tolerance))
    NCRYSTAL_THROW(BadInput,"v and u are parallel in attempt to construct RotMatrix(v,v_trf,u,u_trf)");
  Vector v1 = v.unit();
  Vector u1 = u.unit();
  Vector v2 = v_trf.unit();
  Vector u2 = u_trf.unit();
  double dot1 = v1.dot(u1);
  double dot2 = v2.dot(u2);
  if (ncabs(dot1-dot2) > tolerance)
    NCRYSTAL_THROW(BadInput,"angle(v,u)!=angle(v_trf,u_trf) in attempt to construct RotMatrix(v,v_trf,u,u_trf)");

  Vector cross2 = v2.cross(u2);
  Vector cross1 = v1.cross(u1);

  double M[18]={ v2.x(), u2.x(), cross2.x(), 1, 0, 0,
                 v2.y(), u2.y(), cross2.y(), 0, 1, 0,
                 v2.z(), u2.z(), cross2.z(), 0, 0, 1 };
  Matrix aMat(3,6,M);
  aMat.rref();
  double bmat_data[] = { aMat[0][3], aMat[0][4], aMat[0][5],
                         aMat[1][3], aMat[1][4], aMat[1][5],
                         aMat[2][3], aMat[2][4], aMat[2][5] };
  RotMatrix b_matrix( bmat_data );
  double lab_matrix_data[] = { v1.x(), u1.x(), cross1.x(),
                               v1.y(), u1.y(), cross1.y(),
                               v1.z(), u1.z(), cross1.z() };
  RotMatrix lab_matrix( lab_matrix_data);

  RotMatrix result( lab_matrix * b_matrix );
  m_rowcount=3;
  m_colcount=3;
  std::swap(m_data,result.m_data);
  nc_assert(ncabs(determinant()-1)<tolerance*10);
}

double NCrystal::RotMatrix::determinant() const
{
  nc_assert_always(m_data.size()==9);
  const double * m = &m_data[0];
#if 1
  StableSum sum;
  sum.add(  m[0]*m[4]*m[8] );
  sum.add( -m[0]*m[5]*m[7] );
  sum.add(  m[1]*m[5]*m[6] );
  sum.add( -m[1]*m[3]*m[8] );
  sum.add(  m[2]*m[3]*m[7] );
  sum.add( -m[2]*m[4]*m[6] );
  return sum.sum();
#else
  return ( m[0]*(m[4]*m[8]-m[5]*m[7])
           + m[1]*(m[5]*m[6]-m[3]*m[8])
           + m[2]*(m[3]*m[7]-m[4]*m[6]) );
#endif
}

void NCrystal::rotateToFrame( double sinab, double cosab, const Vector& a, const Vector& b, Vector&v, RNG * rng )
{
  nc_assert(v.isUnitVector(1e-3));
  nc_assert( a.isUnitVector() && b.isUnitVector() );
  nc_assert(ncabs(cosab-a.dot(b))<1e-6);

  if (ncabs(sinab)<1e-10) {
    if (!rng)
      NCRYSTAL_THROW(CalcError,"rotateToFrame called with parallel vectors so rotation is not fully specified.")
    //a and b are essentially parallel, the rotation is not well defined! First,
    //we apply the rotation which takes (0,0,1) into b=(bx,by,bz), by rotating
    //around crossproduct((0,0,1),(bx,by,bz))=(by,-bx,0):
    PhiRot phirot_ztob(b.z(),-std::sqrt(1.0-b.z()*b.z()));
    Vector axis_ztob(b.y(),-b.x(),0.);
    double m2 = axis_ztob.mag2();
    if (m2>1e-12) {
      axis_ztob *= 1.0/std::sqrt(m2);
      v = phirot_ztob.rotateVectorAroundAxis( v, axis_ztob );
      nc_assert((phirot_ztob.rotateVectorAroundAxis( Vector(0,0,1), axis_ztob )-b).mag2()<1e-8);
    } else {
      //Doubly degenerate! vector b is parallel to z. If b=(0,0,1) leave v
      //alone, otherwise b=(0,0,-1) reflect in zy plane.
      if (b.z()<0.0)
        v.z() *= -1.0;
    }

    //Next, perform a random rotation around b:
    double rand_cosphi, rand_sinphi;
    std::tie(rand_cosphi, rand_sinphi) = randPointOnUnitCircle( *rng );
    PhiRot phirot_random(rand_cosphi,rand_sinphi);
    v = phirot_random.rotateVectorAroundAxis( v, b );

    //remove numerical errors in normalisation:
    nc_assert(v.isUnitVector(1e-3));
    v.normalise();
    return;

  }

  //We are dealing with the typical non-degenerate case of non-parallel vectors
  //a and b. It is possible to show that the rotation matrix, R, is given by the
  //three column vectors:
  //
  //  R = {   (a-cosab*b)/sinab   |   crossproduct(b,a)/sinab   |   b   }
  //
  //R can for instance be derived by demanding R*(0,0,1)=b, R*(sinab,0,cosab)=a
  //and R*(0,sinab,0)=crossproduct(b,a).

  double s = 1.0/sinab;

  //construct the column vectors, arranging operations so as to hopefully avoid
  //needless temporaries:
  Vector col1(b); col1 *= (-cosab); col1 += a; col1 *= s;
  Vector col2(b.cross(a)); col2 *= s;
  //col3 is just b

  v.set( v.x()*col1.x()+v.y()*col2.x()+v.z()*b.x(),
         v.x()*col1.y()+v.y()*col2.y()+v.z()*b.y(),
         v.x()*col1.z()+v.y()*col2.z()+v.z()*b.z() );

  //remove numerical errors in normalisation:
  nc_assert(v.isUnitVector(1e-3));
  v.normalise();
}

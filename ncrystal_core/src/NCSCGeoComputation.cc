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

#include "NCSCGeoComputation.hh"
#include "NCVector.hh"
#include "NCRotMatrix.hh"
#include "NCMatrix.hh"
#include <cmath>

NCrystal::SCGeoComputation::SCGeoComputation(double a, double b, double c,
                                             double alpha, double beta, double gamma)
  : m_cry2lab(0),m_reci_lattice(0)
{
  alpha *= M_PI/180.;
  beta  *= M_PI/180.;
  gamma *= M_PI/180.;

  double uv = a*b*c*sqrt( 1.0 - cos(alpha)*cos(alpha)
                               - cos(beta)*cos(beta)
                               - cos(gamma)*cos(gamma)
                               + 2.0*cos(alpha)*cos(beta)*cos(gamma) );
  std::vector<double> lattice(9,0.);
  lattice[0]=a;
  lattice[2]=c*cos(beta);
  lattice[3]= b*cos(gamma);
  lattice[4]= b*sin(gamma);
  lattice[5]=-c*(cos(beta)*cos(gamma)-cos(alpha))/sin(gamma);
  lattice[8]=uv/(a*b*sin(gamma));


  m_reci_lattice= new RotMatrix(lattice);
  m_reci_lattice->inv();
  *m_reci_lattice*=2.*M_PI;

}

NCrystal::SCGeoComputation::~SCGeoComputation()
{
  delete m_cry2lab;
  delete m_reci_lattice;
}

NCrystal::Vector NCrystal::SCGeoComputation::getReciDir(const NCrystal::Vector & hkl) const
{
  return *m_reci_lattice*hkl;
}

NCrystal::Vector NCrystal::SCGeoComputation::getReciDir(double h, double k, double l) const
{
  return getReciDir(Vector(h,k,l));
}

NCrystal::Vector NCrystal::SCGeoComputation::getReciVecInRotCry(const NCrystal::Vector& cry) const
{
  nc_assert(m_cry2lab);
  return (*m_cry2lab*cry).unit();
}

NCrystal::RotMatrix * NCrystal::SCGeoComputation::calcTransform( NCrystal::Vector lab1,
                                                                 NCrystal::Vector lab2,
                                                                 NCrystal::Vector crystal1,
                                                                 NCrystal::Vector crystal2)
{
  crystal1.unit();
  crystal2=crystal2.unit();
  Vector crystal3=crystal1.cross(crystal2);

  lab1=lab1.unit();
  lab2=lab2.unit();
  Vector lab3= lab1.cross(lab2);

  double M[18]={crystal1.x(), crystal2.x(), crystal3.x(), 1, 0, 0,
      crystal1.y(), crystal2.y(), crystal3.y(), 0, 1, 0,
      crystal1.z(), crystal2.z(), crystal3.z(), 0, 0, 1};
  Matrix aMat(3,6,M);
  aMat.rref();
  double bmat_data[] = {aMat[0][3],aMat[0][4],aMat[0][5],
      aMat[1][3],aMat[1][4],aMat[1][5],
      aMat[2][3],aMat[2][4],aMat[2][5]};
  RotMatrix b_matrix(bmat_data );
  double lab_matrix_data[] = {lab1.x(), lab2.x(), lab3.x(),
      lab1.y(), lab2.y(), lab3.y(),
      lab1.z(), lab2.z(), lab3.z()};
  RotMatrix lab_matrix( lab_matrix_data);

  delete m_cry2lab;
  m_cry2lab = new RotMatrix(lab_matrix *b_matrix);
  return m_cry2lab;
}

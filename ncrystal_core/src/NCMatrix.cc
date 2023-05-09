////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2023 NCrystal developers                                   //
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

#include "NCrystal/internal/NCMatrix.hh"
#include <iostream>

NCrystal::Matrix::~Matrix()
{
}

NCrystal::Matrix NCrystal::operator*(const NCrystal::Matrix& m1, const NCrystal::Matrix& m2)
{
  nc_assert( m1.m_colcount==m2.m_rowcount) ;
  Matrix res;
  res.m_rowcount = m1.m_rowcount;
  res.m_colcount = m2.m_rowcount;
  res.m_data.reserve_hint(res.m_rowcount*res.m_colcount);
  for (unsigned i = 0; i < m1.m_rowcount; ++i) {
    for (unsigned j = 0; j < m2.m_colcount; ++j) {
      double temp = 0.0;
      for (unsigned k = 0; k < m1.m_colcount; ++k)
        temp += m1.m_data[i * m1.m_colcount + k] * m2.m_data[k * m2.m_colcount + j];
      res.m_data.push_back(temp);
    }
  }
  return res;
}

void NCrystal::Matrix::rref(double epsilon)
{
  unsigned r = 0;
  for (unsigned c = 0; c < m_colcount && r < m_rowcount; ++c) {
    int j = r;
    for (unsigned i = r+1; i < m_rowcount; ++i) {
      if (ncabs((*this)[i][c]) > ncabs((*this)[j][c]))
        j = i;
    }
    if (ncabs((*this)[j][c]) < epsilon)
      continue;
    std::swap_ranges(m_data.begin() + j*m_colcount,
        m_data.begin() + j*m_colcount + m_colcount,
        m_data.begin() + r*m_colcount);
    double s = 1.0 / (*this)[r][c];
    for (unsigned jj = 0; jj < m_colcount; ++jj)
      m_data[r*m_colcount + jj] *= s;
    for (unsigned i = 0; i < m_rowcount; ++i) {
      if (i != r) {
        double t = (*this)[i][c];
        for (unsigned jj = 0; jj < m_colcount; ++jj)
          m_data[i*m_colcount + jj] -= t * (*this)[r][jj];
      }
    }
    ++r;
  }
}


void NCrystal::Matrix::inv(double epsilon  )
{
  if( m_colcount!= m_rowcount)
    NCRYSTAL_THROW(CalcError,"inv: asking inverse matrix for a non-square matrix.");

  unsigned tmp_colcount = m_colcount*2;
  decltype(m_data) new_data;
  const auto newsize = m_rowcount*tmp_colcount;
  new_data.reserve_hint(newsize);
  for ( unsigned i = 0; i < newsize; ++i )
    new_data.emplace_back();

  for (unsigned i = 0; i < m_rowcount; ++i) {
    for (unsigned j = 0; j < m_colcount ; ++j) {
      new_data[i * tmp_colcount + j] = m_data[i * (m_colcount) + j];
    }
  }

  //make identity matrix
  for (unsigned j = 1; j <= m_colcount; ++j) {
    new_data[(m_rowcount-j) * tmp_colcount + tmp_colcount-j] = 1.;
  }

  new_data.swap(m_data);
  m_colcount *=2;
  rref(epsilon);
  m_colcount /=2;
  new_data.swap(m_data);

  for (unsigned i = 0; i < m_rowcount; ++i) {
    for (unsigned j = 0; j < m_colcount ; ++j) {
      m_data[i * (m_colcount) + j]= new_data[i * (tmp_colcount) + j+ m_colcount] ;
    }
  }

}


NCrystal::Matrix NCrystal::Matrix::getInv(double epsilon) const
{
  Matrix mat( MatrixAllowCopy, *this );
  mat.inv(epsilon);
  return mat;
}


std::ostream& NCrystal::operator << (std::ostream & o, const NCrystal::Matrix& matrix)
{
  o << " {\n";
  for(unsigned i=0;i<matrix.m_rowcount;i++) {
    o << " {";
    for(unsigned j=0;j<matrix.m_colcount;j++)
      o << " " << matrix[i][j]  ;
    o << " }\n";
  }
  o << " }\n";
  return o;
}

void NCrystal::Matrix::print() const
{
  std::cout<< *this <<std::flush;
}

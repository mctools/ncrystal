#ifndef NCrystal_Matrix_hh
#define NCrystal_Matrix_hh

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

#include "NCrystal/NCDefs.hh"
#include "NCrystal/internal/NCMath.hh"
#include <cstring>
#include <ostream>

namespace NCrystal {
  class Matrix {

  public:

    Matrix() ;

    Matrix(unsigned row, unsigned col, const double* data= 0);
    Matrix(unsigned row, unsigned col, const VectD& data);

    virtual ~Matrix();

    //Copy semantics:
    Matrix(const Matrix &) ;
    Matrix& operator=(const Matrix&);

    //Move semantics:
    Matrix(Matrix && o);
    Matrix& operator=(Matrix&& o);

    Matrix& operator*=(const double&);
    const double* operator[](unsigned irow) const;
    friend Matrix operator*(const Matrix&, const Matrix&);
    Matrix operator* (double );
    Matrix operator/ (double );
    bool operator== (const Matrix& );
    bool operator!= (const Matrix& );
    friend std::ostream& operator << (std::ostream &, const Matrix&);

    //Matrix transposition:
    Matrix operator~() const;
    void transpose();

    void print();
    void unit();
    void rref(double epsilon = 1e-5);//bring to reduced row echelon form
    void inv(double epsilon = 1e-5);
    Matrix getInv(double epsilon = 1e-5);


    unsigned nRows() const;
    unsigned nCols() const;
    double mag2() const;//sum of elements squared
    double mag() const { return sqrt(mag2()); }

  protected:

    VectD m_data;
    unsigned m_rowcount;
    unsigned m_colcount;
    //set(..) should not be public, otherwise NCRotMatrix won't be able to
    //guarantee it is always 3x3 (unless we make set(..) virtual):
    void set(unsigned row, unsigned col, const double* data);
  };

  std::ostream& operator << (std::ostream &, const Matrix&);
  Matrix operator*(const Matrix& m1, const NCrystal::Matrix& m2);

}


////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::Matrix::Matrix()
  : m_rowcount(0), m_colcount(0)
{
}

inline NCrystal::Matrix::Matrix(const NCrystal::Matrix & o)
{
  *this = o;
}

inline NCrystal::Matrix::Matrix(NCrystal::Matrix && o)
  : m_data(std::move(o.m_data)),
    m_rowcount(o.m_rowcount),
    m_colcount(o.m_colcount)
{
}

inline NCrystal::Matrix& NCrystal::Matrix::operator=(NCrystal::Matrix&& o)
{
  m_data = std::move(o.m_data);
  m_rowcount = o.m_rowcount;
  m_colcount = o.m_colcount;
  return *this;
}

inline void NCrystal::Matrix::set(unsigned row, unsigned col, const double* data)
{
  unsigned n = row*col;
  if ( !n ) {
    //must allow this, since the default constructor also gives row=col=0.
    if ( row || col )
      NCRYSTAL_THROW(BadInput,
                     "NCMatrix number of rows and columns must both be positive or both zero");
    m_rowcount = m_colcount = 0;
    m_data.clear();
  } else {
    nc_assert(data);
    m_data.resize(n);
    std::memcpy(&m_data[0],data,n*sizeof(double));
    m_rowcount = row;
    m_colcount = col;
  }
}

inline NCrystal::Matrix::Matrix(unsigned row, unsigned col, const double* data)
{
  set( row, col, data );
}

inline NCrystal::Matrix::Matrix(unsigned row, unsigned col, const VectD& data) {
  if ( !row*col == data.size() )
    NCRYSTAL_THROW(BadInput,"NCMatrix constructor got inconsistent data length");
  set( row, col, data.empty() ? 0 : &data[0] );
}

inline NCrystal::Matrix& NCrystal::Matrix::operator=(const NCrystal::Matrix& o)
{
  set(o.m_rowcount,o.m_colcount,o.m_data.empty() ? 0 : &(o.m_data[0]));
  return *this;
}

inline const double* NCrystal::Matrix::operator[](unsigned irow) const
{
  nc_assert( irow < m_rowcount) ;
  return &(m_data[0]) + irow * m_colcount;
}

inline NCrystal::Matrix NCrystal::Matrix::operator* (double factor)
{
  Matrix res;
  res.m_colcount = m_colcount;
  res.m_rowcount = m_rowcount;
  if (m_data.empty())
    return res;
  res.m_data.resize(m_data.size(),factor);
  double * itres = &res.m_data[0];
  double * it = &m_data[0];
  double * itresE = itres + res.m_data.size();
  for (;itres!=itresE;++it,++itres)
    *itres *= *it;
  return res;
}

inline NCrystal::Matrix NCrystal::Matrix::operator/ (double factor)
{
  return *this * ( 1.0 / factor );
}

inline bool NCrystal::Matrix::operator== (const NCrystal::Matrix& o)
{
  nc_assert(m_data.size()==o.m_data.size());
  return m_data == o.m_data;
}

inline bool NCrystal::Matrix::operator!= (const NCrystal::Matrix& o)
{
  return !( *this == o );
}

inline NCrystal::Matrix NCrystal::Matrix::operator~() const
{
  Matrix res;
  res.m_colcount = m_colcount;
  res.m_rowcount = m_rowcount;
  if (m_data.empty())
    return res;
  res.m_data.reserve(m_data.size());
  for (unsigned i = 0; i < m_colcount; i++)
    for (unsigned j = 0; j < m_rowcount; j++)
      res.m_data.push_back(m_data[j * m_colcount + i]);
  return res;
}

inline void NCrystal::Matrix::transpose()
{
  //TODO: Should be able to do this without malloc.
  VectD v;
  std::swap(m_data, v);
  m_data.reserve(v.size());
  for (unsigned i = 0; i < m_colcount; i++)
    for (unsigned j = 0; j < m_rowcount; j++)
      m_data.push_back(v[j * m_colcount + i]);
  std::swap(m_rowcount, m_colcount);
}

inline unsigned NCrystal::Matrix::nRows() const
{
  return m_rowcount;
}

inline unsigned NCrystal::Matrix::nCols() const
{
  return m_colcount;
}

inline double NCrystal::Matrix::mag2() const
{
  double sum(0.0);
  VectD::const_iterator it = m_data.begin();
  VectD::const_iterator itE = m_data.end();
  for (;it!=itE;++it)
    sum += (*it) * (*it);
  return sum;
}

inline NCrystal::Matrix& NCrystal::Matrix::operator*=(const double& f)
{
  VectD::iterator it = m_data.begin();
  VectD::iterator itE = m_data.end();
  for (;it!=itE;++it)
    *it *= f;
  return *this;
}

inline void NCrystal::Matrix::unit()
{
  double themag2 = mag2();
  if (themag2==1.0)
    return;
  if (!themag2)
    NCRYSTAL_THROW(CalcError,"NCMatrix::unit(): Can't scale null-matrix.");
  *this *= 1.0/sqrt(themag2);
}

#endif

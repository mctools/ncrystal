#ifndef NCrystal_Matrix_hh
#define NCrystal_Matrix_hh

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

#include "NCrystal/internal/utils/NCMath.hh"

namespace NCRYSTAL_NAMESPACE {


  struct MatrixAllowCopy_t {};
  constexpr MatrixAllowCopy_t MatrixAllowCopy = MatrixAllowCopy_t{};

  class Matrix : private MoveOnly {
  public:

    Matrix() ;

    // Matrix(unsigned row, unsigned col, const double* data= 0);
    Matrix(unsigned nrows, unsigned ncols, Span<const double> );

    virtual ~Matrix();

    //Protected copy semantics:
    Matrix( MatrixAllowCopy_t, const Matrix &);

    //Move semantics:
    Matrix(Matrix && o);
    Matrix& operator=(Matrix&& o);

    Matrix& operator*=(const double&);
    const double* operator[](unsigned irow) const;
    friend Matrix operator*(const Matrix&, const Matrix&);
    Matrix operator* (double );
    Matrix operator/ (double );

    bool isEqual(const Matrix&) const;
    friend std::ostream& operator << (std::ostream &, const Matrix&);

    //Matrix transposition:
    Matrix operator~() const;
    void transpose();

    void unit();
    void rref(double epsilon = 1e-5);//bring to reduced row echelon form
    void inv(double epsilon = 1e-5);
    Matrix getInv(double epsilon = 1e-5) const;

    unsigned nRows() const;
    unsigned nCols() const;
    double mag2() const;//sum of elements squared
    double mag() const { return std::sqrt(mag2()); }

    const SmallVector<double,9>& rawData() const { return m_data; }

  protected:

    //Use small vector optimisation for data (NSMALL=9 is optimised for usage as
    //3x3 RotMatrix):
    SmallVector<double,9> m_data;
    unsigned m_rowcount;
    unsigned m_colcount;
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

inline NCrystal::Matrix::Matrix(MatrixAllowCopy_t, const NCrystal::Matrix & o)
  : m_data(SVAllowCopy,o.m_data), m_rowcount(o.m_rowcount), m_colcount(o.m_colcount)
{
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

inline NCrystal::Matrix::Matrix(unsigned row, unsigned col, Span<const double> data)
{
  const auto n = row*col;
  if ( data.size() != n )
    NCRYSTAL_THROW(BadInput,
                   "NCMatrix number of rows and columns not consistent with data length");

  if ( !n ) {
    //must allow this, since the default constructor also gives row=col=0.
    if ( row || col )
      NCRYSTAL_THROW(BadInput,
                     "NCMatrix number of rows and columns must both be positive or both zero");
    m_rowcount = m_colcount = 0;
  } else {
    m_data.setByCopy(data.begin(),data.end());
    m_rowcount = row;
    m_colcount = col;
  }
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
  res.m_data.reserve_hint(m_data.size());
  for ( auto e : m_data )
    res.m_data.push_back( e * factor );
  return res;
}

inline NCrystal::Matrix NCrystal::Matrix::operator/ (double factor)
{
  return *this * ( 1.0 / factor );
}

inline bool NCrystal::Matrix::isEqual(const NCrystal::Matrix& o) const
{
  if ( m_rowcount != o.m_rowcount )
    return false;
  if ( m_colcount != o.m_colcount )
    return false;
  nc_assert( m_data.size() == o.m_data.size() );
  auto it = m_data.begin();
  auto it2 = o.m_data.begin();
  auto itE = it + m_data.size();
  nc_assert(itE == m_data.end());
  for ( ; it!=itE; ++it, ++it2 ) {
    if ( *it != *it2 )
      return false;
  }
  return true;
}

inline NCrystal::Matrix NCrystal::Matrix::operator~() const
{
  Matrix res;
  res.m_colcount = m_colcount;
  res.m_rowcount = m_rowcount;
  if (m_data.empty())
    return res;
  res.m_data.reserve_hint(m_data.size());
  for (unsigned i = 0; i < m_colcount; i++)
    for (unsigned j = 0; j < m_rowcount; j++)
      res.m_data.push_back(m_data[j * m_colcount + i]);
  return res;
}

inline void NCrystal::Matrix::transpose()
{
  //TODO: Should be able to do this without malloc and/or extra copies of all
  //elements state..
  decltype(m_data) v;
  std::swap(m_data, v);
  m_data.reserve_hint(v.size());
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
  for ( auto e : m_data )
    sum += e*e;
  return sum;
}

inline NCrystal::Matrix& NCrystal::Matrix::operator*=(const double& f)
{
  for (auto& e : m_data)
    e *= f;
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

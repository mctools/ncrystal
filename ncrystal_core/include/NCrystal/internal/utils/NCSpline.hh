#ifndef NCrystal_Spline_hh
#define NCrystal_Spline_hh

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

  class CubicSpline {
  public:
    CubicSpline();//default constructed instance is invalid until ::set() is called.
    CubicSpline( const VectD& y,
                 double derivative_y_left = 0.0,
                 double derivative_y_right = 0.0 );
    ~CubicSpline();
    void set( const VectD& y,
              double derivative_y_left = 0.0,
              double derivative_y_right = 0.0 );
    double evalWithAssert(double x) const;
    double evalUnbounded(double x) const;
    void swap(CubicSpline&o);
  private:
    friend class SplinedLookupTable;
    std::size_t m_nm2;
    std::vector<PairDD > m_data;
  };

  class SplinedLookupTable {
  public:
    SplinedLookupTable();//default constructed instance is invalid until ::set() is called.

    //Setup splined lookuptable. The parameters name and description are
    //optional, and will only be used if the environment variable
    //NCRYSTAL_DEBUG_SPLINES is set, in which case debug files containing info
    //about the spline will be created (supposedly for later inspection by
    //NCrystal developers). Avoid spaces and special characters in name.

    SplinedLookupTable( const VectD& fvals,double a,double b,double fprime_a, double fprime_b,
                        const std::string& name="", const std::string& description="" );
    SplinedLookupTable( const Fct1D* thefct,double a,double b,double fprime_a, double fprime_b,unsigned npts = 1000,
                        const std::string& name="", const std::string& description="" );
    ~SplinedLookupTable();
    void set( const VectD& fvals,double a,double b,double fprime_a, double fprime_b,
              const std::string& name="", const std::string& description="" );
    void set( const Fct1D* thefct,double a,double b,double fprime_a, double fprime_b,unsigned npts = 1000,
              const std::string& name="", const std::string& description="" );
    double eval(double x) const;//<-- query the resulting lookup table
    void swap(SplinedLookupTable&o);
    double getLower() const { return m_a; }
    double getUpper() const { return m_b; }
    double getInvDelta() const { return m_invdelta; }
  private:
    double m_a;
    double m_invdelta;
    CubicSpline m_spline;
    double m_b;//only used in getUpper
    void producefile( const Fct1D* thefct,
                      double fprime_a, double fprime_b,
                      const std::string& name,const std::string& descr ) const;

  };

  class PiecewiseLinearFct1D {
    //A simple piecewise linear function, with controllabe behaviour for
    //over/underflow requests (an exception or a fixed value).
  public:
    struct OutOfBoundsYValues { Optional<double> underflowYValue, overflowYValue; };

    //Construct by {xi,yi} pairs, where the {xi} values form a valid grid.
    PiecewiseLinearFct1D( VectD x, VectD y, OutOfBoundsYValues = { NullOpt, NullOpt } );
    double eval( double x ) const;
    double operator()( double x ) const { return eval(x); }

    //Dump points to text file for debugging purposes:
    void dumpToFile( const std::string& filename ) const;

    const VectD& xValues() const noexcept { return m_x; }
    const VectD& yValues() const noexcept { return m_y; }
    const OutOfBoundsYValues& outOfBoundsYValues() const noexcept { return m_ofVals; }
    std::size_t size() const { return m_x.size(); }

  private:
    double evalEdgeCase( VectD::const_iterator, double x ) const;
    VectD m_x, m_y;
    OutOfBoundsYValues m_ofVals;
  };


}


////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::CubicSpline::CubicSpline() : m_nm2(0) {}
inline NCrystal::CubicSpline::CubicSpline( const VectD& y, double ypa, double ypb ) { set(y,ypa,ypb); }
inline NCrystal::CubicSpline::~CubicSpline(){}
inline double NCrystal::CubicSpline::evalWithAssert(double x) const {
  nc_assert(x>=0.0&&x<=m_nm2+1);
  return evalUnbounded(x);
}
inline double NCrystal::CubicSpline::evalUnbounded(double x) const {
  nc_assert(m_nm2>0);//will fail if default constructed and set() was never called
  std::size_t idx = ncmin(static_cast<std::size_t>(x),m_nm2);
  double b = x-idx;//fraction inside bin
  double a = 1.0-b;
  nc_assert(idx+1<m_data.size());
  std::vector<PairDD >::const_iterator it = m_data.begin()+idx;
  double tmp = a * (it->first);
  double tmp2 = (a*a*a-a) * (it->second);
  ++it;
  tmp += b*(it->first);
  tmp2 += (b*b*b-b) * (it->second);
  return tmp + 0.166666666666666666666666666666666666666666666666666667 * tmp2;
}


inline NCrystal::SplinedLookupTable::SplinedLookupTable() : m_a(0), m_invdelta(0), m_b(0) {}
inline NCrystal::SplinedLookupTable::SplinedLookupTable( const VectD& fvals,
                                                         double a,double b,
                                                         double fprime_a, double fprime_b,
                                                         const std::string& name,
                                                         const std::string& desc )
{
  set(fvals,a,b,fprime_a,fprime_b,name,desc);
}
inline NCrystal::SplinedLookupTable::~SplinedLookupTable(){}
inline NCrystal::SplinedLookupTable::SplinedLookupTable( const Fct1D* thefct,
                                                         double a,double b,
                                                         double fprime_a, double fprime_b,
                                                         unsigned npts,
                                                         const std::string& name,
                                                         const std::string& desc )
{
  set(thefct,a,b,fprime_a,fprime_b,npts,name,desc);
}

inline double NCrystal::SplinedLookupTable::eval(double x) const {
  return m_spline.evalUnbounded((x-m_a)*m_invdelta);
}

inline void NCrystal::SplinedLookupTable::swap(NCrystal::SplinedLookupTable&o) {
  std::swap(m_a,o.m_a);
  std::swap(m_b,o.m_b);
  std::swap(m_invdelta,o.m_invdelta);
  m_spline.swap(o.m_spline);
}

inline void NCrystal::CubicSpline::swap(NCrystal::CubicSpline&o) {
  std::swap(m_nm2,o.m_nm2);
  std::swap(m_data,o.m_data);
}

inline NCrystal::PiecewiseLinearFct1D::PiecewiseLinearFct1D( VectD x,
                                                             VectD y,
                                                             OutOfBoundsYValues ofVals )
  : m_x(std::move(x)),
    m_y(std::move(y)),
    m_ofVals(std::move(ofVals))
{
  nc_assert_always( m_x.size() >= 2 );
  nc_assert_always( m_x.size() == m_y.size() );
  nc_assert( nc_is_grid(m_x) );
#ifndef NDEBUG
  for ( auto& e : m_y ) {
    nc_assert( !ncisinf(e) );
    nc_assert( !ncisnan(e) );
  }
#endif
  m_x.shrink_to_fit();
  m_y.shrink_to_fit();
}

inline double NCrystal::PiecewiseLinearFct1D::eval( double x ) const {
  nc_assert(!ncisnan(x));
  auto it = std::lower_bound( m_x.begin(), m_x.end(), x );
  if ( it == m_x.end() || it == m_x.begin() )
    return evalEdgeCase( it, x );
  auto idx = std::distance(m_x.begin(),it);
  auto itY = std::next(m_y.begin(),idx);
  const double x1 = *it--;
  const double x0 = *it;
  const double y1 = *itY--;
  const double y0 = *itY;
  return y0 + (y1-y0 ) * ( x - x0 ) / ( x1 - x0 );
}

#endif

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

#include "NCrystal/internal/utils/NCSpline.hh"
#include "NCrystal/internal/utils/NCFastSearch.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCFileUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include <fstream>
#include <iomanip>
#include <sstream>

namespace NC = NCrystal;

//TODO: For convenience, CubicSpline derivatives should default to
//special values (e.g. inf), in which case derivatives will be approximated,
//either linearly from the edge points or using a numerical derivative
//estimation.

void NC::CubicSpline::set( const VectD& y,
                           double derivative_y_left,
                           double derivative_y_right )
{
  const std::size_t n = y.size();
  nc_assert_always(n>3);
  const std::size_t nm1 = n-1;
  VectD u(nm1,0.0);
  VectD y2(n,0.0);
  y2[0] = -0.5;
  double dy_left = y[1]-y[0];
  u[0] = 3.0*( dy_left-derivative_y_left );
  for (std::size_t i = 1; i<nm1;++i) {
    double p = 2.0 + 0.5 * y2[i-1];
    nc_assert(p!=0.0);
    p = 1.0/p;
    y2[i] = -0.5*p;
    double uu = y[i+1] - 2.0*y[i] + y[i-1];
    u[i] = (3.0*uu -0.5*u[i-1])*p;
  }
  double qn = 0.5;
  double un = 3.0 * (derivative_y_right -(y[nm1]-y[n-2]));
  nc_assert((qn*y2[n-2]+1.0) != 0.0);
  y2[nm1] = (un-qn*u[n-2]) / (qn*y2[n-2]+1.0);
  for (std::size_t k = nm1; k; --k) {
    std::size_t km1 = k-1;
    y2[km1] *= y2[k];
    y2[km1] += u[km1];
  }
  std::vector<PairDD > data;
  data.reserve(y.size());
  for (std::size_t i = 0; i < y.size(); ++i) {
    data.push_back(std::make_pair(y[i],y2[i]));
    nc_assert(!ncisnan(data.back().first));
    nc_assert(!ncisnan(data.back().second));
  }
  //all good, set:
  std::swap(m_data,data);
  m_nm2 = n-2;
}

void NC::SplinedLookupTable::set( const VectD& fvals,
                                  double a,double b,
                                  double fprime_a, double fprime_b,
                                  const std::string& name,
                                  const std::string& description  )
{
  nc_assert(b>a);
  nc_assert(fvals.size()>3);
  m_a = a;
  m_b = b;
  double delta = (b-a)/(fvals.size()-1);
  m_spline.set( fvals, fprime_a*delta, fprime_b*delta );
  nc_assert(delta>0.0);
  m_invdelta = 1.0 / delta;
  if (ncgetenv_bool("DEBUG_SPLINES"))
    producefile( 0/*true function not known*/,fprime_a, fprime_b,name,description );

}

void NC::SplinedLookupTable::set( const Fct1D* thefct,
                                  double a,double b,
                                  double fprime_a, double fprime_b,
                                  unsigned npts,
                                  const std::string& name,
                                  const std::string& description )
{
  nc_assert(!ncisnan(fprime_a));
  nc_assert(!ncisnan(fprime_b));
  nc_assert(thefct);
  nc_assert(b>a);
  nc_assert(npts>3);
  m_a = a;
  m_b = b;
  VectD fvals;
  fvals.reserve(npts);
  double delta = (b-a)/(npts-1);
  std::size_t nm1 = npts-1;
  for (unsigned i = 0; i < nm1; ++i) {
    fvals.push_back(thefct->eval(a+i*delta));
    nc_assert(!ncisnan(fvals.back()));
  }
  fvals.push_back(thefct->eval(b));
  nc_assert(!ncisnan(fvals.back()));
  m_spline.set( fvals, fprime_a*delta, fprime_b*delta );
  nc_assert(delta>0.0);
  m_invdelta = 1.0 / delta;


  if (ncgetenv_bool("DEBUG_SPLINES"))
    producefile( thefct,fprime_a, fprime_b,name,description );
}

void NC::SplinedLookupTable::producefile( const Fct1D* thefct,
                                          double fprime_a, double fprime_b,
                                          const std::string& username,
                                          const std::string& userdesc ) const
{

  std::string name = (username.empty()?std::string("unknownspline"):username);
  std::string description = (userdesc.empty()?std::string("no description"):userdesc);
  std::string filename;
  for (unsigned ifile = 0;ifile<=99;++ifile) {
    std::stringstream s;
    s << "ncrystal_splinedlt_"<<name;
    if (ifile>0)
      s<<"_"<<ifile;
    s<<".txt";
    filename = s.str();
    if (!file_exists(filename))
      break;
  }

  if (file_exists(filename)) {
    NCRYSTAL_WARN("Aborted writing of spline data file ncrystal_splinedlt_"
                  <<username<<"_??.txt - too many files in current dir.");
    return;
  }

  std::ofstream ofs (filename.c_str(), std::ofstream::out);
  ofs << std::setprecision(20);
  ofs << "#ncrystal_splined_lookuptable\n";
  ofs << "#name = "<<name<<"\n";
  ofs << "#description = "<<description<<"\n";
  ofs << "#a = "<<getLower()<<"\n";
  ofs << "#b = "<<getUpper()<<"\n";
  ofs << "#fprime_a = "<<fprime_a<<"\n";
  ofs << "#fprime_b = "<<fprime_b<<"\n";
  ofs << "#input_fvals = ";
  for (std::size_t i = 0; i < m_spline.m_data.size(); ++i)
    ofs<<" "<<m_spline.m_data[i].first;
  ofs << "\n#data_colums = x,spline_of_x";
  if (thefct)
    ofs <<",truefct_of_x";
  ofs<<"\n";
  std::size_t numpts = 100*m_spline.m_data.size();
  if (numpts>1000000)
    numpts = std::max<std::size_t>(numpts/10,1000000);
  double delta = (m_b-m_a)/(numpts-1.0);
  for (std::size_t i = 0; i<=numpts;++i) {
    double x = (i==numpts?m_b:m_a+i*delta);
    ofs << x << " "<<eval(x);
    if (thefct)
      ofs<<" "<<thefct->eval(x);
    ofs<<"\n";
  }
  NCRYSTAL_MSG("Wrote "<<filename<<" (since NCRYSTAL_DEBUG_SPLINE is set).");
}

void NC::PiecewiseLinearFct1D::dumpToFile( const std::string& filename ) const {
  std::ofstream ofs (filename.c_str(), std::ofstream::out);
  ofs << std::setprecision(20);
  ofs << "#colnames=x,y\n";
  ofs << "#plotstyle=*-\n";
  ofs << "#overflow=";
  if  ( m_ofVals.overflowYValue.has_value() )
    ofs << "none\n";
  else
    ofs << m_ofVals.overflowYValue.value()<< "\n";
  ofs << "#underflow=";
  if  ( m_ofVals.underflowYValue.has_value() )
    ofs << "none\n";
  else
    ofs << m_ofVals.underflowYValue.value()<< "\n";

  for ( auto i : ncrange(m_x.size()) ) {
    ofs << m_x.at(i) << " "<<m_y.at(i) << "\n";
  }
  ofs << std::flush;
}

double NC::PiecewiseLinearFct1D::evalEdgeCase( VectD::const_iterator it,
                                               double x ) const
{
  if ( it == m_x.end() ) {
    if ( !m_ofVals.overflowYValue.has_value() )
      NCRYSTAL_THROW2(CalcError,"PiecewiseLinearFct1D: Out of bounds: x>xmax"
                      " and no overflow value supplied (x="<<x<<", xmax="<<m_x.back()<<").");
    return m_ofVals.overflowYValue.value();
  };
  if ( x >= m_x.front() )
    return m_y.front();
  if ( !m_ofVals.underflowYValue.has_value() )
    NCRYSTAL_THROW2(CalcError,"PiecewiseLinearFct1D: Out of bounds: x<xmin"
                    " and no underflow value supplied (x="<<x<<", xmin="<<m_x.front()<<").");
  return m_ofVals.underflowYValue.value();
}

NC::PCHIPInterp::PCHIPInterp( Span<const double> x, Span<const double> y,
                              Optional<double> opt_ls, Optional<double> opt_rs )
{
#ifndef NDEBUG
  nc_assert(x.size() == y.size());
  nc_assert(nc_is_grid(x));
  nc_assert(std::isfinite(opt_ls.value_or(1.0)));
  nc_assert(std::isfinite(opt_rs.value_or(1.0)));
  for ( auto e : y )
    nc_assert( std::isfinite(e) );
#endif
  m_data.resize(5u * x.size() - 2u);
  m_n = x.size();

  const std::size_t oh = 2 * m_n;
  const std::size_t oi = 3 * m_n - 1;
  const std::size_t om = 4 * m_n - 2;

  for (std::size_t i = 0; i < m_n; ++i) {
    vectAt(m_data, i) = x[i];
    vectAt(m_data, m_n + i) = y[i];
  }

  for (std::size_t i = 0; i + 1 < m_n; ++i) {
    const double h = x[i + 1] - x[i];
    nc_assert( std::isfinite(h) );
    vectAt(m_data, oh + i) = h;
    vectAt(m_data, oi + i) = 1.0 / h;
  }

  auto sec = [this,oi](std::size_t i) {
    return ( ( vectAt(m_data, m_n + i + 1) - vectAt(m_data, m_n + i) )
             * vectAt(m_data, oi + i) );
  };

  auto limit = [](double a, double b) {
    if (b == 0.0 || a * b <= 0.0)
      return 0.0;
    if (ncabs(a) > 3.0 * ncabs(b))
      return 3.0 * b;
    return a;
  };

auto edgeSlope = [this,&sec](std::size_t i0,
                               std::size_t i1) {
    if (m_n == 2)
      return sec(0);

    const double h0 = vectAt(m_data, 2 * m_n + i0);
    const double h1 = vectAt(m_data, 2 * m_n + i1);
    const double d0 = sec(i0);
    const double d1 = sec(i1);

    double a = ((2.0 * h0 + h1) * d0 - h0 * d1)
             / (h0 + h1);

    if (a * d0 <= 0.0)
      a = 0.0;
    else if (d0 * d1 < 0.0 && ncabs(a) > 3.0 * ncabs(d0))
      a = 3.0 * d0;

    return a;
  };

  const double ls = ( opt_ls.has_value() ?
                      opt_ls.value() :
                      edgeSlope(0, 1) );
  const double rs = ( opt_rs.has_value() ?
                      opt_rs.value() :
                      edgeSlope(m_n - 2, m_n - 3) );

  vectAt(m_data, om) = limit(ls, sec(0));

  for (std::size_t i = 1; i + 1 < m_n; ++i) {
    const double a = sec(i - 1);
    const double b = sec(i);
    if (a * b <= 0.0) {
      vectAt(m_data, om + i) = 0.0;
    } else {
      const double w0
        = 2.0 * vectAt(m_data, oh + i) + vectAt(m_data, oh + i - 1);
      const double w1
        = vectAt(m_data, oh + i) + 2.0 * vectAt(m_data, oh + i - 1);
      vectAt(m_data, om + i) = (w0 + w1) / (w0 / a + w1 / b);
    }
  }

  vectAt(m_data, om + m_n - 1) = limit(rs, sec(m_n - 2));
}

double NC::PCHIPInterp::eval(double x) const
{
#ifndef NDEBUG
  nc_assert_always( m_n >= 2 );
  nc_assert(std::isfinite(x));
#endif
  const double* p = m_data.data();
  const double* xp = p;
  const double* yp = p + m_n;
  const double* hp = p + 2 * m_n;
  const double* ip = p + 3 * m_n - 1;
  const double* mp = p + 4 * m_n - 2;
#ifndef NDEBUG
  nc_assert(x >= xp[0]);
  nc_assert(x <= xp[m_n - 1]);
#endif
  if (x == xp[m_n - 1])
    return yp[m_n - 1];
  const std::size_t i = NC::fastUpperBoundIdx(xp, m_n, x) - 1;
  const double t = (x - xp[i]) * ip[i];
  const double t2 = t * t;
  const double u = 1.0 - t;
  const double u2 = u * u;
  const double h00 = (1.0 + 2.0 * t) * u2;
  const double h10 = t * u2;
  const double h01 = t2 * (3.0 - 2.0 * t);
  const double h11 = t2 * (t - 1.0);
  return ( h00*yp[i] + h10*hp[i]*mp[i] + h01*yp[i + 1] + h11*hp[i]*mp[i + 1] );
}

NC::PairDD NC::PCHIPInterp::endPointSlopes() const
{
  nc_assert_always(m_n >= 2);
  return { vectAt(m_data, 4 * m_n - 2), vectAt(m_data, 5 * m_n - 3 ) };
}

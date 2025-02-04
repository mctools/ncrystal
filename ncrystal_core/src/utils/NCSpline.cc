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

#include "NCrystal/internal/utils/NCSpline.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <fstream>
#include <iomanip>

//TODO: For convenience, CubicSpline derivatives should default to
//special values (e.g. inf), in which case derivatives will be approximated,
//either linearly from the edge points or using a numerical derivative
//estimation.

void NCrystal::CubicSpline::set( const VectD& y,
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

void NCrystal::SplinedLookupTable::set( const VectD& fvals,
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

void NCrystal::SplinedLookupTable::set( const Fct1D* thefct,
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


#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCFileUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include <fstream>
#include <sstream>
#include <iomanip>

void NCrystal::SplinedLookupTable::producefile( const Fct1D* thefct,
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

void NCrystal::PiecewiseLinearFct1D::dumpToFile( const std::string& filename ) const {
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

double NCrystal::PiecewiseLinearFct1D::evalEdgeCase( VectD::const_iterator it, double x ) const
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

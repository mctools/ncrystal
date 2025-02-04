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

#include "NCrystal/internal/utils/NCPointwiseDist.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include <cstdio>

namespace NC = NCrystal;

NC::PointwiseDist::PointwiseDist(const VectD &xvals, const VectD &yvals)
  : PointwiseDist( VectD(xvals), VectD(yvals) )
{
}

NC::PointwiseDist::PointwiseDist( VectD&& xvals, VectD&& yvals )
  : m_x(xvals), m_y(yvals)
{
  if(m_x.size()!=m_y.size() || m_y.size()<2 )
    NCRYSTAL_THROW(CalcError, "input vector size error.");

  if(!std::is_sorted(m_x.begin(),m_x.end()))
    NCRYSTAL_THROW(CalcError, "points of the distribution are not sorted.");

  for(std::size_t i=0;i<m_y.size();i++)
  {
    if(m_y[i]<0.)
      NCRYSTAL_THROW(CalcError, "function is negative.");

    if(ncisinf(m_y[i]))
      NCRYSTAL_THROW(CalcError, "function is not finite.");
  }

  m_cdf.reserve(m_y.size());
  StableSum totalArea;

  m_cdf.push_back(0.);
  for(std::size_t i=1;i<m_y.size();i++)
  {
    double area = (m_x[i]-m_x[i-1])*0.5*(m_y[i]+m_y[i-1]);
    if(area<0)
      NCRYSTAL_THROW(CalcError, "Negative probability density");
    totalArea.add( area );
    m_cdf.push_back( totalArea.sum() );
  }

  double totalAreaVal = totalArea.sum();
  if ( !(totalAreaVal>0.0) )
    NCRYSTAL_THROW(CalcError, "No area in distribution.");

  double normfact = 1.0/totalAreaVal;
  for ( auto& e : m_cdf )
    e *= normfact;
  for ( auto& e : m_y )
    e *= normfact;
  nc_assert( ncabs(1.0-m_cdf.back()) < 1.0e-14 );
  m_cdf.back() = 1.0;
}

std::pair<double,unsigned> NC::PointwiseDist::percentileWithIndex(double p ) const
{
  nc_assert(p>=0.&&p<=1.0);
  nc_assert(m_x.size()<std::numeric_limits<unsigned>::max());
  if(p==1.)
    return std::pair<double,unsigned>(m_x.back(),
                                      static_cast<unsigned>(m_x.size()-2));

  std::size_t i = std::max<std::size_t>(std::min<std::size_t>(std::lower_bound(m_cdf.begin(), m_cdf.end(), p)-m_cdf.begin(),m_cdf.size()-1),1);
  nc_assert( i>0 && i < m_x.size() );
  double dx = m_x[i]-m_x[i-1];
  double c = (p-m_cdf[i-1]);
  double a = m_y[i-1];
  double d = m_y[i] - a;
  double zdx;
  if (!a) {
    zdx = d>0.0 ? std::sqrt( ( 2.0 * c * dx ) / d ) : 0.5*dx;//(a,d)=(0,0) should not really happen?
  } else {
    double e = d * c / ( dx * a * a );
    if (ncabs(e)>1e-7) {
      //apply formula:
      zdx = ( std::sqrt( 1.0 + 2.0 * e ) - 1.0 ) * dx * a / d;
    } else {
      //calculate via expansion (solves numerical issue when d is near zero):
      zdx = ( 1 + 0.5 * e * ( e - 1.0 ) ) * c / a;
    }
  }
  return std::pair<double,unsigned>( ncclamp(m_x[i-1] + zdx,m_x[i-1],m_x[i]),
                                     static_cast<unsigned>( i-1 ) );
}

double NC::PointwiseDist::commulIntegral( double x ) const
{
  //Above or below edges is easy:
  if ( x <= m_x.front() )
    return 0.0;
  if ( x >= m_x.back() )
    return 1.0;

  //Find bin with binary search:
  auto it = std::upper_bound( m_x.begin(), m_x.end(), x );
  nc_assert( it != m_x.end() );
  nc_assert( it != m_x.begin() );

  //We are in the interval [std::prev(it),it], find parameters of this last bin:
  auto i1 = std::distance(m_x.begin(),it);
  nc_assert(i1>0);
  auto i0 = i1 - 1;
  const double x1 = vectAt(m_x,i0);
  const double y1 = vectAt(m_y,i0);
  const double x2 = vectAt(m_x,i1);
  const double y2 = vectAt(m_y,i1);

  //Find contribution in this bin as as
  //<length in bin>*<average height in bin over used part>:
  nc_assert( x2 - x1 > 0.0 );
  const double dx = x-x1;
  const double slope = ( y2-y1 ) / (x2-x1);
  const double last_bin_contrib = dx * ( y1 + 0.5 * dx * slope );

  //Combine with preceding bins from m_cdf:
  return vectAt(m_cdf,i0) + last_bin_contrib;
}

double NC::PointwiseDist::sampleBelow( RNG& rng, double xtrunc ) const
{
  //Above or below edges is easy:
  if ( xtrunc <= m_x.front() ) {
    if ( xtrunc == m_x.front() )
      return m_x.front();
    NCRYSTAL_THROW2(BadInput,"PointwiseDist::sampleBelow asked to sample point below distribution");
  }
  if ( xtrunc >= m_x.back() )
    return sample(rng);

  return percentile( rng.generate() * commulIntegral( xtrunc ) );
}

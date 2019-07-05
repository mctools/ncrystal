////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2019 NCrystal developers                                   //
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

#include "NCPointwiseDist.hh"
#include "NCMath.hh"
#include <algorithm>
#include <cstdio>

NCrystal::PointwiseDist::PointwiseDist(const std::vector<double> &xvals, const std::vector<double> &yvals, double iw)
  : m_x(xvals), m_y(yvals), m_iweight(iw)
{
  if(m_x.size()!=m_y.size() || m_y.size()<2 )
    NCRYSTAL_THROW(CalcError, "input vector size error.");

  if(!ncis_sorted(m_x.begin(),m_x.end()))
    NCRYSTAL_THROW(CalcError, "points of the distribution are not sorted.");

  for(std::size_t i=0;i<m_y.size();i++)
  {
    if(m_y[i]<0.)
      NCRYSTAL_THROW(CalcError, "function is negative.");

    if(ncisinf(m_y[i]))
      NCRYSTAL_THROW(CalcError, "function is not finite.");
  }

  m_cdf.reserve(m_y.size());
  double totalArea = 0.;

  m_cdf.push_back(0.);
  for(std::size_t i=1;i<m_y.size();i++)
  {
    double area = (m_x[i]-m_x[i-1])*0.5*(m_y[i]+m_y[i-1]);
    if(area<0)
      NCRYSTAL_THROW(CalcError, "Negative probability density");
    totalArea+=area;
    m_cdf.push_back(totalArea);
  }

  if (!totalArea)
    NCRYSTAL_THROW(CalcError, "No area in distribution.");

  double normfact = 1.0/totalArea;
  for(std::size_t i=0;i<m_cdf.size();i++)
  {
    m_cdf[i] *= normfact;
    m_y[i] *= normfact;
  }
}

NCrystal::PointwiseDist::~PointwiseDist()
{
}

double NCrystal::PointwiseDist::percentile(double p ) const
{
  nc_assert(p>=0.&&p<=1.0);
  if(p==1.)
    return m_x.back();

  std::size_t i = std::max<std::size_t>(std::min<std::size_t>(std::lower_bound(m_cdf.begin(), m_cdf.end(), p)-m_cdf.begin(),m_cdf.size()-1),1);
  double dx = m_x[i]-m_x[i-1];
  double c = (p-m_cdf[i-1]);
  double a = m_y[i-1];
  double d = m_y[i] - a;
  double zdx;
  if (!a) {
    zdx =  d>0.0 ? std::sqrt( ( 2.0 * c * dx ) / d ) : 0.5*dx;//a=0 and d=0 should not really happen...
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
  return ncclamp(m_x[i-1] + zdx,m_x[i-1],m_x[i]);
}

void NCrystal::PointwiseDist::setIntegralWeight(double iw)
{
  m_iweight=iw;
}


NCrystal::PointwiseDist& NCrystal::PointwiseDist::operator+=(const NCrystal::PointwiseDist& right)
{
  if(right.m_x.size()!=this->m_x.size() || right.m_x.size()!=this->m_y.size() || right.m_x.size()!=this->m_cdf.size())
    NCRYSTAL_THROW(CalcError,"PointwiseDist objects are not compatible (grid-sizes differs).");
  for(unsigned i=0;i<right.m_x.size();++i)
  {
    if(this->m_x[i]!=right.m_x[i])
      NCRYSTAL_THROW(CalcError,"Can not add distributions with different grid values.");
  }
  
  double totweight =  this->m_iweight + right.m_iweight;
  double ratiothis = this->m_iweight /totweight;
  double ratioright = right.m_iweight /totweight;

  for(unsigned i=0;i<right.m_x.size();i++)
  {
    this->m_y[i] = this->m_y[i]*ratiothis + right.m_y[i]*ratioright;
    this->m_cdf[i] = this->m_cdf[i]*ratiothis + right.m_cdf[i]*ratioright;
  }

  this->m_iweight = totweight;

  return *this;
}

NCrystal::PointwiseDist& NCrystal::PointwiseDist::operator*= (double frac)
{
  this->m_iweight *= frac;
  return *this;
}


void NCrystal::PointwiseDist::print() const
{
  for(unsigned i=0;i<m_cdf.size();i++)  {
    printf("cdf idx %u, val %e\n",i, m_cdf[i]);
  }
}

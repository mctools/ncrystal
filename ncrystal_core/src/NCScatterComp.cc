////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2018 NCrystal developers                                   //
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

#include "NCrystal/NCScatterComp.hh"
#include "NCrystal/NCDefs.hh"
#include "NCRandUtils.hh"
#include <algorithm>
#include <limits>

NCrystal::ScatterComp::ScatterComp(const char * calculator_type_name)
  : Scatter(calculator_type_name), m_threshold_lower(0.0), m_threshold_upper(infinity), m_isOriented(-1)
{
}

NCrystal::ScatterComp::~ScatterComp()
{
  std::vector<Component>::const_iterator it = m_calcs.begin();
  std::vector<Component>::const_iterator itE = m_calcs.end();
  for (;it!=itE;++it)
    it->scatter->unref();
}

bool NCrystal::ScatterComp::Component::operator<(const NCrystal::ScatterComp::Component& o) const
{
  return o.threshold_lower > threshold_lower;
}

void NCrystal::ScatterComp::addComponent(Scatter* scat, double thescale )
{
  RCGuard guard(scat);//ensure we always ref/unref scat even in case of exceptions.
  if (!scat)
    NCRYSTAL_THROW(BadInput,"ScatterComp::addComponent Got NULL scatter.");
  if (thescale<0.0)
    NCRYSTAL_THROW(BadInput,"ScatterComp::addComponent Component scale is negative.");
  std::vector<Component>::const_iterator it(m_calcs.begin()), itE(m_calcs.end());
  for (;it!=itE;++it) {
    if (it->scatter == scat)
        NCRYSTAL_THROW(BadInput,"ScatterComp::addComponent got same scatter multiple times.");
  }
  m_calcs.reserve(m_calcs.size()+1);
  scat->validate();
  Component c;
  c.scale = thescale;
  c.scatter = scat;
  scat->domain(c.threshold_lower,c.threshold_upper);
  if (m_calcs.empty() || c.threshold_lower < m_threshold_lower)
    m_threshold_lower = c.threshold_lower;
  if (m_calcs.empty() || c.threshold_upper > m_threshold_upper)
    m_threshold_upper = c.threshold_upper;
  registerSubCalc(scat);
  m_calcs.push_back(c);
  scat->ref();

  //Sort by threshold so lowers comes first, preserving the original order when
  //thresholds are equal:
  std::stable_sort(m_calcs.begin(),m_calcs.end());

  m_isOriented = -1;//invalidate (scatter might be incomplete ScatterComp so we
                    //can't always know already if it is oriented or not).

  validate();
}

double NCrystal::ScatterComp::crossSection(double ekin, const double (&indir)[3] ) const
{
  double c(0);
  std::vector<Component>::const_iterator it = m_calcs.begin();
  std::vector<Component>::const_iterator itE = m_calcs.end();
  if (it==itE)
    NCRYSTAL_THROW(BadInput,"ScatterComp::crossSection queried with no components added.");
  for (;it!=itE;++it) {
    if (ekin<it->threshold_lower)
      break;
    if (ekin>it->threshold_upper)
      continue;
    c += it->scatter->crossSection(ekin,indir) * it->scale;
  }
  return c;
}

void NCrystal::ScatterComp::generateScattering( double ekin, const double (&indir)[3],
                                                double (&outdir)[3], double& de ) const
{
  double rand_choice = getRNG()->generate() * crossSection(ekin,indir);
  double c(0);
  std::vector<Component>::const_iterator it = m_calcs.begin();
  std::vector<Component>::const_iterator itE = m_calcs.end();
  if (it==itE)
    NCRYSTAL_THROW(BadInput,"ScatterComp::generateScattering queried with no components added.");
  for (;it!=itE;++it) {
    if (ekin<it->threshold_lower)
      break;
    if (ekin>it->threshold_upper)
      continue;
    c += it->scatter->crossSection(ekin,indir) * it->scale;
    if (rand_choice <= c) {
      it->scatter->generateScattering(ekin, indir, outdir, de);
      return;
    }
  }
  //Should get here only in case of rounding errors or if called outside
  //domain(). No cross-section means no action:
  outdir[0] = indir[0];
  outdir[1] = indir[1];
  outdir[2] = indir[2];
  de = 0;
}

bool NCrystal::ScatterComp::isOriented() const {
  if (m_isOriented==-1)
    checkIsOriented();
  return (bool)m_isOriented;
}

void NCrystal::ScatterComp::checkIsOriented() const
{
  m_isOriented = 0;
  std::vector<Component>::const_iterator it = m_calcs.begin();
  std::vector<Component>::const_iterator itE = m_calcs.end();
  for (;it!=itE;++it) {
    if (it->scatter->isOriented()) {
      m_isOriented = 1;
      break;
    }
  }
}


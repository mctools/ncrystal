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

#include "NCLCRefModels.hh"
#include "NCRandUtils.hh"

namespace NC = NCrystal;

NC::LCBraggRef::LCBraggRef(Scatter* scb, Vector lcaxis_lab, unsigned nsample)
  : Scatter("LCBraggRef"),
    m_sc(scb),
    m_lcaxislab(lcaxis_lab.unit()),
    m_nsample(nsample),
    m_nsampleprime(nsample)
{
  registerSubCalc(scb);
  while (!isPrime(m_nsampleprime))
    ++m_nsampleprime;
}

NC::LCBraggRef::~LCBraggRef()
{
}

void NC::LCBraggRef::domain(double& ekin_low, double& ekin_high) const
{
  return m_sc.obj()->domain(ekin_low,ekin_high);
}

double NC::LCBraggRef::crossSection( double ekin, const double (&indirraw)[3] ) const
{
  Vector indir = asVect(indirraw).unit();
  Vector lccross = m_lcaxislab.cross(indir);
  double lcdot = m_lcaxislab.dot(indir);
  double sumxs = 0.0;
  for (unsigned i = 0; i<m_nsampleprime; ++i) {
    PhiRot phirot(i*((2.0*M_PI)/m_nsampleprime)-M_PI);
    Vector v = phirot.rotateVectorAroundAxis( indir, m_lcaxislab, lccross, lcdot );
    sumxs += m_sc.obj()->crossSection(ekin,NC_CVECTOR_CAST(v));
  }
  return sumxs/m_nsampleprime;
}

void NC::LCBraggRef::generateScattering( double ekin,
                                         const double (&indirraw)[3],
                                         double (&outdir)[3],
                                         double& delta_ekin ) const
{
  Vector indir = asVect(indirraw).unit();
  Vector lccross = m_lcaxislab.cross(indir);
  double lcdot = m_lcaxislab.dot(indir);

  std::vector<double> xs;
  std::vector<PhiRot> pr;
  xs.reserve(m_nsample);
  pr.reserve(m_nsample);

  double sumxs = 0.0;

  //Get cross-sections at nsample random phi rotations:
  RandomBase * rand = getRNG();
  for (unsigned i = 0; i<m_nsample; ++i) {
    double cosphi,sinphi;
    randPointOnUnitCircle( rand, cosphi, sinphi );
#if __cplusplus >= 201103L
    pr.emplace_back(cosphi,sinphi);
#else
    pr.push_back(PhiRot(cosphi,sinphi));
#endif
    Vector v = pr.back().rotateVectorAroundAxis( indir, m_lcaxislab, lccross, lcdot );
    xs.push_back( sumxs += m_sc.obj()->crossSection(ekin,NC_CVECTOR_CAST(v)) );
  }

  if (!sumxs) {
    //no xs, do nothing.
    asVect(outdir)=indir;
    delta_ekin = 0;
    return;
  }
  //Select one phi rotation at random:
  PhiRot& phirot = pr.at(pickRandIdxByWeight(rand,xs));

  //Scatter!
  Vector v = phirot.rotateVectorAroundAxis( indir, m_lcaxislab, lccross, lcdot );
  Vector outdir_rot;
  m_sc.obj()->generateScattering(ekin, NC_CVECTOR_CAST(v),
                                 NC_VECTOR_CAST(outdir_rot), delta_ekin);
  asVect(outdir) = phirot.rotateVectorAroundAxis( outdir_rot, m_lcaxislab, true/*reverse*/);
}

NC::LCBraggRndmRot::LCBraggRndmRot(Scatter* scb, Vector lcaxis_lab, unsigned nsample)
  : Scatter("LCBraggRndmRot"),
    m_sc(scb),
    m_lcaxislab(lcaxis_lab.unit()),
    m_nsample(nsample)
{
  registerSubCalc(scb);
  nc_assert_always(nsample>0);
  cache.rotations.reserve(nsample);
  cache.xscommul.reserve(nsample);
}

NC::LCBraggRndmRot::~LCBraggRndmRot()
{
}

void NC::LCBraggRndmRot::domain(double& ekin_low, double& ekin_high) const
{
  return m_sc.obj()->domain(ekin_low,ekin_high);
}

double NC::LCBraggRndmRot::crossSection( double ekin, const double (&indirraw)[3] ) const
{
  //We always regenerate directions on each cross-section call!
  cache.rotations.clear();
  cache.xscommul.clear();
  Vector indir = asVect(indirraw).unit();
  Vector lccross = m_lcaxislab.cross(indir);
  double lcdot = m_lcaxislab.dot(indir);
  double sumxs = 0.0;
  RandomBase * rand = getRNG();
  for (unsigned i = 0; i<m_nsample; ++i) {
    double cosphi, sinphi;
    randPointOnUnitCircle( rand, cosphi, sinphi );
#if __cplusplus >= 201103L
    cache.rotations.emplace_back(cosphi, sinphi);
#else
    cache.rotations.push_back(PhiRot(cosphi, sinphi));
#endif
    Vector v = cache.rotations.back().rotateVectorAroundAxis( indir, m_lcaxislab, lccross, lcdot );
    cache.xscommul.push_back(sumxs += m_sc.obj()->crossSection(ekin,NC_CVECTOR_CAST(v)));
  }
  return cache.xscommul.back()/m_nsample;
}

void NC::LCBraggRndmRot::generateScattering( double ekin,
                                             const double (&indirraw)[3],
                                             double (&outdir)[3],
                                             double& delta_ekin ) const
{
  delta_ekin = 0;

  if (cache.rotations.empty())
    crossSection(ekin,indirraw);//trigger generation of random directions.
  nc_assert(!cache.xscommul.empty());

  if (!cache.xscommul.back()) {
    //no xs, do nothing.
    asVect(outdir) = asVect(indirraw);
    return;
  }

  //Select one phi rotation at random:
  PhiRot& phirot = cache.rotations.at(pickRandIdxByWeight(getRNG(),cache.xscommul));

  //Scatter!
  nc_assert(asVect(indirraw).isUnitVector());
  Vector v = phirot.rotateVectorAroundAxis( asVect(indirraw), m_lcaxislab);
  Vector outdir_rot;
  m_sc.obj()->generateScattering(ekin, NC_CVECTOR_CAST(v),
                                 NC_VECTOR_CAST(outdir_rot), delta_ekin);
  asVect(outdir) = phirot.rotateVectorAroundAxis( outdir_rot, m_lcaxislab, true/*reverse*/);
}

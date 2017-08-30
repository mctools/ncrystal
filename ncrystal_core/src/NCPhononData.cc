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

#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>
#include "NCPhononData.hh"
#include "NCrystal/NCException.hh"
#include "NCMath.hh"
#include "NCFastConvolve.hh"

//first element in m_phvec is the one-phonon data
NCrystal::PhononNumDyns::PhononNumDyns(const std::vector<double> &spec, double btime, double dt,double kt)
: m_spec(spec.begin(),spec.end())
{
  m_btime = btime;
  m_dt = dt;
  nc_assert(m_dt>0.0);
  m_inv_dt = 1.0/dt;
  m_kt = kt;
  nc_assert(m_spec.size()>2);
  m_spec_size_minus_1 = m_spec.size() - 1;
  m_etime = m_btime+(m_spec.size()-1)*m_dt;

  nc_assert_always(!m_spec.empty());
  std::vector<double>::const_iterator max_it = std::max_element(m_spec.begin(), m_spec.end());
  m_asym_sigma = 1.0 / (*max_it * 2.0*M_PI);

  int center_asym_index = (int)(max_it-m_spec.begin());
  m_asym_centre = ( center_asym_index- (int)(m_spec.size())/2 )*m_dt;

  //asym form is a pdf. The area must be scaled to unity.
  //This integration looks just like the zeroth moment sum rule for S.
  double asym_area = 0.;
  for(unsigned i=0;i<m_spec.size();i++)
  {
    asym_area +=  m_spec[i];
  }

  asym_area *= dt;
  double inv_asym_area(1.0/asym_area);
  for(unsigned i=0;i<m_spec.size();i++)
  {
    m_spec[i] *= inv_asym_area;
  }

}

NCrystal::PhononNumDyns::~PhononNumDyns(){}

double NCrystal::PhononNumDyns::getPhononDensity(double time) const
{
  if(time < m_btime || time > m_etime)
    return 0.;
  double a = (time-m_btime)*m_inv_dt;
  double floor_a = std::floor(a);
  size_t index = floor_a;
  double f = a - floor_a;//a-index instead would mix int and double => slower.
  nc_assert(index<m_spec_size_minus_1);
  return m_spec[index]*(1.0-f) + f*m_spec[index+1];
}

const std::vector<double> &NCrystal::PhononNumDyns::getSpectrum() const
{
  return m_spec;
}



NCrystal::PhononCalculator::PhononCalculator(const std::vector<double> &spec, double btime, double dt, unsigned phnum, double kt)
{
  m_kt=kt;
  m_phvec.reserve(phnum+1);
  m_phvec.push_back(PhononNumDyns(spec, btime, dt, kt));
  m_fftnum=phnum;
  calcfftconv(phnum);
}

NCrystal::PhononCalculator::~PhononCalculator()
{
}


void NCrystal::PhononCalculator::worker_conv(const FastConvolve& fc, unsigned order)
{
  //unsigned idx1 = 0;
  //unsigned idx2 = order-1;

  unsigned idx1 = order/2;
  unsigned idx2 = order - idx1-1;

  if(idx1>=m_phvec.size() || idx2>=m_phvec.size())
  {
    NCRYSTAL_THROW(CalcError,"NCrystal::PhononCalculator::conv input data is not available.");
  }

  double dt = m_phvec.at(idx1).getDeltaTime();
  if(dt !=  m_phvec.at(idx2).getDeltaTime())
    NCRYSTAL_THROW(CalcError,"NCrystal::PhononCalculator::conv time step of the spectra are different.");

  std::vector<double> phonon_spe;
  double start_time = m_phvec.at(idx1).getBeginTime() + m_phvec.at(idx2).getBeginTime();
  fc.fftconv(m_phvec.at(idx1).getSpectrum(),
             m_phvec.at(idx2).getSpectrum(),
             phonon_spe,
             m_phvec.at(idx2).getDeltaTime());


  m_phvec.push_back(PhononNumDyns(phonon_spe, start_time, dt, m_kt));
}

void NCrystal::PhononCalculator::calcfftconv(unsigned phnum)
{
  if(m_phvec.size()!=1)
    NCRYSTAL_THROW(CalcError,"NCrystal::PhononCalculator::calcfftconv single phonon data should be present in the m_phvec to run calcfftconv. ");

  FastConvolve f (m_phvec.at(0).getSpectrum().size()*(phnum+1));
  //zero phonon spectrum is empty
  for(unsigned i=1;i<phnum;i++)  {
    worker_conv(f, i);
  }
}


#ifndef NCrystal_PhononData_hh
#define NCrystal_PhononData_hh

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

#include "NCrystal/NCDefs.hh"
#include <vector>

namespace NCrystal {

  class PhononNumDyns {
  public:
    PhononNumDyns( const std::vector<double> &spec, double btime, double dt,
                   double kt );

    virtual ~PhononNumDyns();
    double getPhononDensity(double time) const;
    const std::vector<double> &getSpectrum() const;

    double getBeginTime() const {return m_btime;}
    double getEndTime() const {return m_etime;}
    double getDeltaTime() const {return m_dt;}
    double getSigma() const {return m_asym_sigma;}
    double getAsymCentre() const {return m_asym_centre;}

  private:
    std::vector<double> m_spec;
    size_t m_spec_size_minus_1;
    double m_btime;
    double m_etime;
    double m_dt;
    double m_inv_dt;
    double m_asym_centre;
    double m_asym_sigma;
    double m_kt;
  };

  typedef std::vector<PhononNumDyns > PhononVectorLow;

  class FastConvolve;

  class PhononCalculator {
  public:
    PhononCalculator(const std::vector<double> &spec, double btime, double dt, unsigned phnum, double kt);
    ~PhononCalculator();
    void calcfftconv(unsigned phnum);
    double interpolate(unsigned phonnum, double energytime) const
    {
      //inlined for efficiency
      nc_assert(phonnum<m_phvec.size());
      return m_phvec[phonnum].getPhononDensity(energytime);
    }
  private:
    double m_kt;
    PhononVectorLow m_phvec;
    unsigned m_fftnum;
    void worker_conv(const FastConvolve& fc, unsigned order);
  };

}

#endif

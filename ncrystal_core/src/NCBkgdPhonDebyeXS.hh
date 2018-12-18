#ifndef NCrystal_BkgdPhonDebyeXS_hh
#define NCrystal_BkgdPhonDebyeXS_hh

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

#include "NCrystal/NCDefs.hh"
#include <vector>
#include <utility>

namespace NCrystal {

  class Info;
  class BkgdPhonDebyeXS : public RCBase {
  public:

    //Use PhononDebye calculations to estimate background (non-Bragg)
    //contributions to the scattering cross-section. The initialisation (only to
    //be performed via createBkgdPhonDebyeXS) incurs a significant overhead, so
    //caching strategies might be relevant.

    //Energy range from 1e-5eV to ncmin(1,log(DBL_MAX)*kt*2)).  Outside this
    //interval, appropriate extrapolations (at high E, only if
    //extrapolate_from_peak is set).

    double getXS(const double& ekin) const;
    virtual ~BkgdPhonDebyeXS(){}

    //data and methods used just by initialisation code:
    BkgdPhonDebyeXS(double kt, bool extrapolate_from_peak );
    void accumInelastic(const std::vector<double> & xs, double frac);
    void setSaturatedXSAndInit(double xs);
    std::vector<double>& getEnergyVector() { return m_en_inel; }
  private:
    std::vector<double> m_en_inel, m_xs_inel;
    std::vector<std::pair<double,double> > m_k3;
    double m_saturated_xs;
    double m_k;
    double m_k2;
    bool m_extrapolate_from_peak;
    void cropIncompleteHighEBins();
  };

  RCHolder<const BkgdPhonDebyeXS> createBkgdPhonDebyeXS( const Info* ci,
                                                         unsigned nphonon = 0,
                                                         bool include_phonzeroinco = true,
                                                         bool only_phonzeroinco = false,
                                                         bool extrapolate_from_peak = true );
}

#endif

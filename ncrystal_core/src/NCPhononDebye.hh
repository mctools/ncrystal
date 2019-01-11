#ifndef NCrystal_PhononDebye_hh
#define NCrystal_PhononDebye_hh

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

#include "NCrystal/NCDefs.hh"
#include <vector>
#include <string>

namespace NCrystal {

  class PhononCalculator;

  //This class calculates the total inelastic scattering cross section in the
  //Debye approximation, harmonic incoherent approximation and phonon expansion approximation

  //Implementation of a highly efficient but less accurate Gaussian approximation [1]
  //to estimate the saturated expansion order in the class constructor
  //Sjolander, Arkiv for Fysik., Bd 14, nr 21, 1958
  //automatic order determination is enabled by setting max_phonon_order in the constructor to zero
  //
  //If phonzeroinco is 1, incoherent contributions from the zero'th order phonon
  //will be included (the coherent part is Bragg diffraction, which is handled
  //by PCBragg, SCBragg, or LCBragg). If it is 0, they will be excluded and if
  //it is 2, *only* those contributions will be included (must be used with
  //max_phonon_order=1).

  class NCRYSTAL_API PhononDebye  {
  public:
    PhononDebye(double debye_energy, double kt, const std::string & ele_name,
                unsigned max_phonon_order=0, int phonzeroinco = 1 );
    virtual ~PhononDebye();
    void doit(const std::vector<double> &ekin_vec, std::vector<double> &xs_vec,
              unsigned alpha_grid_size=50, unsigned beta_sym_grid_size=100);
    double getMSD() const;

    //Get max_phonon_order from constructor, unless it was 0, in which case the
    //estimated order can be retrieved here:
    unsigned getMaxPhononNum() const { return m_max_phononnum; }

  private:
    double getS(unsigned beta_index, double alpha) const;
    double interpolate(double a, double fa, double b, double fb, double x) const;
    double integrateAlphaInterval(double a1, double s1, double a2, double s2 ) const;
    void getSecondarySpectrum(double kin, double* spec) const;
    void getAlphaLimits(double kin, double beta, double &lower, double& upper) const;
    double g0knl(double eps) const;
    double g0barknl(double eps) const;
    double d0bar_power_knl(double eps) const;
    void calcSymSab(const PhononCalculator&cal,const std::vector<double> &alpha,
                    const std::vector<double> &beta, std::vector<double> &result) const;
    struct para{
      double x, z,factor, delta0_bar, kt;
    };

    double sigma_1 (const para& p, double n) const;
    double sigma_2 (const para& p, double n) const;


    //////////
    double m_debye, m_kt, m_dt;
    double m_msd, m_gamma0 , m_gamma0_bar, m_delta0_bar;
    std::string m_ele_name;
    std::vector<double> m_alpha;
    std::vector<double> m_beta;
    std::vector<double> m_sab;
    unsigned m_max_phononnum;
    double m_max_wl2ekin;
    int m_phonzeroinco;
  };
}

#endif

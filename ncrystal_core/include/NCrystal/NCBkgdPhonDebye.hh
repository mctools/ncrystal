#ifndef NCrystal_BkgdPhonDebye_hh
#define NCrystal_BkgdPhonDebye_hh

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

#include "NCrystal/NCScatterXSCurve.hh"

namespace NCrystal {

  class Info;
  class BkgdPhonDebyeXS;

  class NCRYSTAL_API BkgdPhonDebye : public ScatterXSCurve {
  public:

    //Provides the total inelastic scattering cross section in the Debye
    //approximation, harmonic incoherent approximation and phonon expansion
    //approximation.
    //
    //The parameter thermalise is passed on to the ScatterXSCurve base class and
    //determines how energy transfers are modelled when generating scatterings
    //and the nphonon parameter control the number of terms used in the phonon
    //expansion. Higher numbers potentially provide increased precision at the
    //cost of added initialisation time, but only to a certain point due to
    //numerical uncertainties. If the value is 0, the code will attempt to
    //automatically select a reasonable number of terms for the given crystal.
    //Finally, if include_phonzeroinco is true, incoherent contributions from
    //the zero'th order phonon will be included (the coherent part being Bragg
    //diffraction, which is handled by other classes). The
    //no_highe_extrapolation_from_peak flag is used to disable the high-energy
    //cropping and extrapolation to saturated XS.

    enum GenScatterMode { elastic, thermalise, modeldeltae };

    BkgdPhonDebye( const Info*,
                   GenScatterMode genscattermode = thermalise,
                   unsigned nphonon = 0,
                   bool include_phonzeroinco = true,
                   bool only_phonzeroinco = false,
                   bool no_highe_extrapolation_from_peak = false );

    virtual double crossSectionNonOriented(double ekin) const;
    virtual void generateScatteringNonOriented( double ekin, double& angle, double& de ) const;
    virtual void generateScattering( double ekin, const double (&neutron_direction)[3],
                                     double (&resulting_neutron_direction)[3], double& delta_ekin ) const;

    //Check if BkgdPhonDebye can be created from the info:
    static bool hasSufficientInfo(const Info*);

  protected:
    virtual ~BkgdPhonDebye();
    const BkgdPhonDebyeXS* m_debphonmodel;
    GenScatterMode m_genscatmode;
  };
}

#endif

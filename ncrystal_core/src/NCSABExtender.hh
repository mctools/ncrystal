#ifndef NCrystal_SABExtender_hh
#define NCrystal_SABExtender_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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
#include "NCFreeGasUtils.hh"

namespace NCrystal {

  namespace SAB {

    //Class providing models used to extend cross-sections and sampling
    //capabilities of an S(alpha,beta) table to large energies. In the simplest
    //scenario it would wrap a single-target free-gas model, but other models
    //might be more suitable depending on the kernel in question. It could be
    //alternative models (e.g. short collision time approximation), or
    //multi-target models for when the kernel is provided for an entire molecule
    //or unit cell, describing more than a single element, or the isotopic
    //weight differences in the target elements are non-negligible.

    class SABExtender : private MoveOnly {
    public:
      SABExtender() = default;
      virtual ~SABExtender();

      //Get cross-section or sample S(alpha,beta) over allowed kinematic region
      //at given incident neutron energy:
      virtual double crossSection(double ekin) const = 0;
      virtual PairDD sampleAlphaBeta(RandomBase&, double ekin) const = 0;
    };

    class SABFGExtender : public SABExtender {
    public:
      //Extend with single target type FreeGas model.
      SABFGExtender( double temp_kelvin, double target_mass_amu, SigmaBound );
      SABFGExtender( double temp_kelvin, double target_mass_amu, SigmaFree );
      virtual ~SABFGExtender();
      double crossSection(double ekin) const override;
      PairDD sampleAlphaBeta(RandomBase&, double ekin) const override;
    private:
      FreeGasXSProvider m_xsprovider;
      double m_t, m_m;
    };

  }

}

#endif

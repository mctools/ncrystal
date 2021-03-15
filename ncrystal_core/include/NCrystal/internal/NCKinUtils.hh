#ifndef NCrystal_KinUtils_hh
#define NCrystal_KinUtils_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2021 NCrystal developers                                   //
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

#include "NCrystal/NCTypes.hh"

//Various utilities related to kinematics, including those related to the
//S(alpha,beta) formalism.

namespace NCrystal {

  //Get deltaE and mu=cos(scattering_angle) from (alpha,beta). It is recommended
  //to check muIsotropicAtBeta before calling this (for numerical safety, but it
  //might also save the effort of sampling alpha).
  PairDD convertAlphaBetaToDeltaEMu(double alpha, double beta, NeutronEnergy ekin, double kT );
  inline PairDD convertAlphaBetaToDeltaEMu(PairDD alphabeta, NeutronEnergy ekin, double kT )
  {
    return convertAlphaBetaToDeltaEMu(alphabeta.first,alphabeta.second,ekin,kT);
  }

  //Get kinematically accessible alpha region for given ekin/kT and beta:
  PairDD getAlphaLimits( double ekin_div_kT, double beta );

  //Near kinematic endpoint. alpha-, alpha, and alpha+ might become numerically
  //indistinguishable at floating point precision, but mu=cos(scattering_angle)
  //should be sampled isotropically in this limit (since
  //S(alpha-,beta)~=S(alpha+,beta)). For consistency, always use this function
  //to detect this degenerate case:
  inline bool muIsotropicAtBeta( double beta, double ekin_div_kT )
  {
    nc_assert(beta>=-ekin_div_kT);
    nc_assert(ekin_div_kT>=0.0);
    constexpr double lim(-1.0+1e-14);
    return beta <= ekin_div_kT*lim;
  }


}

#endif

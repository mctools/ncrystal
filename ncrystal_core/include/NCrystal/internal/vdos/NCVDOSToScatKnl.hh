#ifndef NCrystal_VDOSToScatKnl_hh
#define NCrystal_VDOSToScatKnl_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2025 NCrystal developers                                   //
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

#include "NCrystal/internal/sab/NCScatKnlData.hh"
#include "NCrystal/internal/vdos/NCVDOSGn.hh"

namespace NCRYSTAL_NAMESPACE {

  /////////////////////////////////////////////////////////////////////////////
  // Expand a vibrational density of state (VDOS) spectrum into a full-blown //
  // scattering kernel.                                                      //
  /////////////////////////////////////////////////////////////////////////////

  using ScaleGnContributionFct = std::function<double(unsigned)>;
  ScatKnlData createScatteringKernel( const VDOSData&,
                                      unsigned vdosluxlvl = 3,//0 to 5, affects binning, Emax, etc.
                                      double targetEmax = 0.0,//if 0, will depend on luxlvl. Error if set to unachievable value.
                                      const VDOSGn::TruncAndThinningParams ttpars = VDOSGn::TruncAndThinningChoices::Default,
                                      ScaleGnContributionFct = nullptr,
                                      Optional<unsigned> override_max_order = NullOpt );

  //The ScaleGnContributionFct argument can be used to apply a scale factor to
  //the Gn contribution, at the point where it is used to add a contribution
  //into S(alpha,beta). This can for instance be used in the scenario where the
  //coherent single-phonon contribution is modelled elsewhere, and therefore the
  //G1 contribution should be reduced to take out sigma_coherent. In this case,
  //the scaling function should return sigma_incoh/(sigma_coh+sigma_incoh) for
  //the argument n==1, and 1.0 for n>1.

  //Internal functions, exposed here for testing:
  VectD setupAlphaGrid( double kT, double msd, double alphaMax, unsigned npts );
  VectD setupBetaGrid( const VDOSGn& Gn, double betaMax, unsigned luxlvl, unsigned override_nbins );
  PairDD rangeXNexpMX(unsigned n, double eps, double accuracy = 1e-13 );
  PairDD findExtremeSABPointWithinAlphaPlusCurve(double E_div_kT, PairDD alphaRange, PairDD betaRange);
  bool sabPointWithinAlphaPlusCurve(double E_div_kT, double alpha, double beta );

}

#endif

#ifndef NCrystal_VDOSToScatKnl_hh
#define NCrystal_VDOSToScatKnl_hh

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

#include "NCrystal/internal/NCScatKnlData.hh"
#include "NCrystal/internal/NCVDOSGn.hh"

namespace NCrystal {

  /////////////////////////////////////////////////////////////////////////////
  // Expand a vibrational density of state (VDOS) spectrum into a full-blown //
  // scattering kernel.                                                      //
  /////////////////////////////////////////////////////////////////////////////

  ScatKnlData createScatteringKernel( const VDOSData&,
                                      unsigned vdosluxlvl = 3,//0 to 5, affects binning, Emax, etc.
                                      double targetEmax = 0.0,//if 0, will depend on luxlvl. Error if set to unachievable value.
                                      const VDOSGn::TruncAndThinningParams ttpars = VDOSGn::TruncAndThinningChoices::Default );

  //Internal functions, exposed here for testing:
  VectD setupAlphaGrid( double kT, double msd, double alphaMax, unsigned npts );
  VectD setupBetaGrid( const VDOSGn& Gn, double betaMax, unsigned luxlvl, unsigned override_nbins );
  PairDD rangeXNexpMX(unsigned n, double eps, double accuracy = 1e-13 );
  PairDD findExtremeSABPointWithinAlphaPlusCurve(double E_div_kT, PairDD alphaRange, PairDD betaRange);
  bool sabPointWithinAlphaPlusCurve(double E_div_kT, double alpha, double beta );

}

#endif

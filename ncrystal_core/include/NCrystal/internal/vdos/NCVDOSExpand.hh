#ifndef NCrystal_VDOSExpand_hh
#define NCrystal_VDOSExpand_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

#include "NCrystal/internal/vdos/NCVDOSGn.hh"
#include "NCrystal/internal/vdos/NCVDOSLux.hh"
#include "NCrystal/internal/utils/NCRect.hh"

namespace NCRYSTAL_NAMESPACE {

  class VDOSData;

  namespace VDOS {

    struct GnExpansion final {
      //Class representing the result of a multiphonon expansion into Gn
      //functions. Not just the Gn functions themselves, but also other
      //information which is needed in order to encode them in actual
      //S(alpha,beta) kernels.

      // The Gn functions of the expanded phonon orders:
      VDOSGn Gn;

      // Sjolander's 2W factor divided by alpha. This determines the scale of
      // phonon weights along alpha, and is proportional to both kT and
      // mean-squared-displacements:
      double alpha2x;

      //The Emax value is what guided how many phonon orders were needed in the
      //expansion. In some cases it might actually be lower than the value aimed
      //for. The value here is what is actually covered by the Gn functions in
      //this expansion:
      NeutronEnergy suggestedEmax;

      //If encoding the phonon expansions into an S(alpha,beta) table, this is
      //the minimum rectangular range of (alpha,beta) values that table must
      //cover (Rectangle x is alpha, Rectangle y is beta):
      Rectangle sabRange;

      //Diagnostic only (see doc/claude_session_vdos_fma_reprod.md): the
      //full (alpha,beta) range of each individual phonon order n=1..maxOrder
      //(i.e. before intersecting with the kinematically accessible region),
      //indexed abRanges.at(n-1). Used to bisect cross-platform/compiler
      //reproducibility divergences in the order-growth stopping criterion
      //in expandVDOSToGnFcts, by comparing per-order ranges instead of only
      //the final union in sabRange:
      std::vector<Rectangle> abRanges;
    };

    //Perform Sjolander expansion. The targetEmax default value is given by
    //vdoslux, but can be overridden:
    GnExpansion
    expandVDOSToGnFcts( const VDOSData&,
                        VDOSLux vdoslux = VDOSLux(),
                        Optional<NeutronEnergy> targetEmax = NullOpt);

  }

}

#endif

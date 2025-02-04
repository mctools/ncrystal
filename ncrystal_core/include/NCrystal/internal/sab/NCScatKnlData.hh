#ifndef NCrystal_ScatKnlData_hh
#define NCrystal_ScatKnlData_hh

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

#include "NCrystal/interfaces/NCSABData.hh"
#include "NCrystal/internal/utils/NCSpan.hh"

///////////////////////////////////////////////////////////////////////////////////
//                                                                               //
// Infrastructure for specifying, and performing basic validation on, scattering //
// kernels in a variety of formats: S(alpha,beta), S(q,w), symmetric/asymmetric, //
// with or without ("scaled") detailed balance factors.                          //
//                                                                               //
// The ScatKnlData struct is similar to SABData (cf. NCSABData.hh), but with     //
// more options for how the scattering kernel can be defined. For simplicity,    //
// algorithms working on scattering kernels are in general only required to work //
// with SABData objects.                                                         //
//                                                                               //
///////////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  struct ScatKnlData : private MoveOnly {
    VectD alphaGrid, betaGrid, sab;
    Temperature temperature = Temperature{-1.0};
    SigmaBound boundXS = SigmaBound{-1.0};
    AtomMass elementMassAMU = AtomMass{-1.0};
    enum class KnlType { SAB,        //Standard S(alpha,beta), fields have same meaning as on SABData.
                         SCALED_SAB, //Values in sab table are actually S'(alpha,beta)=S(alpha,beta)*exp(beta/2)
                         SCALED_SYM_SAB,//Same + S'(alpha,beta) is an even function in beta, and the kernel is
                         //only specified directly for non-negative beta values.
                         SQW,//alpha/beta/sab values are actually Q/omega/S(Q,omega) values.
                         Unspecified };
    KnlType knltype = KnlType::Unspecified;
    double suggestedEmax = 0;//optional, 0 means not available.
    bool betaGridOptimised = false;//set to true if beta grid is known to not require thickening before sampling.
  };

  struct ScatKnlDataView {
    //Similar to ScatKnlData, but a read-only non-owning view which is able to
    //wrap both ScatKnlData and SABData objects to look like ScatKnlData
    //objects. That way, both ScatKnlData and SABData objects can be analysed
    //using the same functions.
    ScatKnlDataView(const ScatKnlData&);
    ScatKnlDataView(const SABData&);
    const Span<const double> alphaGrid, betaGrid, sab;
    const Temperature temperature;
    const SigmaBound boundXS;
    const AtomMass elementMassAMU;
    const ScatKnlData::KnlType knltype;
    double suggestedEmax;
  };

  //Validate consistency of kernel's data fields. Throws BadInput exceptions
  //if not ok. Skips most expensive checks if cheapOnly.
  void validateScatKnlData( const ScatKnlDataView& );

}


////////////////////////////
// Inline implementations //
////////////////////////////

inline NCrystal::ScatKnlDataView::ScatKnlDataView(const NCrystal::ScatKnlData& o)
  : alphaGrid(o.alphaGrid),
    betaGrid(o.betaGrid),
    sab(o.sab),
    temperature(o.temperature),
    boundXS(o.boundXS),
    elementMassAMU(o.elementMassAMU),
    knltype(o.knltype),
    suggestedEmax(o.suggestedEmax)
{
}

inline NCrystal::ScatKnlDataView::ScatKnlDataView(const NCrystal::SABData& o)
  : alphaGrid(o.alphaGrid()),
    betaGrid(o.betaGrid()),
    sab(o.sab()),
    temperature(o.temperature()),
    boundXS(o.boundXS()),
    elementMassAMU(o.elementMassAMU()),
    knltype(ScatKnlData::KnlType::SAB),
    suggestedEmax(o.suggestedEmax())
{
}

#endif

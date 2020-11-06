#ifndef NCrystal_LCRefModels_hh
#define NCrystal_LCRefModels_hh

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

#include "NCrystal/NCScatter.hh"
#include "NCrystal/internal/NCVector.hh"
#include "NCrystal/internal/NCRotMatrix.hh"

namespace NCrystal {

  class LCBraggRef : public Scatter {
  public:
    //Simple but very slow implementation of layered crystals. Mainly provided
    //as a reference (should give increasingly better result with higher
    //nsample).
    LCBraggRef(Scatter* scbragg, Vector lcaxis_lab, unsigned nsample = 1000);
    virtual ~LCBraggRef();
    virtual void domain(double& ekin_low, double& ekin_high) const;
    virtual double crossSection( double ekin, const double (&indirraw)[3] ) const;
    virtual void generateScattering( double ekin,
                                     const double (&indirraw)[3],
                                     double (&outdir)[3],
                                     double& delta_ekin ) const;
  private:
    RCHolder<const Scatter> m_sc;
    Vector m_lcaxislab;
    unsigned m_nsample;
    unsigned m_nsampleprime;
  };

  class LCBraggRndmRot : public Scatter {
  public:
    //Like LCBraggRef, but using random crystallite rotations even for
    //crossSection calls - and reusing the same orientations in a subsequent
    //call to generateScattering. Again, this is mainly provided as a reference
    //and not really recommended for general usage.
    LCBraggRndmRot(Scatter* scbragg, Vector lcaxis_lab, unsigned nsample = 1);
    virtual ~LCBraggRndmRot();
    virtual void domain(double& ekin_low, double& ekin_high) const;
    virtual double crossSection( double ekin, const double (&indirraw)[3] ) const;
    virtual void generateScattering( double ekin,
                                     const double (&indirraw)[3],
                                     double (&outdir)[3],
                                     double& delta_ekin ) const;
  private:
    RCHolder<const Scatter> m_sc;
    Vector m_lcaxislab;
    unsigned m_nsample;
    mutable struct Cache {
      std::vector<PhiRot> rotations;//rotations sampled
      VectD xscommul;//cross-sections at the sampled rotations.
    } cache;
  };

}

#endif


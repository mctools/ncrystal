#ifndef NCrystal_PCBragg_hh
#define NCrystal_PCBragg_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2023 NCrystal developers                                   //
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

#include "NCrystal/NCInfo.hh"
#include "NCrystal/NCProcImpl.hh"

namespace NCrystal {

  class PCBragg final : public ProcImpl::ScatterIsotropicMat {
  public:

    //Calculates Bragg diffraction in a powdered (or non-textured
    //polycrystalline) material. Does not account for texture, grain size,
    //dspacing-deviations and other similar effects.

    const char * name() const noexcept final { return "PCBragg"; }

    //Constructor:
    PCBragg( const Info& );

    //Specialised constructors taking (dspacing,fsquared*multiplicity) pairs.
    //Either needs structure info, or just v0*n_atoms, unit cell volume in Aa^3
    //and number atoms per unit cell:

    using VectDFM = std::vector<PairDD>;
    PCBragg( const StructureInfo&, VectDFM&& );
    PCBragg( double v0_times_natoms, VectDFM&& );

    //There is a maximum wavelength at which Bragg diffraction is possible, so
    //lower energy bound will reflect this (upper bound is infinity):
    EnergyDomain domain() const noexcept final;

    CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ) const final;
    ScatterOutcomeIsotropic sampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const final;

    //Two PCBragg instances can be merged by merging the plane lists:
    std::shared_ptr<Process> createMerged( const Process& other,
                                           double scale_self,
                                           double scale_other ) const override;


    //Empty, no planes:
    PCBragg( no_init_t ) {}

  protected:
    Optional<std::string> specificJSONDescription() const override;
  private:
    CosineScatAngle genScatterMu(RNG&, NeutronEnergy ekin) const;
    std::size_t findLastValidPlaneIdx( NeutronEnergy ekin) const;
    NeutronEnergy m_threshold = NeutronEnergy{kInfinity};
    VectD m_2dE;
    VectD m_fdm_commul;
    void init( const StructureInfo&, VectDFM&& );
    void init( double v0_times_natoms, VectDFM&& );
  };

}

#endif

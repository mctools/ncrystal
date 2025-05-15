#ifndef NCrystal_PowderBragg_hh
#define NCrystal_PowderBragg_hh

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

#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/interfaces/NCProcImpl.hh"
#include "NCrystal/internal/utils/NCExtraTypes.hh"

namespace NCRYSTAL_NAMESPACE {

  class PowderBragg final : public ProcImpl::ScatterIsotropicMat {
  public:

    //Calculates Bragg diffraction in a powdered (or non-textured
    //polycrystalline) material. Does not account for texture, grain size,
    //dspacing-deviations and other similar effects.

    const char * name() const noexcept override { return "PowderBragg"; }

    //Constructor:
    PowderBragg( PreparedPowderInputData&& );

    //Other constructors are in principle already covered by
    //PreparedPowderInputData, but keeping them for now for backwards
    //compatiblity:
    using VectDFM = std::vector<PairDD>;
    PowderBragg( double v0_times_natoms, VectDFM&& );
    PowderBragg( const StructureInfo&, VectDFM&& );
    PowderBragg( const Info& );
    PowderBragg( no_init_t );//Empty, no planes

    //There is a maximum wavelength at which Bragg diffraction is possible, so
    //lower energy bound will reflect this (upper bound is infinity):
    EnergyDomain domain() const noexcept override;

    CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ) const override;
    ScatterOutcomeIsotropic sampleScatterIsotropic(CachePtr&,
                                                   RNG&,
                                                   NeutronEnergy ) const override;

    //Two PowderBragg instances can be merged by merging the plane lists:
    std::shared_ptr<Process> createMerged( const Process& other,
                                           double scale_self,
                                           double scale_other ) const override;



#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
    bool isPureElasticScatter() const override { return true; }
    std::pair<CrossSect,ScatterOutcome>
    evalXSAndSampleScatter( CachePtr&,
                            RNG&,
                            NeutronEnergy,
                            const NeutronDirection& ) const override;
    std::pair<CrossSect,ScatterOutcomeIsotropic>
    evalXSAndSampleScatterIsotropic(CachePtr&,
                                    RNG&,
                                    NeutronEnergy ) const override;
    void evalManyXSIsotropic( CachePtr&,
                              const double* ekin,
                              std::size_t N,
                              double* out_xs ) const override;
#endif

  protected:
    Optional<std::string> specificJSONDescription() const override;
  private:
    CosineScatAngle genScatterMu(RNG&,
                                 NeutronEnergy ekin,
                                 std::size_t idx) const;
    std::size_t findLastValidPlaneIdx( NeutronEnergy ekin) const;
    NeutronEnergy m_threshold = NeutronEnergy{kInfinity};
    VectD m_2dE;
    VectD m_fdm_commul;
    void init( PreparedPowderInputData&& );
  };

}

#endif

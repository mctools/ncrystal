#ifndef NCrystal_ElIncScatter_hh
#define NCrystal_ElIncScatter_hh

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

#include "NCrystal/interfaces/NCProcImpl.hh"

namespace NCRYSTAL_NAMESPACE {

  class ElIncXS;
  class Info;

  struct ElIncScatterCfg {
    //Outside due to https://stackoverflow.com/questions/17430377/
    bool use_sigma_incoherent = true;
    bool use_sigma_coherent = false;
    double scale_factor = 1.0;
  };

  class ElIncScatter final : public ProcImpl::ScatterIsotropicMat {
  public:

    //Model elastic-incoherent scatterings based on the Debye model with
    //isotropic atomic displacements. For each element in the crystalline
    //material, the mean-squared-displacement must be provided, along with the
    //corresponding bound incoherent scattering cross-section and the scale. The
    //scale will often be the fraction of the element (by count). For more
    //details, see section 2.3 of https://doi.org/10.1016/j.cpc.2019.07.015 and
    //comments in the NCElIncXS.hh header file.

    const char * name() const noexcept override { return "ElIncScatter"; }
    virtual ~ElIncScatter();

    //Constructor similar to the ElIncXS constructor:
    ElIncScatter( const VectD& elements_meanSqDisp,
                  const VectD& elements_boundincohxs,
                  const VectD& elements_scale );

    //Construct from Info. Using the cfg parameters, it is possible to apply an
    //overall cross section scale-factor (>0.0), and to add either of
    //sigma_incoherent and/or sigma_coherent. The latter is only appropriate
    //under the incoherent approximation, and by default only the former is
    //used. The Debye Waller factors will be determined via AtomInfo objects if
    //available, otherwise from DynInfo objects. In that case, elements whose
    //DynInfo is unable to provide Debye Waller factors, e.g. sterile/freegas
    //and for now scatknl, will be ignored (if all elements are ignored, an
    //error is raised).
    ElIncScatter( const Info&, const ElIncScatterCfg& cfg = ElIncScatterCfg() );

    CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ) const override;
    ScatterOutcomeIsotropic sampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const override;

    //Simple additive merge:
    std::shared_ptr<Process> createMerged( const Process& other,
                                           double scale_self,
                                           double scale_other ) const override;


    //Specialised constructor providing internal state directly:
    ElIncScatter( std::unique_ptr<ElIncXS> );

    //Helper function for determining if ElincScatter(info,...) can proceed
    //(i.e. if there is enough info to estimate Debye-Waller factors):
    static bool hasSufficientInfo( const Info&, const ElIncScatterCfg& cfg = ElIncScatterCfg() );

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
    bool isPureElasticScatter() const override { return true; }
    std::pair<CrossSect,ScatterOutcomeIsotropic>
    evalXSAndSampleScatterIsotropic(CachePtr&, RNG&, NeutronEnergy ) const override;
    void evalManyXSIsotropic( CachePtr&, const double* ekin, std::size_t N, double* out_xs ) const override;
#endif

  protected:
    Optional<std::string> specificJSONDescription() const override;
  private:
    std::unique_ptr<ElIncXS> m_elincxs;
  };
}

#endif

#ifndef NCrystal_SCBragg_hh
#define NCrystal_SCBragg_hh

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
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/interfaces/NCSCOrientation.hh"

namespace NCRYSTAL_NAMESPACE {

  class PlaneProvider;

  class SCBragg final : public ProcImpl::ScatterAnisotropicMat {
  public:

    //Calculates Bragg diffraction in a single crystals with given orientation
    //and mosaicity, assuming a Gaussian distribution of crystallite
    //orientations.
    //
    //Parameter delta_d is the deviation of d-spacing in a non-ideal crystal.
    //
    //plane_provider optional object which allows one to override the planes
    //otherwise present in the Info object. It will only be used during
    //initialisation in the constructor, and the SCBragg instance will *not*
    //assume ownership of it.
    //
    //For a description of the prec and ntrunc parameters, see NCGaussMos.hh.
    SCBragg( const Info&,
             const SCOrientation&,
             MosaicityFWHM,
             double delta_d = 0,
             PlaneProvider * plane_provider = nullptr,
             double prec = 1e-3, double ntrunc = 0.0 );

    const char * name() const noexcept override { return "SCBragg"; }

    virtual ~SCBragg();

    //There is a maximum wavelength at which Bragg diffraction is possible, so
    //lower bound will reflect this (upper bound is infinity):
    EnergyDomain domain() const noexcept override;

    CrossSect crossSection(CachePtr&, NeutronEnergy, const NeutronDirection& ) const override;
    ScatterOutcome sampleScatter(CachePtr&, RNG&, NeutronEnergy, const NeutronDirection& ) const override;

#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
    bool isPureElasticScatter() const override { return true; }
#endif

  protected:
    Optional<std::string> specificJSONDescription() const override;
  private:
    struct pimpl;
    std::unique_ptr<pimpl> m_pimpl;
  };

}

#endif

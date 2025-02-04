#ifndef NCrystal_SABSampler_hh
#define NCrystal_SABSampler_hh

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
#include "NCrystal/internal/sab/NCSABExtender.hh"

namespace NCRYSTAL_NAMESPACE {

  class SABSamplerAtE : private NoCopyMove {
    //For sampling (alpha,beta) values at a given energy, Ei, or lower. This
    //is intended to be implemented with the rejection method, taking
    //advantage of the fact that the kinematically accessible region of the
    //(alpha,beta)-plane, is a strictly expanding region with neutron
    //energy. Due to acceptance-rate consideration of the MC rejection method,
    //implementations are expected to be better performing when ekin is not
    //too much lower than Ei (which is why we need a whole grid of Ei values,
    //rather than just one with Ei=Emax).
  public:
    virtual PairDD sampleAlphaBeta(double ekin_div_kT, RNG&) const = 0;
    virtual ~SABSamplerAtE() = default;
  };

  class SABSampler final : private MoveOnly {
  public:

    struct EGridMargin {
      //The value (typically 1 or a value a bit larger than 1, e.g. 1.05),
      //specifies a safety margin for what egrid point to actually sample the
      //overlay distribution. A value of 1.05 would indicate that Egrid must be
      //at least 1.05*ekin. A value of 1.0 indicates no additional margin is
      //required.
      double value = 1.0;
      constexpr EGridMargin(double vv = 1.0) noexcept : value(vv) {}//for C++11
    };


    //Wanted to use std::vector, but because of (older?) VSCode not supporting
    //move-only types in std::vector, we use SmallVector for this:
    using SABSamplerAtEList = SmallVector<std::unique_ptr<SABSamplerAtE>,1>;
    //using SABSamplerAtEList = std::vector<std::unique_ptr<SABSamplerAtE>>;

    void setData( Temperature temperature,
                  VectD&& egrid,
                  SABSamplerAtEList&&,
                  std::shared_ptr<const SAB::SABExtender>,
                  double xsAtEmax,
                  EGridMargin );

    SABSampler( Temperature temperature,
                VectD&& egrid,
                SABSamplerAtEList&&,
                std::shared_ptr<const SAB::SABExtender>,
                double xsAtEmax,
                EGridMargin );

    SABSampler() = default;//invalid instance.
    ~SABSampler();

    //Sample (alpha,beta) values directly:
    PairDD sampleAlphaBeta( NeutronEnergy, RNG&) const;

    //Convenience (calls sampleAlphaBeta, then converts):
    PairDD sampleDeltaEMu( NeutronEnergy, RNG& rng) const;

    //Move ok:
    SABSampler( SABSampler&& ) = default;
    SABSampler& operator=( SABSampler&& ) = default;

  private:
    VectD m_egrid;
    SABSamplerAtEList m_samplers;
    double m_kT = 0.0;
    std::shared_ptr<const SAB::SABExtender> m_extender;
    double m_xsAtEmax = 0.0, m_k1 = 0.0, m_k2 = 0.0;
    PairDD sampleHighE(NeutronEnergy, RNG&) const;
    EGridMargin m_egridMargin;
  };
}

#endif

#ifndef NCrystal_MMC_Source_hh
#define NCrystal_MMC_Source_hh

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

#include "NCrystal/internal/minimc/NCMMC_Defs.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    class Geometry;

    struct SourceMetaData {

      //Some sources are safe to use concurrently:
      bool concurrent = false;

      //Set if the source provides a well defined "nominal" (mean) energy or
      //direction of emitted neutrons. Note that if user specified a wavelength,
      //rather than energy, it will be the mean wavelength rather than mean
      //energy that is represented in the meanEnergy variable:
      Optional<NeutronDirection> meanDirection;
      Optional<NeutronEnergy> meanEnergy;

      //Set if the source provides the same fixed energy or direction of all
      //emitted neutrons:
      Optional<NeutronDirection> fixedDirection;
      Optional<NeutronEnergy> fixedEnergy;

      //Neutron energies formatted in a string (e.g. "1Aa", "25meV", "2eV",
      //"(0.1-1)eV", etc.) suitable for e.g. plot labels):
      std::string energyDescription;

      //Approximate range of energies, needed for tally ranges. Some Sources
      //will have infinite range (e.g. thermal spectrum sources), but should
      //then try to capture "most" generated neutrons. Currently this is not
      //defined rigourously, but is roughly "3 sigma".
      using ERange = std::pair<NeutronEnergy,NeutronEnergy>;
      Optional<ERange> approxERange;
    };

    class Source : NoCopyMove {
    public:
      //Sources (or "particle guns") are responsible for filling up baskets of
      //pending neutrons. When needed, a source is given a basket, in which it
      //should fill more neutrons. Note that the basket passed to this method
      //might not be empty, and the existing entries should be left alone! If
      //the call does not completely fill the basket, it will be assumed that
      //the source has run out (which should never happen for Infinite sources).
      virtual void fillBasket( RNG&, NeutronBasket& ) = 0;

      //Provide source metaData:
      virtual const SourceMetaData& metaData() const = 0;

      //Query if generated particles might be outside the provided geometry (and
      //therefore not need an initial propagation step to the volume):
      virtual bool particlesMightBeOutside( const Geometry& ) const = 0;

      //Query total weight of all particles already provided by source:
      virtual ParticleCountSum particlesProvided() const = 0;

      //Serialisation of source configuration, as original source cfg-string, or
      //as a JSON object with more direct access to individual values. The JSON
      //object is a dictionary like (where one of the params should always be
      //the name of the type of source (i.e. "name":"constant").  Specifically
      //the layout is: { "cfgstr" : "...", "decoded" : { "k" : "v", ...} }

      virtual void toJSON(std::ostream&) const = 0;
      virtual void toString(std::ostream&) const = 0;
    };

    //NB: SourcePtr's are not const pointers:
    using SourcePtr = shared_obj<Source>;

    //For convenience + abstraction, the following function is used to transform
    //a source description string into an actual source object.
    SourcePtr createSource( const StrView& );

    //Documentation of all options as JSON dictionary:
    void sourceOptsDocsToJSON( std::ostream& );


  }
}

#endif

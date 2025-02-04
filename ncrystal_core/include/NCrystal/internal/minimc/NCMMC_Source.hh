#ifndef NCrystal_MMC_Source_hh
#define NCrystal_MMC_Source_hh

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

#include "NCrystal/internal/minimc/NCMMC_Basket.hh"
#include "NCrystal/core/NCTypes.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    class Geometry;

    struct SourceMetaData {
      //A source can be infinite or finite, and a finite source may or may not
      //know the total number of neutrons it can provide.
      std::string description;
      bool isInfinite = true;
      Optional<std::size_t> totalSize;//for finite sources of known size
      bool concurrent = false;//If it is safe to use in multithreaded context.
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

      //Provide source metaData. This function will only be called a few times,
      //and is not performance critical.
      virtual SourceMetaData metaData() const = 0;

      //Query if generated particles might be outside the provided geometry (and
      //therefore not need an initial propagation step to the volume):
      virtual bool particlesMightBeOutside( const Geometry& ) const = 0;

    };

    //NB: SourcePtr's are not const pointers:
    using SourcePtr = shared_obj<Source>;

    //For convenience + abstraction, the following function is used to transform
    //a source description string into an actual source object.
    SourcePtr createSource( const char * );

  }
}

#endif

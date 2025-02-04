#ifndef NCrystal_MMC_Geom_hh
#define NCrystal_MMC_Geom_hh

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
#include "NCrystal/internal/utils/NCSpan.hh"
#include "NCrystal/internal/utils/NCVector.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    class Geometry : NoCopyMove {
    public:
      //Geometry interface works on entire NeutronBasket's for efficiency,
      //writing results into provided data arrays.

      //Finds the distance a particle needs to propagate forward in space in
      //order to enter the volume. Returns 0.0 for neutrons already inside the
      //volume, and -1.0 for neutrons that are outside the volume and whose path
      //does not intersect the volume.
      virtual void distToVolumeEntry( const NeutronBasket&, Span<double> ) const = 0;

      //Finds the distance out of a volume. This is usually the most performance
      //critical geometry function. It is undefined behaviour to invoke this
      //function unless all neutrons in the basket are already inside the
      //volume.
      virtual void distToVolumeExit( const NeutronBasket&, Span<double> ) const = 0;

      //Check if a particular point is inside (scalar method meant to be used by
      //source-geometry checks, not neutron baskets):
      virtual bool pointIsInside( const Vector& ) const = 0;

    };

    using GeometryPtr = shared_obj<const Geometry>;

    //For convenience + abstraction, the following function is used to transform
    //a geometry description string into an actual geometry object.
    GeometryPtr createGeometry( const char * );

  }
}

#endif

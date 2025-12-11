#ifndef NCrystal_MMC_Geom_hh
#define NCrystal_MMC_Geom_hh

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
      //does not intersect the volume. As a special case, a neutron skirting
      //along the surface of the volume should return 0.0 if the surface locally
      //has a flat edge and the neutron direction is in the plane of that edge,
      //or -1.0 if the surface locally has a curved edge and the neutron is
      //about to leave the volume.
      //
      //If offset>0, the first offset entries in both the basket and the
      //destination buffer are ignored. Fixme: do we need to use span here?
      virtual void distToVolumeEntry( const NeutronBasket&,
                                      Span<double>,
                                      size_t offset ) const = 0;

      //Finds the distance out of a volume. This is usually the most performance
      //critical geometry function. It is undefined behaviour to invoke this
      //function unless all neutrons in the basket are already inside the
      //volume. If the neutron is propagated this distance forward along its
      //direction, it should afterwards be heading out of the shape (in other
      //words, a neutron sitting directly on the surface and a direction
      //pointing into the shape should not have a distToVolumeExit of 0).
      //FIXME: Double check the above is exactly valid also for a sphere, then
      //our scenario for pencil beam into a sphere does not have to move src_z
      //to (1-1e-14) times sphere_radius.
      virtual void distToVolumeExit( const NeutronBasket&, Span<double> ) const = 0;

      //Check if a particular point is inside (scalar method meant to be used by
      //source-geometry checks, not neutron baskets):
      virtual bool pointIsInside( const Vector& ) const = 0;

      //Serialisation, as original geometry cfg-string, or as a JSON object with
      //more direct access to individual values. The JSON object is a dictionary
      //like (where one of the params should always be the name of the type of
      //volume (i.e. "name":"sphere"):
      //{ "cfgstr" : "...", "decoded" : { "k" : "v", ...} }

      virtual void toJSON(std::ostream&) const = 0;
      virtual void toString(std::ostream&) const = 0;

      //fixme: Would it be useful to have a method returning a bounding sphere?
      //And perhaps also an inscribed sphere?
    };

    using GeometryPtr = shared_obj<const Geometry>;

    //For convenience + abstraction, the following function is used to transform
    //a geometry description string into an actual geometry object.
    GeometryPtr createGeometry( const char * );

  }
}

#endif

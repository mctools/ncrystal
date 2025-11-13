#ifndef NCrystal_ExtraTypes_hh
#define NCrystal_ExtraTypes_hh

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

#include "NCrystal/core/NCTypes.hh"

///////////////////////////////////////////////////
//                                               //
//  Extra types not intended for the public API  //
//                                               //
///////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace PowderBraggInput {

    //Common data structures which should contain everything (and not much more)
    //needed for Bragg diffraction powder processes. The MergedData struct
    //contains |F|^2 and plane multiplicity already combined by multiplication,
    //while the Data struct keeps them separate. A MergedData object has enough
    //information for an ideal Bragg powder model without extinction, while a
    //model including extinction will need a Data object instead.

    struct CellData {
      //For Powder Bragg processes, not much knowledge is needed about the cell,
      //since the structure factors already encapsulate most of the information
      //about the atoms:
      double volume = 0.0;//Aa^3
      unsigned n_atoms = 0.0;//Number of atoms per unit cell
    };

    template<class TPlane>
    struct DataImpl : private MoveOnly {
      using plane_t = TPlane;
      using PlaneList = std::vector<TPlane>;
      PlaneList planes;
      CellData cell;
    };

    struct Plane {
      double dsp ;//d-spacing [Angstrom]
      double fsq;//structure factor squared [barn]
      double mult;//multiplicity [integral, in double for efficiency]
    };
    using Data = DataImpl<Plane>;

    struct MergedPlane {
      double dsp;//d-spacing [Angstrom]
      double fsqmult;//structure factor squared times multiplicity [barn]
    };
    using MergedData = DataImpl<MergedPlane>;
  }

  namespace Extn {

    //Various types related to extinction models, to be used both in
    //configuration and extinction components. Here BC is shorthand of
    //Becker-Coppens.

    enum class BC_RecipeVersion { Std2025 = 0,
                                  Lux2025 = 1,
                                  Classic1974 = 2,
                                  Default = Std2025 };
    enum class BC_Component{ P = 0, G = 1, L = 2, F = 3,
                             Primary = P, ScndGauss = G,
                             ScndLorentz = L, ScndFresnel = F };
    enum class BC_ScndXDef { Standard,  //Standard BC X_G/X_L
                             TypeI,     //Standard, but dropping terms to enforce TypeI
                             TypeII,    //Standard, but dropping terms to enforce TypeII
                             Decoupled, //Decoupling X_G/X_L from the domain
                                        //size, reducing domain-size effects to
                                        //a point-like effect absorbed on F^2
                                        //through primary extinction
                             Default = Decoupled };

    enum class BC_ScndComponent { G = 1, L = 2, F = 3,
                                  ScndGauss = G,
                                  ScndLorentz = L,
                                  ScndFresnel = F,
                                  Gauss = G,
                                  Lorentz = L,
                                  Fresnel = F,
                                  Default = G };//fixme: decide if L is better.
    inline constexpr BC_Component as_BC_Component( BC_ScndComponent sc )
    {
      return ( sc==BC_ScndComponent::G? BC_Component::G
               : ( sc == BC_ScndComponent::L
                   ? BC_Component::L : BC_Component::F ) );
    }

  }
}

#endif

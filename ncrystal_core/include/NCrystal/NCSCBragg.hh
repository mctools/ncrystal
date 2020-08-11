#ifndef NCrystal_SCBragg_hh
#define NCrystal_SCBragg_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2020 NCrystal developers                                   //
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

#include "NCrystal/NCScatter.hh"
#include "NCrystal/NCSCOrientation.hh"

namespace NCrystal {

  class Info;
  class PlaneProvider;

  class NCRYSTAL_API SCBragg : public Scatter {
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
    SCBragg(const Info*,
            const SCOrientation&,
            double mosaicity,
            double delta_d = 0,
            PlaneProvider * plane_provider = 0,
            double prec = 1e-3, double ntrunc = 0.0 );

    //The cross-section (in barns):
    virtual double crossSection( double ekin,
                                 const double (&neutron_direction)[3] ) const;

    //There is a maximum wavelength at which Bragg diffraction is possible,
    //so ekin_low will be set to reflect this (ekin_high will be set to infinity):
    virtual void domain(double& ekin_low, double& ekin_high) const;

    //Generate scatter angle according to Bragg diffraction (defaulting to
    //isotropic if Bragg diffraction is not possible for the provided wavelength
    //and direction). This is elastic scattering and will always result in
    //delta_ekin=0:
    virtual void generateScattering( double ekin,
                                     const double (&neutron_direction)[3],
                                     double (&resulting_neutron_direction)[3],
                                     double& delta_ekin ) const ;

  private:
    virtual ~SCBragg();
    struct pimpl;
    pimpl * m_pimpl;
  };

}

#endif

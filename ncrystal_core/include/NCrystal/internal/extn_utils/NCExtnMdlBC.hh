#ifndef NCrystal_ExtnMdlBC_hh
#define NCrystal_ExtnMdlBC_hh

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

#include "NCrystal/internal/utils/NCExtraTypes.hh"
#include "NCrystal/internal/utils/NCMsg.hh"//fixme
#include "NCrystal/internal/extn_utils/NCExtnUtils.hh"
#include "NCrystal/internal/extn_utils/NCExtnBC.hh"

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////////////////////////////////////////////////////
  //
  // Implementation of models described in:
  //
  // Sabine, T. M. "The flow of radiation in a real crystal." (2006): 609-616.
  // doi: 10.1107/97809553602060000603
  //
  //////////////////////////////////////////////////////////////////////////////

  namespace Extn {

    //Any model for which y is a function of (x,sintheta) and for which it is
    //possible to evaluate x and y via template classes. The initialisation of
    //ModelData must be done in the calling code.

    template<class TXCalc,class TYCalc>
    struct GenericModel final {
      using ModelData = typename TXCalc::ModelData;

      struct NeutronData {
        typename TXCalc::NeutronXCalc xcalc;
        NeutronWavelength wavelength;
      };

      static NeutronData initNeutronData( const ModelData& md,
                                          NeutronEnergy ekin )
      {
        NeutronWavelength wl{ ekin };
        //TXCalc::NeutronXCalc foo{ md, wl };
        return NeutronData{ typename TXCalc::NeutronXCalc{ md, wl }, wl };
      }

      struct PlaneData {
        //FIXME: Delegate this to TXcalc???
        //Fields required by all models:
        double dsp;
        double fdm;
        //Other fields, as needed by this particular model:
        double fsq;//structure factor squared [barn]
        double inv2dsp; // = 1/(2*dspacing)
      };

      static PlaneData initPlaneData( const PowderBraggInput::Plane& p )
      {
        PlaneData res;
        res.dsp = p.dsp;
        res.fdm = p.fsq * p.dsp * p.mult;
        res.fsq = p.fsq;
        nc_assert( p.dsp > 0.0 );
        res.inv2dsp = 1.0 / ( 2.0 * p.dsp );
        return res;
      }

      static double extinctionFactor( const ModelData& /*used through neutron data*/,
                                      const NeutronData& n,
                                      const PlaneData& p )
      {
        const double sintheta = n.wavelength.get() * p.inv2dsp;
        const double x = n.xcalc.calcX( p, sintheta );
        return TYCalc::eval( x, sintheta );
      }
    };
  }
}

#endif

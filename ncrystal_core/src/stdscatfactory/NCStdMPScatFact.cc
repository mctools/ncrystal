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

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/interfaces/NCProcImpl.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/extd_utils/NCProcCompBldr.hh"

namespace NC = NCrystal;

//////////////////////////////////////////////////////////////////
//                                                              //
// The standard Scatter factory handling multi-phase            //
// materials. Note that multiple phases might be defined at the //
// MatCfg level, but can also appear only at the Info level     //
// (e.g. if a single NCMAT file produced multiple phases). This //
// factory deals with the second kind.                          //
//                                                              //
//////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  class StdMPScatFact : public FactImpl::ScatterFactory {
  public:
    const char * name() const noexcept final { return "stdmpscat"; }

    MultiPhaseCapability multiPhaseCapability() const override
    {
      return MultiPhaseCapability::MPOnly;
    }

    Priority query( const FactImpl::ScatterRequest& request ) const final
    {
      nc_assert_always( request.info().isMultiPhase() );
      return Priority{ 100 };
    }

    ProcImpl::ProcPtr produce( const FactImpl::ScatterRequest& request ) const final
    {
      //NB: Nested create-scatter calls in this method use
      //FactImpl::createScatter, and not globalCreateScatter typically used in
      //other plugins. This is because the individual phases in a multiphase
      //material might themselves be multiphase, meaning that we actually WANT
      //to chain the call back to ourselves (which globalCreateScatter is
      //designed to prevent).

      //Note that in order to calculate the scale factors, it is not enough to
      //use the phase-volume fractions. We must also take the differences in
      //number densities of each phase into account (since we are providing
      //cross sections per-atom, not per-volume).

      nc_assert( request.isMultiPhase() );

      Utils::ProcCompBldr scatter_phases;

      const double totnd = request.info().getNumberDensity().dbl();
      if ( ! ( totnd > 0.0 ) )
        return ProcImpl::getGlobalNullScatter();//should not be possible but just to be safe

      for ( auto&& info_ph : enumerate( request.info().getPhases() ) ) {

        //Create corresponding scatter object and add:
        const double phase_fraction = info_ph.val.first * ( info_ph.val.second->getNumberDensity().dbl() /  totnd ) ;
        if ( ! phase_fraction )
          continue;

        auto child_request = request.createChildRequest( info_ph.idx );
        scatter_phases.addfct_cl( [child_request,phase_fraction]()
        {
          ProcImpl::ProcComposition::ComponentList complist;
          complist.emplace_back( phase_fraction, FactImpl::createScatter( child_request ) );
          return complist;
        } );

      }

      //NB: When we add support for SANS physics in NCMAT files (not via
      //@CUSTOM_ sections), we will handle it here.
      return scatter_phases.finalise_scatter();
    }
  };

}

//Finally, a function which can be used to enable the above factory. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_stdmpscat_factory)()
{
  NC::FactImpl::registerFactory(std::make_unique<NC::StdMPScatFact>());
}

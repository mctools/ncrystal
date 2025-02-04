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

#include "NCrystal/internal/sanshardsphere/NCSANSSphScat.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;

//////////////////////////////////////////////////////////////////
//                                                              //
// A small test factory providing a simple SANS hard spheres    //
// model for files with @CUSTOM_HARDSPHERESANS sections.        //
//                                                              //
// This model is mainly intended to test and demonstrate the    //
// feasibility of SANS models in NCrystal.                      //
//                                                              //
//////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace {

    class SansHardSphereFact : public FactImpl::ScatterFactory {
    private:
      static constexpr auto custom_section_name = "HARDSPHERESANS";
    public:
      static constexpr auto the_factory_name = "hardspheresans";
      const char * name() const noexcept final { return the_factory_name; }

      MultiPhaseCapability multiPhaseCapability() const override
      {
        return MultiPhaseCapability::MPOnly;
      }

      Priority query( const FactImpl::ScatterRequest& request ) const final
      {
        //NB: In addition to the request.sans() check here, the
        //extractCustomDataForSANSPlugin helper function also looks for sans=0
        //entries in the phase list itself, as appropriate.
        return ( ( request.get_sans() && hasCustomDataForSANSPlugin( request.info(), custom_section_name ) )
                 ? Priority{ 200 }
                 : Priority::Unable );
      }

      ProcImpl::ProcPtr produce( const FactImpl::ScatterRequest& request ) const final
      {
        nc_assert( request.get_sans() );
        auto datas = extractCustomDataForSANSPlugin( request.info(), custom_section_name );
        nc_assert( !datas.empty() );

        ProcImpl::ProcComposition::ComponentList process_list;

        //First add all the physics provided by other code:
        process_list.emplace_back( globalCreateScatter(request) );

        //Now add our SANS model as appropriate (this might be more than one
        //instance in case of complicated material setups):
        for ( auto& data : datas ) {

          //Decode and check syntax in custom section data:
          const auto& lines = data.customData;
          double sphere_radius;
          if ( lines.size()!=1
               || lines.at(0).size()!=1
               || !safe_str2dbl(lines.at(0).at(0),sphere_radius)
               || !(sphere_radius>0.0)
               || !(sphere_radius<1e6) ) {
            NCRYSTAL_THROW2(BadInput,"Syntax error in @CUSTOM_"<< custom_section_name
                            <<" section. Expects a single positive number (the sphere radius in Angstrom).");
          }

          //Our SANS process
          process_list.emplace_back( makeSO<SANSSphereScatter>( data.scale, SANSSphereScatter::sphere_radius{ sphere_radius } ) );
        }

        return ProcImpl::ProcComposition::consumeAndCombine(std::move(process_list));
      }
    };

  }
}

//Finally, a function which can be used to enable the above factory. Note that
//this function is forward declared elsewhere or might be dynamically invoked
//(hence the C-mangling), and its name should not be changed just here:

extern "C" void NCRYSTAL_APPLY_C_NAMESPACE(register_experimentalscatfact)()
{
  NC::FactImpl::registerFactory(std::make_unique<NC::SansHardSphereFact>());
}

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

#include "NCrystal/internal/dynapi/NCDynAPI_Type1_v1.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCDynAPIUtils.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace DynAPI {

    class Type1_v1_Impl final : public ::NCrystalDynamicAPI::DynAPI_Type1_v1 {
    public:
      using PubScatterProcess = ::NCrystalDynamicAPI::DynAPI_Type1_v1::ScatterProcess;

      struct ScatterProcess
      {
        ScatterProcess( const char * cfgstr )
          : procptr( FactImpl::createScatter( cfgstr ) )
        {
        }
        ProcImpl::ProcPtr procptr;
      };

      const PubScatterProcess * createScatter( const char * cfgstr ) const override
      {
        return reinterpret_cast<PubScatterProcess*>( new ScatterProcess(cfgstr) );
      }

      void deallocateScatter( const PubScatterProcess * sp ) const override
      {
        delete reinterpret_cast<const ScatterProcess*>(sp);
      }

      double crossSectionUncached( const PubScatterProcess& pub_sp,
                                   double neutron_ekin_eV,
                                   double neutron_dir_ux,
                                   double neutron_dir_uy,
                                   double neutron_dir_uz ) const override
      {
        auto sp = reinterpret_cast<const ScatterProcess*>(&pub_sp);
        CachePtr dummycache;//<--- Fully MT safe, fully inefficient. To be
                            //revisited in a future api version!
        return sp->procptr->crossSection( dummycache,
                                          NeutronEnergy{ neutron_ekin_eV },
                                          NeutronDirection( neutron_dir_ux,
                                                            neutron_dir_uy,
                                                            neutron_dir_uz ) ).dbl();


      }

      void sampleScatterUncached( const PubScatterProcess& pub_sp,
                                  std::function<double()>& rng_fct,
                                  double& neutron_ekin_eV,
                                  double& neutron_dir_ux,
                                  double& neutron_dir_uy,
                                  double& neutron_dir_uz ) const override
      {
        auto sp = reinterpret_cast<const ScatterProcess*>(&pub_sp);
        CachePtr dummycache;//<--- Fully MT safe, fully inefficient. To be
                            //revisited in a future api version!
        DynAPIUtils::RNGWrapper rng( &rng_fct );
        auto out = sp->procptr->sampleScatter( dummycache, rng,
                                               NeutronEnergy{ neutron_ekin_eV },
                                               NeutronDirection( neutron_dir_ux,
                                                                 neutron_dir_uy,
                                                                 neutron_dir_uz ) );
        neutron_ekin_eV = out.ekin.dbl();
        neutron_dir_ux = out.direction[0];
        neutron_dir_uy = out.direction[1];
        neutron_dir_uz = out.direction[2];
      }
    };
  }
}

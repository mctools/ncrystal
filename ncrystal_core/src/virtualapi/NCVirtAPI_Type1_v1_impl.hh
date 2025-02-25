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

#include "NCrystal/virtualapi/NCVirtAPI_Type1_v1.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCVirtAPIUtils.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace VirtAPI {

    class Type1_v1_Impl final : public ::NCrystalVirtualAPI::VirtAPI_Type1_v1 {
    public:
      using PubScatterProcess = ::NCrystalVirtualAPI::VirtAPI_Type1_v1::ScatterProcess;

      struct ScatterProcess
      {
        ScatterProcess( const char * cfgstr )
          : procptr( FactImpl::createScatter( cfgstr ) ) {}
        ScatterProcess( ProcImpl::ProcPtr pp )
          : procptr( std::move(pp) ) {}
        ProcImpl::ProcPtr procptr;
      };

      const PubScatterProcess * createScatter( const char * cfgstr ) const override
      {
        return reinterpret_cast<PubScatterProcess*>( new ScatterProcess(cfgstr) );
      }

      const PubScatterProcess * cloneScatter( const PubScatterProcess * psp ) const override
      {
        return reinterpret_cast<PubScatterProcess*>
          ( new ScatterProcess
            ( reinterpret_cast<const ScatterProcess*>( psp )->procptr) );
      }

      void deallocateScatter( const PubScatterProcess * sp ) const override
      {
        delete reinterpret_cast<const ScatterProcess*>(sp);
      }

      double crossSectionUncached( const PubScatterProcess& pub_sp,
                                   const double* n ) const override
      {
        auto sp = reinterpret_cast<const ScatterProcess*>(&pub_sp);
        CachePtr dummycache;//<--- Fully MT safe, fully inefficient. To be
                            //revisited in a future api version!
        return sp->procptr->crossSection( dummycache,
                                          NeutronEnergy{ n[0] },
                                          NeutronDirection( n[1], n[2], n[3] )
                                          ).dbl();


      }

      void sampleScatterUncached( const PubScatterProcess& pub_sp,
                                  std::function<double()>& rng_fct,
                                  double* n ) const override
      {
        auto sp = reinterpret_cast<const ScatterProcess*>(&pub_sp);
        CachePtr dummycache;//<--- Fully MT safe, fully inefficient. To be
                            //revisited in a future api version!
        VirtAPIUtils::RNGWrapper rng( &rng_fct );
        auto out = sp->procptr->sampleScatter( dummycache, rng,
                                               NeutronEnergy{ n[0] },
                                               NeutronDirection( n[1],
                                                                 n[2],
                                                                 n[3] ) );
        n[0] = out.ekin.dbl();
        n[1] = out.direction[0];
        n[2] = out.direction[1];
        n[3] = out.direction[2];
      }
    };
  }
}

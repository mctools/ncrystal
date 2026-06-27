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

#include "NCrystal/virtualapi/NCVirtAPI_Type1_v1.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCVirtAPIUtils.hh"
#include <unordered_map>

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

      // Reuse a persistent cache per (thread, process) instead of allocating a
      // fresh one on every call. The cache is keyed by process pointer so it is
      // never shared between processes: accessCache does not check process
      // identity, and m_nHistory (the only reset trigger) is not unique per
      // process, so a single shared cache could be applied to the wrong process.
      // The cross-section and scatter paths use separate maps, matching the two
      // original CachePtr objects. thread_local keeps it thread-safe without
      // locks, and replaces the per-call allocation the original flagged as
      // "Fully MT safe, fully inefficient".
      double crossSectionUncached( const PubScatterProcess& pub_sp,
                                   const double* n ) const override
      {
        auto sp = reinterpret_cast<const ScatterProcess*>(&pub_sp);
        thread_local std::unordered_map<const ScatterProcess*, CachePtr> xs_caches;
        return sp->procptr->crossSection( xs_caches[sp],
                                          NeutronEnergy{ n[0] },
                                          NeutronDirection( n[1], n[2], n[3] )
                                          ).dbl();


      }

      void sampleScatterUncached( const PubScatterProcess& pub_sp,
                                  std::function<double()>& rng_fct,
                                  double* n ) const override
      {
        auto sp = reinterpret_cast<const ScatterProcess*>(&pub_sp);
        thread_local std::unordered_map<const ScatterProcess*, CachePtr> scat_caches;
        VirtAPIUtils::RNGWrapper rng( &rng_fct );
        auto out = sp->procptr->sampleScatter( scat_caches[sp], rng,
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

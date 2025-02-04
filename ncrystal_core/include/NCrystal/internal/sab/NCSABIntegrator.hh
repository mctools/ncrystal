#ifndef NCrystal_SABIntegrator_hh
#define NCrystal_SABIntegrator_hh

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

#include "NCrystal/internal/sab/NCSABScatterHelper.hh"
#include "NCrystal/internal/sab/NCSABExtender.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SAB {

    class SABIntegrator final : private MoveOnly {

      //Class which is used to process SABData and prepare helper classes
      //capable of providing cross-sections and/or sampling
      //scatterings.
      //
      //If a specific energy grid is not can be supplied, one will be
      //automatically determined based on an analysis of the scattering kernel
      //in question (see ncmat_doc.md for details of how to specify an energy
      //grid vector).
      //
      //If a SABExtender is not provided, a default single-target free gas
      //extender will be used.

      //Both constructors and destructors of the SABIntegrator are light-weight,
      //and it is safe and recommended to end the life of SABIntegrator after
      //the desired createXXX method has been called. SABIntegrators are not
      //MT-safe, so should not be shared between threads.
      //

    public:
      ~SABIntegrator();
      SABIntegrator( shared_obj<const SABData>,
                     const VectD* egrid = nullptr,
                     std::shared_ptr<const SABExtender> sabextender = nullptr );

      SABXSProvider createXSProvider() { SABXSProvider o; doit(&o,nullptr); return o; }
      SABSampler createSampler() { SABSampler o; doit(nullptr,&o); return o; }
      SABScatterHelper createScatterHelper()
      {
        SABScatterHelper o;
        doit(&o.xsprovider,&o.sampler,&o.specificJSONDescription);
        return o;
      }

    private:
      struct Impl;
      Pimpl<Impl> m_impl;
      void doit(SABXSProvider *, SABSampler*, Optional<std::string>* json = nullptr);
    };
  }
}

#endif

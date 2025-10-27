#ifndef NCrystal_ExtnScatter_hh
#define NCrystal_ExtnScatter_hh

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

#include "NCrystal/interfaces/NCProcImpl.hh"
#include "NCrystal/internal/extn_utils/NCExtnHelper.hh"

namespace NCRYSTAL_NAMESPACE {

  class Info;

  namespace Extn {

    template<class TModel>
    class ExtnScatter final : public ProcImpl::ScatterIsotropicMat {
    public:

      //Generic Extinction scatter model wrapping the model provided by TModel,
      //which must be suitable for being added to the ExtnHelper.

      template<typename ...Args>
      static shared_obj<ExtnScatter> createSO( const PowderBraggInput::Data& data,
                                               Args&& ...args )
      {
        auto modelData = TModel::initModelData( data.cell,
                                                std::forward<Args>(args)... );
        auto helper = ExtnHelper<TModel>( std::move(modelData),
                                          data.cell,
                                          data.planes );
        return makeSO<ExtnScatter<TModel>>( std::move(helper) );
      }

      template<typename ...Args>
      static shared_obj<ExtnScatter> createSO( typename TModel::ModelData&& modelData,
                                               const PowderBraggInput::Data& data )
      {
        auto helper = ExtnHelper<TModel>( std::move(modelData),
                                          data.cell,
                                          data.planes );
        return makeSO<ExtnScatter<TModel>>( std::move(helper) );
      }

      using TExtnHelper = ExtnHelper<TModel>;
      ExtnScatter( TExtnHelper&& helper )
        : m_helper(std::move(helper))
      {
      }

      const char * name() const noexcept override { return "ExtnScatter"; }
      EnergyDomain domain() const noexcept override
      {
        return m_helper.domain();
      }

      CrossSect crossSectionIsotropic(CachePtr&, NeutronEnergy ekin ) const override
      {
        //FIXME: Use cacheptr!!
        return m_helper.xsect_NoCache( ekin );
      }

      ScatterOutcomeIsotropic sampleScatterIsotropic(CachePtr&,
                                                     RNG& rng,
                                                     NeutronEnergy ekin ) const override
      {
        //FIXME: Use cacheptr!!
        return { ekin, m_helper.sampleScatMu_NoCache( rng, ekin ) };
      }

    protected:
      //FIXME: todo Optional<std::string> specificJSONDescription() const override;
      //fixme: createMerged (only on the exact same type and cell/modeldata, and
      //       for that m_helper needs createMerged functionality. Or just ignore??
    private:
      TExtnHelper m_helper;
    };

  }

}

#endif

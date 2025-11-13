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

#include "NCrystal/internal/extn_scatter/NCExtnFactory.hh"
#include "NCrystal/internal/cfgutils/NCCfgExtn.hh"
#include "NCrystal/internal/extn_scatter/NCExtnScatter.hh"
#include "NCrystal/internal/extn_utils/NCExtnSabine.hh"

namespace NC = NCrystal;
namespace NCE = NCrystal::Extn;

namespace NCRYSTAL_NAMESPACE {
  namespace Extn {
    namespace {

      template<class TXCalc,class TYCalc, typename ...Args>
      inline NC::ProcImpl::ProcPtr
      factCreate( PowderBraggInput::Data&& data, Args&& ...args )
      {
        auto modelData = TXCalc::initModelData( data.cell,
                                                std::forward<Args>(args)... );
        return ExtnScatter<GenericModel_1Comp<TXCalc,TYCalc>>::createSO( std::move(modelData),
                                                                         std::move(data) );
      }

      template<class TXCalc,class TYCalc,
               class TXCalcScnd, class TYCalcScnd, typename ...Args>
      inline NC::ProcImpl::ProcPtr
      factCreate2Comp( PowderBraggInput::Data&& data, Length domainSize, Args&& ...args )
      {
        using GM = GenericModel_2Comp<TXCalc,TYCalc,TXCalcScnd,TYCalcScnd>;
        typename GM::ModelData md;
        md.prim = TXCalc::initModelData( data.cell, domainSize );


        md.scnd = TXCalcScnd::initModelData( data.cell,
#ifdef NCRYSTAL_BCSCNDX_ALA_SHUQI
                                             domainSize,
#endif
                                             std::forward<Args>(args)... );
        return ExtnScatter<GM>::createSO( std::move(md), std::move(data) );
      }

      template<BC_RecipeVersion RecipeVersion>
      NC::ProcImpl::ProcPtr factCreateBCRecipe( PowderBraggInput::Data&& data,
                                                ExtnCfg_BC


    }
  }
}


NC::ProcImpl::ProcPtr NCE::createIsotropicExtnProc( PowderBraggInput::Data&& data,
                                                    const Cfg::ExtnCfg& extncfg_obj )
{
    NCRYSTAL_MSG("TKTEST INIT extn");//fixme
  if ( !extncfg_obj.enabled() )
    NCRYSTAL_THROW(BadInput,"createIsotropicExtnProc called "
                   "without extinction being enabled");

  const Cfg::CfgKeyValMap ecfg_kvmap{ Cfg::Extn::accessInternalVarBuf(extncfg_obj) };

  using Cfg::Extn::Model;
  using Cfg::Extn::ExtnCfg_Base;
  using Cfg::Extn::ExtnCfg_Sabine;
  using Cfg::Extn::ExtnCfg_BC;

  auto mdl_base = ExtnCfg_Base::decode(ecfg_kvmap);
  if ( mdl_base.model == Model::Sabine ) {

    auto mdl_sabine = ExtnCfg_Sabine::decode(ecfg_kvmap);
    // if ( mdl_sabine.tilt != ExtnCfg_Sabine::Tilt::Rectangular )
    //   NCRYSTAL_THROW(BadInput,"Sabine extinction model only implemented"
    //                  " for rectangular tilt model currently");
    if ( !mdl_base.grain.has_value() )
      return ExtnScatter<SabineMdlPurePrimary>::createSO( std::move(data),
                                                          mdl_base.domainSize );

    auto& grain = mdl_base.grain.value();
    if (mdl_sabine.correlation == ExtnCfg_Sabine::Correlation::Correlated ) {
      return ExtnScatter<SabineMdlCorrelatedScnd>::createSO( std::move(data),
                                                             mdl_base.domainSize,
                                                             grain.grainSize,
                                                             grain.angularSpread
                                                             );
    }
    if (mdl_sabine.tilt == ExtnCfg_Sabine::Tilt::Rectangular ) {
      return ExtnScatter<SabineMdlUncorrelatedScnd_Rec>::createSO( std::move(data),
                                                                   mdl_base.domainSize,
                                                                   grain.grainSize,
                                                                   grain.angularSpread
                                                                   );
    } else {
      return ExtnScatter<SabineMdlUncorrelatedScnd_Tri>::createSO( std::move(data),
                                                                   mdl_base.domainSize,
                                                                   grain.grainSize,
                                                                   grain.angularSpread
                                                                   );
    }
  } else if ( mdl_base.model == Model::BC ) {
    NCRYSTAL_MSG("TKTEST INIT BC");//fixme

    auto mdl_bc = ExtnCfg_BC::decode(ecfg_kvmap);

    const bool has_primary = mdl_base.domainSize.get() > 0.0;//fixme should be .has_value()
    const bool has_secondary = mdl_base.grain.has_value();//fixme should be .has_value()
    nc_assert_always( has_primary || has_secondary );

    if ( has_primary && !has_secondary ) {
      //BC primary models, no secondary
      switch( mdl_bc.recipeVersion ) {
      case ExtnCfg_BC::RecipeVersion::Std2025:
        return factCreate<BCXCalc_P,BCEval<BC_M::P,BC_RecipeVersion::Std2025>>(std::move(data),mdl_base.domainSize);
      case ExtnCfg_BC::RecipeVersion::Lux2025:
        return factCreate<BCXCalc_P,BCEval<BC_M::P,BC_RecipeVersion::Lux2025>>(std::move(data),mdl_base.domainSize);
      case ExtnCfg_BC::RecipeVersion::Classic1974:
        return factCreate<BCXCalc_P,BCEval<BC_M::P,BC_RecipeVersion::Classic1974>>(std::move(data),mdl_base.domainSize);
      default:
        nc_assert_always(false);
      };
    } else {
      //BC secondary models.
      nc_assert_always( has_secondary && mdl_base.grain.has_value() );
      auto& gr = mdl_base.grain.value();
      if ( mdl_bc.secondaryModel != ExtnCfg_BC::SecondaryModel::Gauss )
        NCRYSTAL_THROW(BadInput,"BC model only has Gaussian secondary for now");//fixme
      if ( !has_primary ) {
        //Pure secondary:
        switch( mdl_bc.recipeVersion ) {
        case ExtnCfg_BC::RecipeVersion::Std2025:
          return factCreate<BCXCalc_G,BCEval<BC_M::G,BC_RecipeVersion::Std2025>>(std::move(data),
#ifdef NCRYSTAL_BCSCNDX_ALA_SHUQI
                                                                                 mdl_base.domainSize,
#endif
                                                                                 gr.grainSize,
                                                                                 gr.angularSpread);
        case ExtnCfg_BC::RecipeVersion::Lux2025:
          return factCreate<BCXCalc_G,BCEval<BC_M::G,BC_RecipeVersion::Lux2025>>(std::move(data),
#ifdef NCRYSTAL_BCSCNDX_ALA_SHUQI
                                                                                 mdl_base.domainSize,
#endif
                                                                                 gr.grainSize,
                                                                                 gr.angularSpread);
        case ExtnCfg_BC::RecipeVersion::Classic1974:
          return factCreate<BCXCalc_G,BCEval<BC_M::G,BC_RecipeVersion::Classic1974>>(std::move(data),
#ifdef NCRYSTAL_BCSCNDX_ALA_SHUQI
                                                                                     mdl_base.domainSize,
#endif
                                                                                     gr.grainSize,
                                                                                     gr.angularSpread);
        default:
          nc_assert_always(false);
        };
      } else {
        //Both primary and secondary
        //FIXME:
        switch( mdl_bc.recipeVersion ) {
        case ExtnCfg_BC::RecipeVersion::Std2025:
          return factCreate2Comp< BCXCalc_P, BCEval<BC_M::P,BC_RecipeVersion::Std2025>,
                                  BCXCalc_G, BCEval<BC_M::G,BC_RecipeVersion::Std2025>
                                  >(std::move(data),mdl_base.domainSize,gr.grainSize,gr.angularSpread);
        case ExtnCfg_BC::RecipeVersion::Lux2025:
          return factCreate2Comp< BCXCalc_P, BCEval<BC_M::P,BC_RecipeVersion::Lux2025>,
                                  BCXCalc_G, BCEval<BC_M::G,BC_RecipeVersion::Lux2025>
                                  >(std::move(data),mdl_base.domainSize,gr.grainSize,gr.angularSpread);
        case ExtnCfg_BC::RecipeVersion::Classic1974:
          return factCreate2Comp< BCXCalc_P, BCEval<BC_M::P,BC_RecipeVersion::Classic1974>,
                                  BCXCalc_G, BCEval<BC_M::G,BC_RecipeVersion::Classic1974>
                                  >(std::move(data),mdl_base.domainSize,gr.grainSize,gr.angularSpread);
        default:
          nc_assert_always(false);
        };
      }
    }
//     if ( mdl_bc.recipeVersion == ExtnCfg_BC::RecipeVersion::Lux2025 )
//       return factCreate<BCXCalc_P,BCEval<BC_M::P,BC_RecipeVersion::Lux2025>(mdl_base.domainSize);
//     nc_assert_always( mdl_bc.recipeVersion == ExtnCfg_BC::RecipeVersion::Lux2025 )
//       return factCreate<BCXCalc_P,BCEval<BC_M::P,BC_RecipeVersion::Lux2025>(mdl_base.domainSize);

// >::createSO( BC_std::initModelData(mdl_base.domainSize)
//                                             }

// std::move(data), mdl_base.domainSize );
//     if ( mdl_bc.recipeVersion == ExtnCfg_BC::RecipeVersion::Lux2025 )
//       return ExtnScatter<BC_lux>::createSO( std::move(data), mdl_base.domainSize );
//     nc_assert_always( mdl_bc.recipeVersion == ExtnCfg_BC::RecipeVersion::Classic1974 );
//     return ExtnScatter<BC_cls>::createSO( std::move(data), mdl_base.domainSize );
  } else {
    NCRYSTAL_THROW(BadInput,"Unsupported extinction model encountered");
  }
}

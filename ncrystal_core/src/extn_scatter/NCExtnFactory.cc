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

NC::ProcImpl::ProcPtr NCE::createIsotropicExtnProc( PowderBraggInput::Data&& data,
                                                    const Cfg::ExtnCfg& extncfg_obj )
{
  if ( !extncfg_obj.enabled() )
    NCRYSTAL_THROW(BadInput,"createIsotropicExtnProc called "
                   "without extinction being enabled");

  Cfg::CfgKeyValMap ecfg_kvmap{ Cfg::Extn::accessInternalVarBuf(extncfg_obj) };

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
                                                             MosaicityFWHM { grain.angularSpread } //fixme cast
                                                             );
    }
    if (mdl_sabine.tilt == ExtnCfg_Sabine::Tilt::Rectangular ) {
      return ExtnScatter<SabineMdlUncorrelatedScnd_Rec>::createSO( std::move(data),
                                                                   mdl_base.domainSize,
                                                                   grain.grainSize,
                                                                   MosaicityFWHM { grain.angularSpread } //fixme cast
                                                                   );
    } else {
      return ExtnScatter<SabineMdlUncorrelatedScnd_Tri>::createSO( std::move(data),
                                                                   mdl_base.domainSize,
                                                                   grain.grainSize,
                                                                   MosaicityFWHM { grain.angularSpread } //fixme cast
                                                                   );
    }
  } else if ( mdl_base.model == Model::BC ) {
    auto mdl_bc = ExtnCfg_BC::decode(ecfg_kvmap);
    if ( mdl_base.grain.has_value() )
      NCRYSTAL_THROW(BadInput,"BC model for now only supports pure primary extinction");//fixme
    using BC_cls = BCMdlPurePrimary<BC_YpParameterisation::Classic1974>;
    using BC_ucls = BCMdlPurePrimary<BC_YpParameterisation::ClassicUpdated2025>;
    using BC_lux = BCMdlPurePrimary<BC_YpParameterisation::Lux2025>;
    if ( mdl_bc.ypform == ExtnCfg_BC::YpForm::Lux2025 )
      return ExtnScatter<BC_lux>::createSO( std::move(data), mdl_base.domainSize );
    if ( mdl_bc.ypform == ExtnCfg_BC::YpForm::ClassicUpdated2025 )
      return ExtnScatter<BC_ucls>::createSO( std::move(data), mdl_base.domainSize );
    nc_assert_always( mdl_bc.ypform == ExtnCfg_BC::YpForm::Classic1974 );
    return ExtnScatter<BC_cls>::createSO( std::move(data), mdl_base.domainSize );
  } else {
    NCRYSTAL_THROW(BadInput,"Unsupported extinction model encountered");
  }
}

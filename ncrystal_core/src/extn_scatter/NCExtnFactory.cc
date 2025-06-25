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
#include "NCrystal/internal/extn_scatter/NCExtnScatterSimple.hh"
#include "NCrystal/internal/cfgutils/NCCfgExtn.hh"

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

  auto mdl_base = ExtnCfg_Base::decode(ecfg_kvmap);
  //if ( oocfg_base.model == Model::BlaBla ) {
  //  ...
  // } else
  {
    if ( mdl_base.model != Model::Sabine )
      NCRYSTAL_THROW(BadInput,"Unsupported extinction model encountered");
    auto mdl_sabine = ExtnCfg_Sabine::decode(ecfg_kvmap);
    if ( mdl_sabine.tilt != ExtnCfg_Sabine::Tilt::Rectangular )
      NCRYSTAL_THROW(BadInput,"Sabine extinction model only implemented"
                     " for rectangular tilt model currently");
    if ( mdl_base.grain.has_value() )
      NCRYSTAL_THROW(BadInput,"Sabine extinction model does not"
                     " yet implement secondary extinction");
    return makeSO<ExtnScatterSimple>( std::move(data), mdl_base.domainSize );
  }
}

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

namespace NC = NCrystal;
namespace NCE = NCrystal::Extinction;

NC::ProcImpl::ProcPtr NCE::createIsotropicExtnProc( PowderBraggInput::Data&& data,
                                                    const Cfg::ExtinctionCfg& ecfg )
{
  if ( !ecfg.enabled() )
    NCRYSTAL_THROW(BadInput,"createIsotropicExtnProc called "
                   "without extinction being enabled");

  nc_assert_always(ecfg.has_sabine());//fixme

  auto& cfg_sabine = ecfg.get_sabine();
  nc_assert_always(!cfg_sabine.generic.grain.has_value());//fixme

  return makeSO<ExtnScatterSimple>( std::move(data), cfg_sabine.generic.domainSize.value );
}

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

#include "NCrystal/core/NCTypes.hh"
#include "NCrystal/internal/cfgutils/NCExtinctionCfg.hh"
#include "NCrystal/internal/cfgutils/NCCfgTypes.hh"
#include "NCrystal/internal/utils/NCStrView.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    Optional<Length> parseLength( const StrView& sv, bool allow_fail = false )
    {
      auto pv = Cfg::units_length::parse(sv);
      if ( !pv.has_value() ) {
        if (!allow_fail)
          NCRYSTAL_THROW2(BadInput,"Invalid length value in"
                          " extinction cfg string: \""<<sv<<"\"");
        return NullOpt;
      }
      //For simplicity we ignore the returned ValDbl_ShortStrOrigRep.
      return Length{ pv.value().first * Length::angstrom };
    }
  }
}

NC::Cfg::ExtinctionCfg::ExtinctionCfg()
  : ExtinctionCfg( ExtinctionCfgData() )
{
}

NC::Cfg::ExtinctionCfg::ExtinctionCfg( const ExtinctionCfgData& indata)
{
  auto raw = StrView(indata.rawData()).trimmed();
  if ( raw.empty() )
    return;
  {
    StrView::size_type ibad
      = raw.find_first_of(Cfg::forbidden_chars_value_strreps);
    if ( ibad < raw.size() ) {
      NCRYSTAL_THROW2(BadInput,"Forbidden characters in extinction cfg"
                      " string: '"<<raw[ibad]<<"'");
    }
  }
  auto scale = parseLength(raw);
  nc_assert_always(scale.has_value());
  m_scale = scale.value();
}

NC::ExtinctionCfgData NC::Cfg::ExtinctionCfg::encode() const
{
  if (!enabled())
    return {};
  std::ostringstream ss;
  constexpr double conv_meter2angstrom = Length::angstrom / Length::meter;
  ss << fmt( m_scale.value().get() * conv_meter2angstrom );
  return { std::move(ss).str() };
}

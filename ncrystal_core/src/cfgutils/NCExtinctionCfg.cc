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
namespace NCC = NCrystal::Cfg;

namespace NCRYSTAL_NAMESPACE {
  namespace Cfg {
    namespace {

      int detail_cmp_sabine( const ExtnCfg_Sabine& a, const ExtnCfg_Sabine& b )
      {
        nc_assert_always( !a.secondary.has_value()
                          && !b.secondary.has_value()
                          && "fixme: sabine secondary cmp not implemented yet" );
        if ( a.blockSize != b.blockSize )
          return ( a.blockSize < b.blockSize ? -1 : 1 );
        // auto sv_a = a.origstr_blockSize.to_view();
        // auto sv_b = b.origstr_blockSize.to_view();
        // if ( sv_a != sv_b )
        //   return ( sv_a < sv_b ? -1 : 1 );
        return 0;
      }
    }
  }
}
// ExtinctionCfg( const ExtnCfg_Sabine& sb )
//   : m_data(sb)
// {
// }

// NCC::ExtinctionCfg& NCC::ExtinctionCfg::operator=( const ExtnCfg_Sabine& sb)
// {
//   m_data = sb;
//   return *this;
// }

// namespace NCRYSTAL_NAMESPACE {
//   namespace Cfg {
//     namespace {
//       Optional<Length> parseLength( const StrView& sv, bool allow_fail = false )
//       {
//         auto pv = Cfg::units_length::parse(sv);
//         if ( !pv.has_value() ) {
//           if (!allow_fail)
//             NCRYSTAL_THROW2(BadInput,"Invalid length value in"
//                             " extinction cfg string: \""<<sv<<"\"");
//           return NullOpt;
//         }
//         //For simplicity we ignore the returned ValDbl_ShortStrOrigRep (fixme:
//         //that's not good enough for re-encoding, we should store it!)
//         return Length{ pv.value().first * Length::angstrom };
//       }
//     }
//   }
// }

// NC::Cfg::ExtinctionCfg::ExtinctionCfg( const ExtinctionCfgData& indata)
// {
//   auto raw = StrView(indata.rawData()).trimmed();
//   if ( raw.empty() )
//     return;
//   {
//     StrView::size_type ibad
//       = raw.find_first_of(Cfg::forbidden_chars_value_strreps);
//     if ( ibad < raw.size() ) {
//       NCRYSTAL_THROW2(BadInput,"Forbidden characters in extinction cfg"
//                       " string: '"<<raw[ibad]<<"'");
//     }
//   }
//   auto scale = parseLength(raw);
//   nc_assert_always(scale.has_value());
//   m_scale = scale.value();
// }

// NC::ExtinctionCfgData NC::Cfg::ExtinctionCfg::encode() const
// {
//   if (!enabled())
//     return {};
//   std::ostringstream ss;
//   constexpr double conv_meter2angstrom = Length::meter / Length::angstrom;
//   ss << fmt( m_scale.value().get() * conv_meter2angstrom );
//   return { std::move(ss).str() };
// }

NCC::ExtinctionCfg::ExtinctionCfg( const ExtinctionCfgData&& data )
  : ExtinctionCfg( detail_from_varbuf_t(), data.detail_accessRawData() )

{
 }

NCC::ExtinctionCfg::ExtinctionCfg( StrView sv )
{
  //Fixme: For now just support the simplest sabine model of only primary
  //extinction.

  auto pv_blockSize = Cfg::units_length::parse(sv);
  if ( !pv_blockSize.has_value() )
    NCRYSTAL_THROW2(BadInput,"Invalid block size length value in"
                    " extinction cfg string: \""<<sv<<"\"");

  m_data = SabineData();
  SabineData& data = m_data.value<SabineData>();
  data.obj.blockSize = Length{ pv_blockSize.value().first * Length::angstrom };
  data.origstr_blockSize = pv_blockSize.value().second;
}

int NCC::ExtinctionCfg::cmp( const ExtinctionCfg& o ) const {

  //fixme, just considering blocksize for now for simplicity:
  if ( m_data.empty() != o.m_data.empty() )
    return ( m_data.empty() ? 1 : -1 );
  if ( m_data.empty() )
    return 0;//identical (because we only have one value type so far, fixme)
  nc_assert_always( has_sabine() && o.has_sabine() );
  return detail_cmp_sabine( get_sabine(), o.get_sabine() );
}

void NCC::ExtinctionCfg::stream( std::ostream& os ) const
{
  if ( !enabled() )
    return;//stream as empty string
  if ( has_sabine() ) {
    auto& m = m_data.value<SabineData>();
    nc_assert_always( !m.obj.secondary.has_value() && "fixme: sabine secondary streaming not implemented yet" );
    if ( !m.origstr_blockSize.empty() )
      os << m.origstr_blockSize;
    else
      os << m.obj.blockSize.as_wavelength();//as_wavelength to easily get units of Aa.
  } else {
    nc_assert_always( !enabled() );
    //disabled means empty string.
  }

}

void NCC::ExtinctionCfg::streamJSON( std::ostream& os ) const
{
  os << "FIXME-ExtinctionCfg-JSON-streaming-not-implemented";
}

NCC::ExtinctionCfg::ExtinctionCfg( detail_from_varbuf_t,
                                   const detail::VarBuf& buf )
{
  if ( buf.empty() ) {
    m_data.clear();
    return;
  }
  static_assert(std::is_trivially_destructible<data_t>::value,"");
  static_assert(std::is_trivially_copyable<data_t>::value,"");
  static_assert(detail::varbuf_calc::buf_align % alignof(data_t) == 0,"");
  memcpy( (char*)&m_data, buf.data(), sizeof(data_t) );
}

NCC::detail::VarBuf NCC::ExtinctionCfg::detail_to_varbuf( detail::VarId varid ) const
{
  static_assert(std::is_trivially_destructible<data_t>::value,"");
  static_assert(std::is_trivially_copyable<data_t>::value,"");
  static_assert(detail::varbuf_calc::buf_align % alignof(data_t) == 0,"");
  return VarBuf( reinterpret_cast<const char*>(&m_data),
                 sizeof(data_t),
                 varid );
}

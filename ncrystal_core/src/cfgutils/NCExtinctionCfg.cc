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

      using ParsedDbl = std::pair<double,ValDbl_ShortStrOrigRep>;
      using KeyValVect = SmallVector<std::pair<StrView,StrView>,4>;

      template<class TParser>
      ParsedDbl parseDbl( StrView sv, StrView context )
      {
        auto pv = TParser::parse(sv);
        if ( !pv.has_value() )
          NCRYSTAL_THROW2(BadInput,"Invalid "<<context<<" value in"
                          " extinction cfg string: \""<<sv<<"\"");
        return pv.value();
      }

      FPValHolder<Length> parseLength( StrView sv, StrView context )
      {
        auto pd = parseDbl<units_length>(sv,context);
        FPValHolder<Length> res;
        res.value = Length{ pd.first * Length::angstrom };
        res.detail_orig_str_rep = pd.second;
        return res;
      }

      FPValHolder<double> parseAngle( StrView sv, StrView context )
      {
        auto pd = parseDbl<units_angle>(sv,context);
        FPValHolder<double> res;
        res.value = pd.first;
        res.detail_orig_str_rep = pd.second;
        return res;
      }

      // struct CommonValues {
      //   ParsedDbl domainSize;
      //   struct Grain {
      //     ParsedDbl grainSize;
      //     ParsedDbl angularSpread;
      //   };
      //   Optional<Grain> grain;
      // };

      ExtnCfg_Generic parseGenericValues( const SmallVector<StrView,3>& v ) {
        nc_assert(v.size()==1||v.size()==3);
        ExtnCfg_Generic res;
        res.domainSize = parseLength( v.front(), "domain size (length)" );

        if ( v.size() == 3 ) {
          res.grain.emplace();
          res.grain.value().grainSize = parseLength( v.at(1),
                                                     "grain size (length)" );
          res.grain.value().angularSpread = parseAngle( v.at(2),
                                                        "block spread (angle)" );//fixme: fix terminology
        }
        return res;
      }

      int detail_cmp_strrep( const detail::StrOrigRep& a,
                             const detail::StrOrigRep& b )
      {
        if ( a.empty() && b.empty() )
          return 0;
        if ( a.empty() )
          return -1;
        if ( b.empty() )
          return 1;
        auto sva = a.to_view();
        auto svb = b.to_view();
        return ( sva == svb ? 0 : ( sva<svb ? -1 : 1 ) );
      }

      int detail_cmp_strrep( const FPValHolder<Length>& a,
                             const FPValHolder<Length>& b )
      {
        if ( !( a.value == b.value ) )
          return a.value < b.value ? -1 : 1;
        return detail_cmp_strrep( a.detail_orig_str_rep,
                                  b.detail_orig_str_rep );
      }

      int detail_cmp_generic( const ExtnCfg_Generic& a, const ExtnCfg_Generic& b )
      {
        nc_assert_always( !a.grain.has_value()
                          && !b.grain.has_value()
                          && "fixme: secondary cmp not implemented yet" );
        return detail_cmp_strrep( a.domainSize, b.domainSize );
      }

      int detail_cmp_sabine( const ExtnCfg_Sabine& a, const ExtnCfg_Sabine& b )
      {
        int r = detail_cmp_generic( a.generic, b.generic );
        if ( r )
          return r;
        if ( a.tilt_type != b.tilt_type )
          return a.tilt_type < b.tilt_type ? -1 : 1;
        if ( a.use_correlated_model != b.use_correlated_model )
          return a.use_correlated_model ? -1 : 1;
        return 0;
      }


      void detail_stream_generic( std::ostream& os, const ExtnCfg_Generic& m )
      {
        if ( m.domainSize.detail_orig_str_rep.empty() )
          os << fmt(m.domainSize.value.get()/Length::angstrom);
        else
          os << m.domainSize.detail_orig_str_rep;//fixme: check if orig is better? Or already checked?
        nc_assert_always(!m.grain.has_value() && "streaming secondary not implemented fixme");
      }

      ExtnCfg_Sabine initModelSabine( const ExtnCfg_Generic& generic,
                                      const KeyValVect& keyvals )
      {
        //Fixme: For now, we simply do not support any models:
        if ( !keyvals.empty() )
          NCRYSTAL_THROW2(BadInput,
                          "Error in extinction cfg for model \"sabine\". "
                          "Keyword \""<<keyvals.front().first
                          <<"\" is not supported by this model");

        if ( generic.grain.has_value() )
          NCRYSTAL_THROW(BadInput,
                         "Secondary extinction for model \"sabine\" "//fixme: do not hardwire model name everywhere.
                         "is not implemented yet.");
        //        generic.domainSize
        ExtnCfg_Sabine cfg_sabine;
        cfg_sabine.generic = generic;
        //fixme: also tilt_type and use_correlated_model;
        return cfg_sabine;
      }
    }
  }
}

NCC::ExtinctionCfg::ExtinctionCfg( const ExtinctionCfgData&& data )
  : ExtinctionCfg( detail_from_varbuf_t(), data.detail_accessRawData() )
{
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

NCC::ExtinctionCfg::ExtinctionCfg( detail_from_varbuf_t,
                                   const detail::VarBuf& buf )
{
  if ( buf.empty() ) {
    m_data.clear();
    nc_assert_always(!enabled());
    return;
  }
  static_assert(std::is_trivially_destructible<data_t>::value,"");
  static_assert(std::is_trivially_copyable<data_t>::value,"");
  static_assert(detail::varbuf_calc::buf_align % alignof(data_t) == 0,"");
  memcpy( (char*)&m_data, buf.data(), sizeof(data_t) );
  nc_assert_always(enabled());//fixme: not 100% sure here
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



NCC::ExtinctionCfg::ExtinctionCfg( StrView sv )
{
  nc_assert_always(!enabled());
  //Fixme: For now just support the simplest sabine model of only primary
  //extinction. We should probably also support an empty strview giving !enabled()

  //Parse syntax like:
  //
  //   "10um"
  //   "10um/100um/1deg"
  //   "10um/100um/1deg/mdl:sabine/corr:1"
  //
  // * Entries are separated by '/'.  An entry with a colon must be a key-value
  //   entry of the form '<key>:<value>', other entries are pure values.
  //
  // * There must be exactly 1 or 3 pure values. If 1, it is assumed to be the
  //   domain size. If 3, it is domain size followed by grain size and
  //   mosaicity. The meaning and accepted range of these parameters might be
  //   slightly different depending on the model.
  //
  // * The allowed key-value entries will depend on the model which can itself
  //   be determined by a key "mdl". If "mdl" is not set explicitly, it will
  //   default to "sabine".
  //
  // * We allow (and ignore) empty parts between '/' separators.
  //
  // * We do not allow duplicate keys, and all pure value entries must come
  //   before any key-value entries.

  auto parts = sv.splitTrimmedNoEmpty('/');
  SmallVector<StrView,3> parts_pureval;
  KeyValVect parts_keyval;
  for ( auto& e : parts ) {
    if ( e.contains(':') ) {
      //key-value (check for duplicate keys already):
      auto kv = e.splitTrimmed<2>(':');
      if ( kv.size() != 2 || kv.front().empty() || kv.back().empty() ) {
        NCRYSTAL_THROW2(BadInput,
                        "Syntax error in extinction cfg \""<<sv
                        <<"\": Bad syntax in entry \""<<e<<"\".");
        auto key = kv.front();
        auto val = kv.back();
        for ( auto& existing : parts_keyval ) {
          if ( existing.first == key )
            NCRYSTAL_THROW2(BadInput,
                            "Syntax error in extinction cfg \""<<sv
                            <<"\": Repeated key \""<<key<<"\".");
        }
        parts_keyval.emplace_back( key, val );
      }
    } else {
      //pure value:
      if ( !parts_keyval.empty() )
        NCRYSTAL_THROW2(BadInput,
                        "Syntax error in extinction cfg \""<<sv<<"\": Pure"
                        " value entries must come before key:value entries.");
      parts_pureval.push_back( e );
    }
  }

  //Handle the pure values:
  if ( !isOneOf(parts_pureval.size(),1,3) )
    NCRYSTAL_THROW2(BadInput,
                    "Syntax error in extinction cfg \""<<sv
                    <<"\": Must have exactly 1 or 3 value entries");
  auto common = parseGenericValues( parts_pureval );

  //Determine model and extract other key,val parameters:
  StrView model = "sabine";
  decltype(parts_keyval) parts_keyval_nomodel;
  for ( auto& e : parts_keyval ) {
    if ( e.first == "mdl" ) {
      model = e.second;
    } else {
      parts_keyval_nomodel.push_back(e);
    }
  }

  //Fixme: Do we perhaps want to force an explicit "mdl:" specification if any
  //key-value entries are specified??

  //Finally, it is time to initialise a model!
  if ( model == "sabine" ) {
    m_data = initModelSabine( common, parts_keyval_nomodel );


  //   SabineData& data = m_data.value<SabineData>();
  // data.obj.blockSize = Length{ pv_blockSize.value().first * Length::angstrom };
  // data.origstr_blockSize = pv_blockSize.value().second;
  } else {
    NCRYSTAL_THROW2(BadInput,
                    "Syntax error in extinction cfg \""<<sv
                    <<"\": Unknown model \""<<model<<"\".");

  }


  nc_assert_always(enabled());
}

void NCC::ExtinctionCfg::stream( std::ostream& os ) const
{

  if ( !enabled() )
    return;//stream as empty string
  if ( has_sabine() ) {
    auto& m = get_sabine();
    detail_stream_generic( os, m.generic );
    if ( m.tilt_type == ExtnCfg_Sabine::TiltType::Triangular )
      os << "/tilt:tri";//rect is default
    if (!m.use_correlated_model)
      os << "/corr:0";//1 by default
  } else {
    //disabled means empty string.
    nc_assert_always( !enabled() );
  }

}

void NCC::ExtinctionCfg::streamJSON( std::ostream& os ) const
{
  os << "FIXME-ExtinctionCfg-JSON-streaming-not-implemented";
}

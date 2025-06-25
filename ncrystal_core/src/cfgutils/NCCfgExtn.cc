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

#include "NCrystal/internal/cfgutils/NCCfgExtn.hh"

namespace NC = NCrystal;
namespace NCC = NCrystal::Cfg;
namespace NCCE = NCrystal::Cfg::Extn;

namespace NCRYSTAL_NAMESPACE {
  namespace Cfg {
    namespace Extn {
      namespace {

        //fixme: revisit. These are not usually user-visible, but will show up
        //in JSON dumps. So they should still be sensible and correspond to keys
        //in the cfg-string.

        constexpr auto key_model = StrView::make("mdl");//fixme "block_size"
        constexpr auto key_domainsize = StrView::make("domain_size");//fixme "block_size"
        constexpr auto key_grainsize = StrView::make("grain_size");//fixme: grain_size.
        constexpr auto key_angspread = StrView::make("angularspread");//fixme "block_spread"

        constexpr auto default_model = StrView::make("sabine");

        void verify_model_name( const CfgKeyValMap& data, StrView name )
        {
          constexpr auto sv_default = StrView::make("default");
          StrView mdl = data.getValue( key_model, sv_default );
          if ( mdl == name || ( mdl == sv_default && default_model == name ) )
            return;
          NCRYSTAL_THROW2(BadInput,"Invalid model name (expected "<<name<<")");
        }

        //Fixme: the following infrastructure to another (local) .hh file, so we
        //could use it for other complicated cfg-variables in the future?:

        //In this context we are parsing all values to string-views, which might
        //be different than the original view. So for full flexibility we must
        //allow to return both Str-views (useful if the string is actually
        //present in the original)

        constexpr static auto flexstr_max_len = ncconstexpr_max( ShortStrDbl::max_strlen,
                                                                 ValDbl_ShortStrOrigRep::max_strlen );
        using FlexShortStr = ShortStr<flexstr_max_len+1>;

        struct ValParserLength {
          using base_parser = units_length;
          using res_type_t = Length;
          static constexpr auto type_name = "length";
          static Length validateAndEncode( double raw_val, StrView context )
          {
            double val = raw_val * Length::angstrom;
            if ( !std::isfinite(val) || !(val>0.0) )
              NCRYSTAL_THROW2(BadInput,"Invalid "<<context<<' '<<type_name
                              <<" value in extinction cfg: \""<<raw_val<<"\"");
            return Length{ val };
          }
        };

        struct ValParserAngle {
          using base_parser = units_angle;
          using res_type_t = double;
          static constexpr auto type_name = "angle";
          static double validateAndEncode( double val, StrView context )
          {
            if ( !std::isfinite(val) || !(val>0.0) || !(val<kPi) )
              NCRYSTAL_THROW2(BadInput,"Invalid "<<context<<' '<<type_name
                              <<" value in extinction cfg: \""<<val<<"\"");
            return val;
          }
        };

        template<class TParser>
        std::pair<FlexShortStr,typename TParser::res_type_t>
        parseValImpl( StrView sv, StrView context = NullOpt )
        {
          auto pv = TParser::base_parser::parse(sv);
          if ( !pv.has_value() || !std::isfinite(pv.value().first) )
            NCRYSTAL_THROW2(BadInput,"Invalid "<<context<<' '
                            <<TParser::type_name<<" value in"
                            " cfg-string: \""<<sv<<"\"");
          typename TParser::res_type_t val
            = TParser::validateAndEncode( pv.value().first, context );
          if ( !pv.value().second.empty() ) {
            auto& vstr = pv.value().second;
            return { FlexShortStr{ vstr.data(), vstr.size() }, val };
          } else {
            auto vstr = dbl2shortstr(pv.value().first);
            return { FlexShortStr{ vstr.data(), vstr.size() }, val };
          }
        }

        template<class TParser>
        typename TParser::res_type_t parseVal( StrView sv )
        {
          constexpr auto context = StrView::make("(unnamed)");
          return parseValImpl<TParser>(sv,context).second;
        }

        template<class TParser>
        FlexShortStr parseValToStr( StrView sv, StrView context )
        {
          return parseValImpl<TParser>(sv,context).first;
        }
      }
    }
  }
}

NC::Cfg::CfgKeyValMap NCCE::decode_cfgstr( VarId varid, StrView sv )
{
  CfgKeyValMap::DecodedData res;

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
  //   default to "sabine". However, a "mdl" must always be set explicitly if
  //   any key-value entries are specified. Also, to ensure we can switch the
  //   default at a later date, if the model is selected by default, it will not
  //   get re-encoded in a normalised cfg-string, but an explicitly selected
  //   model will (even if it is the default model).
  //
  // * We allow (and ignore) empty parts between '/' separators.
  //
  // * We do not allow duplicate keys, and all pure value entries must come
  //   before any key-value entries.

  auto parts = sv.splitTrimmedNoEmpty('/');
  SmallVector<StrView,3> parts_pureval;
  SmallVector<std::pair<StrView,StrView>,4> parts_keyval;
  for ( auto& e : parts ) {
    if ( e.contains(':') ) {
      //key-value (check for duplicate keys already):
      auto kv = e.splitTrimmed<2>(':');
      if ( kv.size() != 2 || kv.front().empty() || kv.back().empty() ) {
        NCRYSTAL_THROW2(BadInput,
                        "Syntax error in extinction cfg \""<<sv
                        <<"\": Bad syntax in entry \""<<e<<"\".");
      }
      auto key = kv.front();
      auto val = kv.back();
      for ( auto& existing : parts_keyval ) {
        if ( existing.first == key )
          NCRYSTAL_THROW2(BadInput,
                          "Syntax error in extinction cfg \""<<sv
                          <<"\": Repeated key \""<<key<<"\".");
      }
      parts_keyval.emplace_back( key, val );
    } else {
      //pure value:
      if ( !parts_keyval.empty() )
        NCRYSTAL_THROW2(BadInput,
                        "Syntax error in extinction cfg \""<<sv<<"\": Pure"
                        " value entries must come before key:value entries.");
      parts_pureval.push_back( e );
    }
  }

  //Handle the pure values (all for ExtnCfg_Base):
  if ( parts_pureval.size()!=1 && parts_pureval.size()!=3 ) {
    NCRYSTAL_THROW2(BadInput,
                    "Syntax error in extinction cfg \""<<sv
                    <<"\": Must have exactly 1 or 3 value entries");
  }

  FlexShortStr valbuf_ds
    = parseValToStr<ValParserLength>( parts_pureval.at(0), "domain size (length)" );
  res.emplace_back(key_domainsize,valbuf_ds.to_view());
  FlexShortStr valbuf_gs;
  FlexShortStr valbuf_as;
  if ( parts_pureval.size() == 3 ) {
    valbuf_gs = parseValToStr<ValParserLength>( parts_pureval.at(1), "grain size (length)" );
    valbuf_as = parseValToStr<ValParserAngle>( parts_pureval.at(2), "block spread (angle)" );//fixme: fix terminology
    res.emplace_back(key_grainsize,valbuf_gs.to_view());
    res.emplace_back(key_angspread,valbuf_as.to_view());
  }
  parts_pureval.clear();

  //Now look at all the key-value pairs, separating out the model specification.
  StrView model;
  decltype(parts_keyval) parts_keyval_nomodel;

  for ( auto& e : parts_keyval ) {
    if ( e.first == "mdl" ) {
      nc_assert( !model.has_value() );//no duplicate keys
      model = e.second;
    } else {
      parts_keyval_nomodel.push_back(e);
    }
  }

  if ( model.has_value() ) {
    res.emplace_back(key_model,model);
  } else if ( !parts_keyval_nomodel.empty() ) {
    NCRYSTAL_THROW2(BadInput,
                    "Syntax error in extinction cfg \""<<sv<<"\": An entry \""
                    <<parts_keyval_nomodel.front().first
                    <<":"<<parts_keyval_nomodel.front().second<<"\" is not"
                    " allowed without explicitly specifying a model"
                      " using the \"mdl:\" keyword");
  }

  //Finally, it is time to initialise a model, which should parse all the
  //remaining key-value entries!
  if ( model == "sabine" ) {
    StrView val_tilt;
    StrView val_corr;
    for ( auto& kv : parts_keyval_nomodel ) {
      if ( kv.first == "tilt" ) {
        if ( kv.second != "rec" && kv.second != "tri" )
          NCRYSTAL_THROW2( BadInput,
                           "Syntax error in extinction cfg \""<<sv<<"\": "
                           "Value of \"tilt\" must be \"rec\" or \"tri\"." );
        val_tilt = kv.second;
      } else if ( kv.first == "corr" ) {
        if ( kv.second != "1" && kv.second != "0" )
          NCRYSTAL_THROW2( BadInput,
                           "Syntax error in extinction cfg \""<<sv<<"\": "
                           "Value of \"corr\" must be \"1\" or \"0\"." );
        val_corr = kv.second;
      } else {
        NCRYSTAL_THROW2( BadInput,
                         "Syntax error in extinction cfg \""<<sv<<"\": "
                         "Model \""<<model<<"\" does not support a \""
                         <<kv.first<<"\" parameter." );
      }
    }
    if ( val_tilt.has_value() && val_tilt != "rec" )
      res.emplace_back("tilt",val_tilt);

    if ( val_corr.has_value() && val_corr != "1" )
      res.emplace_back("corr",val_corr);

  } else if ( model.has_value() ) {
    NCRYSTAL_THROW2(BadInput,
                    "Syntax error in extinction cfg \""<<sv
                    <<"\": Unknown model \""<<model<<"\".");
  }

  nc_assert( !res.empty() );

  // Return result wrapped in a CfgKeyValMap object. Note that although the res
  // object contains StrView's to temporary objects, it is OK since the
  // CfgKeyValMap will repoint them to storage in the associated VarBuf:
  return CfgKeyValMap( res, varid );
}

void NCCE::stream_to_cfgstr( std::ostream& os, const CfgKeyValMap& data )
{
  auto nleft = data.size();
  os << data.getValue( key_domainsize );
  --nleft;
  StrView model;
  if ( nleft && data.hasValue( key_model ) ) {
    model = data.getValue( key_model );
    --nleft;
  }

  if ( nleft && data.hasValue(key_grainsize) ) {
    nc_assert_always( nleft >=2 );
    os << '/' << data.getValue( key_grainsize )
       << '/' << data.getValue( key_angspread );
    nleft -= 2;
  }
  const bool has_kv_pairs = nleft>0;
  if ( !model.has_value() ) {
    //default model: cfg can be short and sweet!
    nc_assert_always( !has_kv_pairs );
    return;
  }

  //Output the remaining key-value pairs, with model first:
  os << "/mdl:" << model;
  for ( auto& kv : data.decoded() ) {
    if ( kv.first == key_model
         || kv.first == key_domainsize
         || kv.first == key_grainsize
         || kv.first == key_angspread )
      continue;
    os << '/' << kv.first<<':'<<kv.second;
  }
}

NCCE::ExtnCfg_Base NCCE::ExtnCfg_Base::decode( const CfgKeyValMap& data )
{
  //At this point CfgKeyValMap will have already been checked for issues, so we
  //can assert that all will decode correctly.
  ExtnCfg_Base res;
  {
    res.domainSize = parseVal<ValParserLength>( data.getValue( key_domainsize ) );
  }
  if ( data.hasValue(key_grainsize) ) {
    res.grain.emplace();
    auto pv = units_length::parse(data.getValue( key_grainsize ));
    res.grain.value().grainSize = parseVal<ValParserLength>( data.getValue( key_grainsize ) );
    res.grain.value().angularSpread = parseVal<ValParserAngle>( data.getValue( key_angspread ) );
  }
  auto mdl = data.getValue( key_model, NullOpt );
  if ( !mdl.has_value() || mdl == "sabine" ) {
    res.model = Model::Sabine;
  } else {
    nc_assert_always( !mdl.has_value() );//unknown model
  }
  return res;
}

NCCE::ExtnCfg_Sabine NCCE::ExtnCfg_Sabine::decode( const CfgKeyValMap& data )
{
  verify_model_name( data, "sabine" );
  ExtnCfg_Sabine res;
  if ( data.getValue( "tilt", "rec" ) == "tri" )
    res.tilt = ExtnCfg_Sabine::Tilt::Triangular;
  if ( data.getValue( "corr", "1" ) == "0" )
    res.correlation = ExtnCfg_Sabine::Correlation::Uncorrelated;
  return res;
}

std::ostream& NCCE::operator<<(std::ostream& os, const ExtnCfg_Base& cfg )
{
  os << "ExtnCfg_Base(";
  os << cfg.domainSize;
  if ( cfg.grain.has_value() ) {
    os << '/' << cfg.grain.value().grainSize
       << '/' << cfg.grain.value().angularSpread;
  }
  nc_assert_always( cfg.model == Model::Sabine );
  os << "/mdl:sabine)";
  return os;
}

std::ostream& NCCE::operator<<(std::ostream& os, const ExtnCfg_Sabine& cfg )
{
  os << "ExtnCfg_Sabine(tilt:"
     << ( cfg.tilt == ExtnCfg_Sabine::Tilt::Rectangular ? "rec" : "tri" )
     << "/corr:"
     << ( cfg.correlation
          == ExtnCfg_Sabine::Correlation::Correlated ? '1' : '0' )
     << ')';
  return os;
}

//fixme: the json encoding encodes values like "0.1m" as the strings that they
//are. Do we need a model-aware json encoding so we can get the actual values as
//fp numbers? And perhaps also the normalised complete cfg in a string for reference.

class NC::Cfg::detail::ExtnCfgBuilder {
public:
  static ExtnCfg createFromVarBuf( VarBuf&& vb )
  {
    ExtnCfg o{ NullOpt };
    o.m_varbuf = std::move(vb);
    return o;
  }
  static VarBuf& accessVarBuf( ExtnCfg& o ) { return o.m_varbuf; }
  static const VarBuf& accessVarBuf( const ExtnCfg& o ) { return o.m_varbuf; }
};

const NC::Cfg::detail::VarBuf& NCCE::accessInternalVarBuf( const ExtnCfg& ec )
{
  return detail::ExtnCfgBuilder::accessVarBuf(ec);
}


std::string NC::Cfg::ExtnCfg::to_string() const
{
  if (m_varbuf.empty())
    return {};
  std::ostringstream ss;
  Extn::stream_to_cfgstr( ss, CfgKeyValMap{ m_varbuf } );
  return ss.str();
}

//#include "NCrystal/internal/cfgutils/NCCfgVars.hh"
NC::Cfg::ExtnCfg::ExtnCfg( const char * strval )
{
  auto sv = StrView(strval?strval:"").trimmed();
  if ( !sv.empty() )
    m_varbuf = Extn::decode_cfgstr( VarId{0}/*fixme VarId::extn*/, sv ).encoded();
}

NC::Cfg::ExtnCfg NCCE::createExtnCfgFromVarBuf( VarBuf&& vb )
{
  ExtnCfg res{ NullOpt };
  if ( !vb.empty() ) {
    //fixme nc_assert_always( vb.metaData() == VarId::extn );
    detail::ExtnCfgBuilder::accessVarBuf(res) = std::move(vb);
  }
  return res;
}

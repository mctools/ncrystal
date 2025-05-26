////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2026 NCrystal developers                                   //
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

#include "NCrystal/internal/cfgutils/NCCfgTypes.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace Cfg {

    namespace {
      struct ParsedDblValue {
        double value;
        StrView unit;
        //Original str rep must be kept in a dedicated buffer, since internal
        //whitespace must be edited away, e.g. "1.0 Aa" -> "1.0Aa".
        ValDbl_ShortStrOrigRep origStrRep;//empty if did not fit.
      };

      Optional<ParsedDblValue> unitSplit(StrView sv)
      {
        sv = sv.trimmed();
        auto val = sv.toDbl();
        if ( val.has_value() ) {
          //converts cleanly (this includes e.g. "inf" etc.), so we simply
          //assume there is no unit.
          if ( sv.size() < ValDbl_ShortStrOrigRep::max_strlen )
            return ParsedDblValue{ val.value(), StrView::make(""), { sv } };
          else
            return ParsedDblValue{ val.value(), StrView::make(""), { NullOpt } };
        }
        //Go look for a unit, simply collecting all alpha chars (and underscores)
        //from the end. We also collect a version without internal whitespace
        //(e.g. "1 cm" -> "1cm" in buf:
        auto n = sv.size();
        while ( n && ( isAlpha(sv[n-1]) || sv[n-1]=='_' ) )
          --n;
        auto sv_val = sv.substr(0,n).trimmed();
        val = sv_val.toDbl();
        if ( !val.has_value() )
          return NullOpt;
        auto sv_unit = sv.substr(n);

        auto n_val_plus_unit = sv_val.size() + sv_unit.size();
        if ( n_val_plus_unit < ValDbl_ShortStrOrigRep::max_strlen ) {
          //There is room to keep the original str rep in an
          //ValDbl_ShortStrOrigRep:
          char buf[ValDbl_ShortStrOrigRep::max_strlen];
          std::memcpy( &buf[0], sv_val.data(), sv_val.size() );
          std::memcpy( &buf[sv_val.size()], sv_unit.data(), sv_unit.size() );
          buf[n_val_plus_unit] = '\0';
          return ParsedDblValue{ val.value(), sv_unit, ValDbl_ShortStrOrigRep(&buf[0],n_val_plus_unit) };
        } else {
          return ParsedDblValue{ val.value(), sv_unit, {NullOpt} };
        }
      }
    }
  }
}


NC::Optional<std::pair<double,NC::Cfg::ValDbl_ShortStrOrigRep>> NC::Cfg::units_temperature::parse(StrView sv)
{
  auto opt_p = unitSplit(sv);
  if (!opt_p.has_value())
    return NullOpt;
  auto& p = opt_p.value();
  constexpr auto uK = StrView::make("K");
  if ( p.unit.empty() || p.unit == uK ) {
    //already in final unit, peel off any trailing "K" on p.origStrRep.
    auto osv = p.origStrRep.to_view().trimmed();
    if ( osv.endswith(uK) )
      p.origStrRep = ValDbl_ShortStrOrigRep( osv.substr(0,osv.size()-uK.size()) );
  } else if ( p.unit == StrView::make("C") ) {
    p.value += 273.15;
  } else if ( p.unit == StrView::make("F") ) {
    p.value = (273.15-(32/1.8)) + (1/1.8)*p.value;
  } else {
    return NullOpt;
  }
  return std::make_pair(p.value,p.origStrRep);
}

void NC::Cfg::units_temperature::listAvailableUnits(std::ostream& os)
{
  os << "K [default], C, F";
}
void NC::Cfg::units_length::listAvailableUnits(std::ostream& os)
{
  os << "Aa [default], nm, mu, mm, cm, m";
}
void NC::Cfg::units_angle::listAvailableUnits(std::ostream& os)
{
  os << "rad [default], deg, arcmin, arcsec";
}

NC::Optional<std::pair<double,NC::Cfg::ValDbl_ShortStrOrigRep>> NC::Cfg::units_length::parse(StrView sv)
{
  auto opt_p = unitSplit(sv);
  if (!opt_p.has_value())
    return NullOpt;
  auto& p = opt_p.value();
  constexpr auto uAa = StrView::make("Aa");
  if ( p.unit.empty() || p.unit == uAa ) {
    //already in final unit, peel off any trailing "Aa" on p.origStrRep.
    auto osv = p.origStrRep.to_view().trimmed();
    if ( osv.endswith(uAa) )
      p.origStrRep = ValDbl_ShortStrOrigRep( osv.substr(0,osv.size()-uAa.size()) );
  } else if ( p.unit == StrView::make("nm") ) {
    p.value *= 10.0;
  } else if ( p.unit == StrView::make("mu") ) {
    p.value *= 1e4;
  } else if ( p.unit == StrView::make("mm") ) {
    p.value *= 1e7;
  } else if ( p.unit == StrView::make("cm") ) {
    p.value *= 1e8;
  } else if ( p.unit == StrView::make("m") ) {
    p.value *= 1e10;
  } else {
    return NullOpt;
  }
  return std::make_pair(p.value,p.origStrRep);
}

NC::Optional<std::pair<double,NC::Cfg::ValDbl_ShortStrOrigRep>> NC::Cfg::units_purenumberonly::parse(StrView sv)
{
  auto opt_p = unitSplit(sv);
  if ( !opt_p.has_value())
    return NullOpt;
  auto& p = opt_p.value();
  if ( !p.unit.empty() )
    return NullOpt;
  return std::make_pair(p.value,p.origStrRep);
}

NC::Optional<std::pair<double,NC::Cfg::ValDbl_ShortStrOrigRep>> NC::Cfg::units_angle::parse(StrView sv)
{
  auto opt_p = unitSplit(sv);
  if (!opt_p.has_value())
    return NullOpt;
  auto& p = opt_p.value();
  constexpr auto uRad = StrView::make("rad");
  if ( p.unit.empty() || p.unit == uRad ) {
    //already in final unit, peel off any trailing "rad" on p.origStrRep.
    auto osv = p.origStrRep.to_view().trimmed();
    if ( osv.endswith(uRad) )
      p.origStrRep = ValDbl_ShortStrOrigRep( osv.substr(0,osv.size()-uRad.size()) );
  } else if ( p.unit == StrView::make("deg") ) {
    p.value *= kDeg;
  } else if ( p.unit == StrView::make("arcmin") ) {
    p.value *= kArcMin;
  } else if ( p.unit == StrView::make("arcsec") ) {
    p.value *= kArcSec;
  } else {
    return NullOpt;
  }
  return std::make_pair(p.value,p.origStrRep);
}

void NC::Cfg::standardInputStrSanityCheck( const char * parname, StrView strrep )
{
  if ( !strrep.has_value() )
    NCRYSTAL_THROW2(BadInput,"Error - StrView without value provided for parameter \""<<parname<<"\"");
  auto badchar_strrep = findForbiddenChar( strrep, forbidden_chars_value_strreps, ExtraForbidOpt::RequireSimpleASCII );
  if ( badchar_strrep.has_value() )
    NCRYSTAL_THROW2(BadInput,"Forbidden character "<<badchar_strrep.value()
                    <<" in "<<parname<<" parameter value!");
}


namespace NCRYSTAL_NAMESPACE {
  namespace Cfg {
    namespace {

      uint32_t consumeUI32( const char*& buf )
      {
        uint32_t v = *reinterpret_cast<const uint32_t*>(buf);
        buf += sizeof(uint32_t);
        return v;

      }

      StrView consumeSV( const char*& buf )
      {
        uint32_t len = consumeUI32(buf);
        StrView v{ buf, (std::size_t)len  };
        buf += len;
        return v;
      }

      using CKVM_Buf = std::vector<char>;//fixme SmallVector<char,512> gives
                                         //compilation errors (and not just the
                                         //obvious .reserve thing)...
      void streamBuf( CKVM_Buf& buf, const char * it, std::size_t n )
      {
        for ( std::size_t i = 0; i < n; ++i ) {
          char val = *it++;
          buf.push_back( val );
        }
      }

      void streamUI32( CKVM_Buf& buf, uint32_t v )
      {
        const char * b = reinterpret_cast<const char*>(&v);
        streamBuf( buf, b, sizeof(uint32_t) );
      }

      void streamSV( CKVM_Buf& buf, const StrView& sv )
      {
        auto n = sv.size();
        nc_assert_always( n < 1000000000 );
        streamUI32( buf, (uint32_t)n );
        streamBuf( buf, sv.data(), n );
      }

    }

    //NOTICE: Constructors are carefully written to ensure that the StrView's in
    //m_dec refer to data in m_enc, and also that m_dec is sorted!
    CfgKeyValMap::CfgKeyValMap( const EncodedData& d )
      : m_enc(d)
    {
      m_dec = decode( m_enc );//decoding must refer to the object on m_enc, not
                              //the original "d", since there might have been a
                              //copy of the data storage area (if VarBuf does
                              //not use remote storage).
    }

    CfgKeyValMap::CfgKeyValMap( const DecodedData& d, VarId v )
      : CfgKeyValMap( encode( d, v ) )
    {
    }

    CfgKeyValMap::DecodedData CfgKeyValMap::decode( const EncodedData& data )
    {
      DecodedData res;
      const char * buf = data.data();
      uint32_t ndata = consumeUI32(buf);
      for ( uint32_t idata = 0; idata < ndata; ++idata ) {
        auto key = consumeSV(buf);
        auto val = consumeSV(buf);
        res.emplace_back( key, val );
      }
      nc_assert_always( std::is_sorted( res.begin(), res.end() ) );//fixme _always
      return res;
    }

    CfgKeyValMap::EncodedData CfgKeyValMap::encode( const DecodedData& raw,
                                                    VarId varId )
    {
      DecodedData data = raw;
      std::stable_sort( data.begin(), data.end() );
      CKVM_Buf buf;
      buf.reserve(512);
      nc_assert_always( data.size() < 1000000000 );
      streamUI32( buf, (uint32_t)data.size() );
      for ( auto i : ncrange( data.size() ) ) {
        const auto& sv_key = data.at(i).first;
        //sanity check key and ensure uniqueness:
        if ( sv_key.empty()
             ||!sv_key.contains_only("abcdefghijklmnopqrstuvwxyz_0123456789"))
          NCRYSTAL_THROW2(BadInput, "CfgKeyValMap::encode invalid key \""
                          <<sv_key<<"\"");
        for ( auto j : ncrange(i) ) {
          if ( sv_key == data.at(j).first )
            NCRYSTAL_THROW2(BadInput, "CfgKeyValMap::encode duplicate key \""
                            <<data.at(i).first<<"\"");
        }
        //sanity check value:
        const auto& sv_val = data.at(i).second;
        if (!sv_val.isSimpleASCII( StrView::AllowTabs::No,
                                   StrView::AllowNewLine::No ))
          NCRYSTAL_THROW2(BadInput, "CfgKeyValMap::encode invalid chars in"
                          " value for key \""<<data.at(i).first<<"\"");
        //stream key and value:
        streamSV( buf, sv_key );
        streamSV( buf, sv_val );
      }
      return VarBuf( buf.data(), buf.size(), varId );
    }

    void CfgKeyValMap::stream( std::ostream& os ) const
    {
      os << "CfgKeyValMap( varid="
         << static_cast<std::underlying_type<VarId>::type>(varId());
      for ( auto& e : m_dec )
        os<<", \""<<e.first<<"\": \""<<e.second<<"\"";
      os << ")";
    }

    void CfgKeyValMap::streamJSON( std::ostream& os ) const
    {
      os << '{';
      bool first(true);
      for ( auto& e : m_dec ) {
        if ( first )
          first = false;
        else
          os << ',';
        ::NCrystal::streamJSON(os,e.first);
        os << ':';
        ::NCrystal::streamJSON(os,e.second);
      }
      os << '}';
    }

    StrView CfgKeyValMap::getValue( StrView key ) const
    {
      //datasets are normally tiny => just do a linear search.
      for ( auto& e : m_dec )
        if ( e.first == key )
          return e.second;
      NCRYSTAL_THROW2(BadInput,
                      "CfgKeyValMap::getValue invalid key \""<<key<<"\"");
    }

    StrView CfgKeyValMap::getValue( StrView key, StrView defval ) const
    {
      //datasets are normally tiny => just do a linear search.
      for ( auto& e : m_dec )
        if ( e.first == key )
          return e.second;
      return defval;
    }

    bool CfgKeyValMap::hasValue( StrView key ) const
    {
      //datasets are normally tiny => just do a linear search.
      for ( auto& e : m_dec )
        if ( e.first == key )
          return true;
      return false;
    }

  }
}

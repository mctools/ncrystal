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

#include "NCrystal/internal/cfgutils/NCCfgTypes.hh"

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
  os << "Aa [default], nm, mm, cm, m";
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
  } else if ( p.unit == StrView::make("nm") ) {
    p.value *= 10.0;
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

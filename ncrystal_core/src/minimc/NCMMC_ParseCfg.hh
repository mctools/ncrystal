#ifndef NCrystal_MMC_ParseCfg_hh
#define NCrystal_MMC_ParseCfg_hh

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

#include "NCrystal/internal/minimc/NCMMC_Basket.hh"
#include "NCrystal/internal/utils/NCSpan.hh"
#include "NCrystal/internal/utils/NCVector.hh"
#include "NCrystal/internal/utils/NCStrView.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {

    namespace parseMMCCfg {

      //Utilities for parsing simple cfg-strings for MMC Geometry and Source.

      using Tokens = SmallVector<std::pair<StrView,StrView>,8>;
      struct TokenInfo {
        StrView mainName;
        Tokens tokens;
      };
      inline TokenInfo tokenize( StrView srcstr )
      {
        //tokenize semicolon separated strings of "key[=value]" pairs into a
        //list of (key,value) entries. Any repeated keys will only have
        //their last values assigned. A key without a value (i.e. no '='
        //char) will have an invalid (as opposed to empty) value, with
        //.has_value()==false.
        TokenInfo results;
        bool first(true);
        for ( auto& e : srcstr.splitTrimmedNoEmpty(';') ) {
          auto ieq = e.find('=');
          StrView key, value;
          if ( ieq == StrView::npos ) {
            key = e;
          } else {
            key = e.substr(0,ieq).rtrimmed();
            value = e.substr(ieq+1).ltrimmed();
          }
          //The optional mainName must be the first entry, and it must not have
          //a value:
          if ( first ) {
            first = false;
            if ( !value.has_value() ) {
              results.mainName = key;
              continue;
            }
          }

          //Override existing, or append new:
          bool handled(false);
          for ( auto& r : results.tokens ) {
            if ( r.first == key ) {
              handled = true;
              r.second = value;//override
              break;
            }
          }
          if ( !handled )
            results.tokens.emplace_back( key, value );
        }
        return results;
      }

      inline void applyDefaults( Tokens& tokens,
                                 StrView default_values )
      {
        std::size_t N = tokens.size();
        auto defvals_tokeninfo = tokenize(default_values);
        nc_assert(!defvals_tokeninfo.mainName.has_value());
        for ( auto& ed : defvals_tokeninfo.tokens ) {
          bool found(false);
          for ( std::size_t i = 0; i < N; ++ i ) {
            if ( ed.first == tokens.at(i).first ) {
              found = true;
              break;
            }
          }
          if (!found)
            tokens.emplace_back(ed);
        }
      }

      inline Optional<StrView> getValue( const Tokens& tokens, StrView key )
      {
        for ( auto& e : tokens )
          if ( e.first == key )
            return e.second;
        return NullOpt;
      }

      inline bool hasValue( const Tokens& tokens, StrView key )
      {
        for ( auto& e : tokens )
          if ( e.first == key )
            return true;
        return false;
      }

      inline StrView getValue_str( const Tokens& tokens, StrView key )
      {
        auto val = getValue(tokens,key);
        if ( !val.has_value() )
          NCRYSTAL_THROW2(BadInput,"Missing required parameter \""<<key<<"\"");
        //Missing (or empty) values are not allowed:
        if ( val.value().empty() )
          NCRYSTAL_THROW2(BadInput,"Missing value for parameter \""<<key<<"\"");
        return val.value();
      }

      inline double getValue_dbl( const Tokens& tokens, StrView key )
      {
        auto val_str = getValue_str(tokens,key);
        auto val_dbl = val_str.toDbl();
        if ( !val_dbl.has_value() )
          NCRYSTAL_THROW2(BadInput,"Invalid value for parameter \""<<key<<"\"");
        if ( std::isnan(val_dbl.value()) || std::isinf(val_dbl.value()) )
          NCRYSTAL_THROW2(BadInput,"Invalid value for parameter \""<<key<<"\"");
        return val_dbl.value();
      }

      inline std::size_t getValue_sizet( const Tokens& tokens, StrView key )
      {
        double x = getValue_dbl(tokens,key);
        std::size_t res;
        if ( ! (x >= 0.0 && (res=static_cast<std::size_t>(x))==x ) )
          NCRYSTAL_THROW2(BadInput,"Invalid value for parameter \""<<key<<"\"");
        return res;
      }

      inline double getValue_weight( const Tokens& tokens, StrView key )
      {
        double w = getValue_dbl(tokens,key);
        if ( !std::isfinite(w) || !(w>0.0) )
          NCRYSTAL_THROW2(BadInput,"Invalid weight value (\""
                          <<key<<"="<<fmt(w)<<"\")");
        return w;
      }

      inline bool getValue_bool( const Tokens& tokens, StrView key )
      {
        auto val_str = getValue_str(tokens,key);
        if ( val_str == "0" )
          return false;
        if ( val_str != "1" )
          NCRYSTAL_THROW2(BadInput,"Invalid value for parameter \""<<key<<"\"");
        return true;
      }

      struct EParsed {
        NeutronEnergy energy;
        std::string title;
      };

      inline EParsed getValue_Energy( const Tokens& tokens,
                                      Optional<EParsed> def_val = NullOpt )
      {
        EParsed res;
        if ( getValue( tokens, "ekin" ).has_value() ) {
          double val = getValue_dbl(tokens,"ekin");
          res.energy = NeutronEnergy{ val };
          const char * unit = "eV";
          if ( val < 0.1 ) {
            val *= 1000.0;
            unit = "meV";
          } else if ( val > 100 ) {
            val *= 0.001;
            unit = "keV";
          }
          res.title = dbl2shortstr( val ).to_string();
          res.title += unit;
        } else if ( getValue( tokens, "wl" ).has_value() ) {
          double val = getValue_dbl(tokens,"wl");
          res.energy = NeutronWavelength{ val };
          res.title = dbl2shortstr( val ).to_string();
          res.title += "Aa";
        } else {
          if ( !def_val.has_value() )
            NCRYSTAL_THROW2(BadInput,"Missing energy value (set in eV or "
                            "angstrom with \"ekin\" and \"wl\" parameters "
                            "respectively).");
          res = def_val.value();
        }
        return res;
      }

      inline void checkNoUnknown( const Tokens& tokens,
                                  const char * raw_all_accepted,
                                  const char * type_name )
      {
        if ( tokens.empty() )
          return;
        auto all_accepted = StrView(raw_all_accepted).splitTrimmedNoEmpty(';');
        auto it = tokens.begin();
        auto itE = tokens.end();
        for ( ; it != itE; ++it ) {
          bool found = false;
          for ( auto& ea : all_accepted ) {
            if ( ea == it->first ) {
              found = true;
              break;
            }
          }
          if ( !found )
            NCRYSTAL_THROW2(BadInput,"Invalid parameter for chosen "
                            <<type_name<<": \""<<it->first<<"\"");
        }
      }
    }
  }
}

#endif

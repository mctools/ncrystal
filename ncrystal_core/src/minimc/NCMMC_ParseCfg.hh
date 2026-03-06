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

      inline StrView getValue_str_allowempty( const Tokens& tokens, StrView key )
      {
        return getValue(tokens,key).value_or(StrView(""));
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

      struct FlexRangeValue {
        //Structure representing string values for positive quantities like
        // energy or wavelength:
        //  "<singleval>" (Mode::Fixed)
        //  "<val1>-<val2>" (Mode::UniformRange)
        //  "<val1>+-<val2>" (Mode::LogNormal)
        enum class Mode { Fixed, UniformRange, LogNormal };
        Mode mode = Mode::Fixed;
        double value;
        Optional<double> secondary_value;
        SmallVector<std::string,2> strforms;
        void toString( std::ostream& os ) const {
          if ( mode == Mode::Fixed ) {
            nc_assert(strforms.size()==1);
            os << strforms.at(0);
          } else {
            nc_assert(strforms.size()==2);
            os << strforms.at(0)
               << ( mode==Mode::UniformRange ? "-" : "+-" )
               << strforms.at(1);
          }
        }
        void toStringWithUnit( std::ostream& os, StrView unit ) const {
          if ( mode == Mode::Fixed ) {
            toString(os);
            os << unit;
            return;
          }
          os << '(';
          toString(os);
          os << ')' << unit;
        }

        void toJSON( std::ostream& os ) const {
          os << "{\"mode\":";
          if ( mode == Mode::Fixed ) {
            nc_assert(strforms.size()==1);
            nc_assert(!secondary_value.has_value());
            os << "\"fixed\",\"value\":";
            streamJSON(os,value);
          } else {
            nc_assert(strforms.size()==2);
            nc_assert(secondary_value.has_value());
            if ( mode == Mode::UniformRange ) {
              os << "\"uniformrange\",\"range_low\":";
              streamJSON(os,value);
              os << ",\"range_high\":";
            } else {
              nc_assert_always(mode == Mode::LogNormal);
              os << "\"lognormal\",\"mean\":";
              streamJSON(os,value);
              os << ",\"rms\":";
            }
            streamJSON(os,secondary_value.value());
          }
          os << '}';
        }
      };

      inline FlexRangeValue parse_flexrange( StrView val_str, StrView key,
                                             const char * energy_helpstr )
      {
        auto bad = [&key,energy_helpstr]()
        {
          NCRYSTAL_THROW2(BadInput,
                          "Invalid value for parameter \""<<key<<"\". "
                          <<energy_helpstr);
        };
        FlexRangeValue res;
        SmallVector<StrView,2> parts;
        if ( val_str.contains("+-") ) {
          res.mode = FlexRangeValue::Mode::LogNormal;
          auto i = val_str.find_first_of("+-");
          parts.emplace_back( val_str.substr(0,i) );
          nc_assert( val_str.substr(i).startswith("+-") );
          parts.emplace_back( val_str.substr(i+2) );
        } else if ( val_str.contains('-') ) {
          res.mode = FlexRangeValue::Mode::UniformRange;
          parts = val_str.splitTrimmed<2>('-');
        }
        auto getPosVal = [&bad](StrView vstr)
        {
          auto v = vstr.toDbl();
          if (!v.has_value()||!std::isfinite(v.value())||!(v.value()>0))
            bad();
          return v.value();
        };

        auto add_strform = [&res](double v, StrView sv)
        {
          std::string sv2;
          {
            std::stringstream ss;
            ss << fmt(v);
            sv2 = ss.str();
          }
          if ( sv.size() < sv2.size() )
            res.strforms.emplace_back( sv.to_string() );
          else
            res.strforms.emplace_back( std::move(sv2) );
        };

        if ( res.mode == FlexRangeValue::Mode::Fixed ) {
          nc_assert( parts.empty() );
          res.value = getPosVal( val_str );
          add_strform( res.value, val_str );
        } else {
          if ( parts.size() != 2 )
            bad();
          res.value = getPosVal( parts.at(0) );
          res.secondary_value = getPosVal( parts.at(1) );
          add_strform( res.value, parts.at(0) );
          add_strform( res.secondary_value.value(), parts.at(1) );
        }

        if ( res.mode == FlexRangeValue::Mode::UniformRange ) {
          nc_assert( res.secondary_value.has_value() );
          if ( ! ( res.value < res.secondary_value.value() ) )
            bad();
        }

        return res;
      }

      struct EParsed {
        //First two variables, that if set overrides all else:
        Optional<NeutronEnergy> fixed_ekin;//fast way, set for usual code path
        Optional<Temperature> maxwell;//maxwell temperature.
        //for fixed, uniformrange, lognormal:
        enum class Mode { Energy, Wavelength };
        struct FRData {
          Mode mode = Mode::Energy;
          FlexRangeValue fr;
        };
        Optional<FRData> flexrange;
        //general:
        PairDD cachevals = PairDD{0.0,0.0};
        std::string description;
      };

      inline EParsed getValue_Energy( const Tokens& tokens,
                                      Optional<StrView> defval_strwl = NullOpt )
      {
        Optional<StrView> str_ekin = getValue( tokens, "ekin" );
        Optional<StrView> str_wl = getValue( tokens, "wl" );
        if ( !str_ekin.has_value() && !str_wl.has_value() )
          str_wl = defval_strwl;

        EParsed res;
        constexpr static auto energy_helpstr =
          "Examples for how to set:"
          " \"ekin=0.025\" (fixed in eV),"
          " \"wl=1.8\" (fixed in Aa),"
          " \"ekin=0.025-0.05\" (uniform range in eV),"
          " \"wl=1.8-5\" (uniform range in Aa),"
          " \"ekin=0.025+-0.001\" (log-normal of given mean and rms in eV),"
          " \"wl=1.8+-0.01\" (log-normal of given mean and rms in Aa),"
          " \"ekin=thermal:77\" (thermal Maxwell of given temperature in K).";
        if ( !str_ekin.has_value() && !str_wl.has_value() )
          NCRYSTAL_THROW2(BadInput,"Missing energy value. "<<energy_helpstr);
        if ( str_ekin.has_value() ) {
          auto& s_ekin = str_ekin.value();
          if ( s_ekin == "help" )
            NCRYSTAL_THROW(BadInput,energy_helpstr);
          if ( str_wl.has_value() )
            NCRYSTAL_THROW2(BadInput,
                            "Do not set both \"ekin\" and \"wl\" parameters.");
          if ( s_ekin.startswith("thermal") ) {
            //Special hook for Maxwell ("thermal:300K"):
            s_ekin = s_ekin.substr(7).trimmed();
            Optional<double> tempval;
            if (s_ekin.startswith(':')&&s_ekin.endswith('K'))
              tempval = s_ekin.substr(1,s_ekin.size()-2).trimmed().toDbl();
            if ( !tempval.has_value()
                 || !(tempval.value()>1e-99)
                 || !(tempval.value()<1e99) )
              NCRYSTAL_THROW(BadInput,"Invalid format for Maxwell distribu"
                             "tion. Use a format like \"ekin=thermal:300K\".");
            res.maxwell = Temperature{ tempval.value() };
            res.cachevals.first = res.maxwell.value().kT() * 0.5;//kT/2
            std::ostringstream ss;
            ss << "Maxwell("<<fmt(tempval.value())<<"K)";
            res.description = ss.str();
          } else {
            res.flexrange.emplace();
            auto& fr = res.flexrange.value();
            fr.mode = EParsed::Mode::Energy;
            fr.fr = parse_flexrange( s_ekin, "ekin", energy_helpstr );
            std::ostringstream ss;
            fr.fr.toStringWithUnit(ss,"eV");
            res.description = ss.str();
          }
        } else {
          nc_assert_always( str_wl.has_value() );
          if ( str_wl.value() == "help" )
            NCRYSTAL_THROW(BadInput,energy_helpstr);
          if ( str_wl.value().startswith("thermal") ) {
            NCRYSTAL_THROW(BadInput,"For a Maxwell distribution, do not use"
                           " the \"wl\" parameter (use \"ekin\" instead,"
                           " like \"ekin=thermal:300K\").");
          }
          res.flexrange.emplace();
          auto& fr = res.flexrange.value();
          fr.mode = EParsed::Mode::Wavelength;
          fr.fr = parse_flexrange( str_wl.value(), "wl", energy_helpstr );
          std::ostringstream ss;
          fr.fr.toStringWithUnit(ss,"Aa");
          res.description = ss.str();
        }

        if ( !res.maxwell.has_value() ) {
          nc_assert_always(res.flexrange.has_value());
          auto& fr = res.flexrange.value();
          //cache fixed_ekin:
          if ( fr.fr.mode == FlexRangeValue::Mode::Fixed ) {
            if ( fr.mode == EParsed::Mode::Wavelength )
              res.fixed_ekin
                = NeutronWavelength( fr.fr.value ).energy();
            else
              res.fixed_ekin = NeutronEnergy( fr.fr.value );
          }
          //cache lognormal parameters:
          if ( fr.fr.mode == FlexRangeValue::Mode::LogNormal ) {
            //convert desired mu/sigma to mu/sigma of the inner gaussian, from
            //which values are sampled before being fed to exp(..):
            const double mu = fr.fr.value;
            const double sigma = fr.fr.secondary_value.value();
            nc_assert_always( ncsquare(mu) > 0 );
            nc_assert_always( mu > 0 );
            nc_assert_always( sigma > 0 );
            const double tmp = 1+ncsquare(sigma)/ncsquare(mu);
            const double mu_n = std::log( mu / std::sqrt(tmp) );
            const double sigma_n = std::sqrt(std::log(tmp));
            res.cachevals = PairDD{ mu_n, sigma_n };
          }
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

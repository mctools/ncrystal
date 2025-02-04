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

#include "NCrystal/internal/cfgutils/NCCfgVars.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace Cfg {
    struct PseudoVar {
      StrView name;
      StrView description;
      constexpr PseudoVar( const StrView& name_, const StrView& descr_ )
        : name(name_), description(descr_)
      {
      }
    };
    //NB: We had by mistake a TopLvlVar struct here and a different TopLvlVar in
    //NCCfgManip.hh, so we renamed the one here to TopLvlVarDef
    struct TopLvlVarDef {
      StrView name;
      StrView description;
      Optional<StrView> allowedUnits;
      TopLvlVarDef( const StrView& name_,
                    const StrView& descr_,
                    Optional<StrView> au_ = NullOpt )
        : name(name_), description(descr_), allowedUnits(au_)
      {
      }
    };

    SmallVector<TopLvlVarDef,3> getTopLvlVars( )
    {
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!! IMPORTANT                            !!!
      // !!! Any modification to pseudo or toplvl !!!
      // !!! parameter docs must be reflected in  !!!
      // !!! parsing code in NCCfgManip.cc        !!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      static constexpr auto sv_density = StrView::make("density");
      static constexpr auto sv_density_units = StrView::make("gcm3 kgm3 perAa3 x");
      static constexpr auto sv_density_descr
        = StrView::make("Modify the density state, which can be a scale factor (specified"
                        " with the unit \"x\"), or an absolute value (using units \"gcm3\" for"
                        " g/cm^3, \"kgm3\" for kg/m^3, or \"perAa3\" for atoms/angstrom^3)."
                        " When an absolute value is specified, that value is simply"
                        " used. However, when a scale factor is specified"
                        " (e.g. density=1.2x), then the previous value is instead scaled by"
                        " that value. Thus, appending \";density=1.2x\" to a cfg-string will"
                        " always increase the resulting material density by 20%. If"
                        " unspecified, the density state will be \"1x\" (i.e. material"
                        " densities are left as they are). Note that since it could easily"
                        " lead to undesired behaviour, scale factor density assignments"
                        " are not allowed for usage when cfg strings are embedded in input"
                        " data (but absolute density values are always allowed).");
      static constexpr auto sv_phasechoice = StrView::make("phasechoice");
      static constexpr auto sv_phasechoice_descr
        = StrView::make("Specific material sub-phases can be selected by assigning an "
                        "index value to this pseudo-parameter. More precisely, the parameter "
                        "picks out child phases in LOADED materials, not at the configuration "
                        "level. This is an important distinction since a single entry at the "
                        "cfg-level might actually result in multiple phases being loaded. As an example, "
                        "one would typically expect that loading a file called \"my_sans_sample.ncmat\" "
                        "would result in a multiphase material with two phases. Specifying "
                        "\"my_sans_sample.ncmat;phasechoice=0\" would then pick out one of "
                        "these phases, and \"my_sans_sample.ncmat;phasechoice=1\" the other. "
                        "When multi-phase materials are defined recursively with some child-phases "
                        "themselves being multi-phased, the phasechoice parameter can be specified "
                        "more than once to navigate deeper into the sub-phase tree."
                        );
      SmallVector<TopLvlVarDef,3> res;
      res.emplace_back(sv_density, sv_density_descr, sv_density_units);
      res.emplace_back(sv_phasechoice, sv_phasechoice_descr);
      return res;
      //NB: Returning like this gives spurious(?) compilation warnings with gcc7:
      // return {
      //   TopLvlVarDef( sv_density, sv_density_descr, sv_density_units ),
      //   TopLvlVarDef( sv_phasechoice, sv_phasechoice_descr ),
      // };
    }

    SmallVector<PseudoVar,3> getPseudoVars( VarGroupId gr )
    {
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!! IMPORTANT                            !!!
      // !!! Any modification to pseudo or toplvl !!!
      // !!! parameter docs must be reflected in  !!!
      // !!! parsing code in NCCfgManip.cc        !!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ( gr == VarGroupId::ScatterBase ) {
        SmallVector<PseudoVar,3> v;
        static constexpr auto sv_bragg = StrView::make("bragg");
        static constexpr auto sv_bragg_descr
          = StrView::make("This is simply an alias for the \"coh_elas\" parameter (although "
                          "the name does not strictly make sense for non-crystalline solids).");
        static constexpr auto sv_elas = StrView::make("elas");
        static constexpr auto sv_elas_descr
          = StrView::make("Convenience parameter which can be used to assign values to "
                          "all of the  \"coh_elas\", \"incoh_elas\", and \"sans\" parameters "
                          "at once. Thus, \"elas=0\" is a convenient way of disabling "
                          "elastic scattering processes and is equivalent to \"coh_elas=0;incoh_elas=0;sans=0\".");
        static constexpr auto sv_comp = StrView::make("comp");
        static constexpr auto sv_comp_descr
          = StrView::make("Convenience parameter which can be used to disable everything except "
                          " the specified components. Note that this crucially does not re-enable "
                          "the listed components if they have already been disabled. Components "
                          "are listed as a comma separated list, and recognised component names "
                          "are: \"elas\", \"incoh_elas\", \"coh_elas\", \"bragg\", \"inelas\", "
                          "and \"sans\".");
        static constexpr auto sv_bkgd = StrView::make("bkgd");
        static constexpr auto sv_bkgd_descr
          = StrView::make("Obsolete parameter which can be used to disable all physics "
                          "processes except bragg diffraction. It only accepts \"bkgd=0\" "
                          "or \"bkgd=none\", and is equivalent to \"inelas=0;incoh_elas=0;sans=0\".");
        return {
          PseudoVar( sv_bkgd, sv_bkgd_descr ),
          PseudoVar( sv_bragg, sv_bragg_descr ),
          PseudoVar( sv_comp, sv_comp_descr ),
          PseudoVar( sv_elas, sv_elas_descr ),
        };
      }
      return {};
    }
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace Cfg {
    namespace{
      using detail_printvargroupentry_t = std::pair<Optional<VarGroupId>,StrView>;
      using detail_printvargrouplist_t = std::array<detail_printvargroupentry_t,5>;
      const detail_printvargrouplist_t& detail_printvargrouplist()
      {
        static detail_printvargrouplist_t vgl =
          { detail_printvargroupentry_t{ VarGroupId::Info, "Base parameters" },
            detail_printvargroupentry_t{ VarGroupId::ScatterBase, "Basic parameters related to scattering processes"},
            detail_printvargroupentry_t{ VarGroupId::ScatterExtra, "Advanced parameters related to scattering processes (single crystals)"},
            detail_printvargroupentry_t{ VarGroupId::Absorption, "Parameters related to absorption processes"},
            detail_printvargroupentry_t{ NullOpt, "Special parameters" } };
        return vgl;
      }

      void detail_dumpCfgVarListAsJSON(std::ostream& os)
      {
        auto dumpGroup = [&os]( Optional<VarGroupId> optSelectedGroup, const StrView& groupTitle )
        {
          streamJSONDictEntry( os,"group_description", groupTitle, JSONDictPos::FIRST );
          os << ',';
          streamJSON(os,"parameters");
          os << ':' <<'[';
          bool first = true;
          auto addDivider = [&os,&first]()
          {
            if (first)
              first = false;
            else
              os << ',';
          };
          if ( optSelectedGroup.has_value() ) {
            auto selectedGroup = optSelectedGroup.value();
            auto varlist_span = Span<const VarInfo>(&varlist[0],&varlist[0]+(sizeof(varlist)/sizeof(varlist[0])));
            for ( const auto& var : varlist_span ) {
              if ( var.groupId() != selectedGroup )
                continue;
              addDivider();
              streamJSONDictEntry( os,"name",var.nameSV(),JSONDictPos::FIRST);
              streamJSONDictEntry( os,"type",var.valueTypeDescr());
              if ( var.hasUnits() ) {
                std::ostringstream tmp;
                var.streamUnitsDescr(tmp);
                streamJSONDictEntry( os,"allowed_input_units",tmp.str());
              } else {
                streamJSONDictEntry( os,"allowed_input_units",json_null_t{});
              }
              if ( var.actualUnitName() != nullptr ) {
                streamJSONDictEntry( os,"unit",var.actualUnitName() );
              }
              if ( var.hasDefaultValue() ) {
                std::ostringstream tmp;
                var.streamDefaultValue(tmp);
                std::ostringstream tmpjson;
                var.streamDefaultValueJSON(tmpjson);
                os<<',';
                streamJSON(os,"default_value");
                os << ':' << tmpjson.str();
                streamJSONDictEntry( os,"default_value_str",tmp.str());
              } else {
                streamJSONDictEntry( os,"default_value",json_null_t{});
                streamJSONDictEntry( os,"default_value_str",json_null_t{});
              }
              streamJSONDictEntry( os,"description",var.description(),JSONDictPos::LAST);
            }
            for ( auto& pv : getPseudoVars( selectedGroup ) ) {
              addDivider();
              streamJSONDictEntry( os,"name",pv.name,JSONDictPos::FIRST);
              streamJSONDictEntry( os,"type", "pseudo" );
              streamJSONDictEntry( os,"description",pv.description,JSONDictPos::LAST);
            }
          } else {
            for ( auto& tv : getTopLvlVars() ) {
              addDivider();
              streamJSONDictEntry( os,"name",tv.name,JSONDictPos::FIRST);
              streamJSONDictEntry( os,"type", "special" );
              if ( tv.allowedUnits.has_value() )
                streamJSONDictEntry( os,"allowed_input_units",tv.allowedUnits.value());
              streamJSONDictEntry( os,"description",tv.description,JSONDictPos::LAST);
            }
          }
          os << "]}";
        };
        os << '[';
        const auto& grs = detail_printvargrouplist();
        for ( auto i : ncrange(grs.size()) ) {
          dumpGroup( grs[i].first,grs[i].second );
          if ( i+1 != grs.size() )
            os << ',';
        }
        os << ']';
      }

    }
  }
}

void NC::Cfg::dumpCfgVarList( std::ostream& os,
                              CfgVarListMode mode,
                              const char * line_prefix_cstr )
{
  if ( mode == CfgVarListMode::JSON ) {
    detail_dumpCfgVarListAsJSON(os);
    return;
  }
  nc_assert_always(isOneOf(mode,CfgVarListMode::TXT_SHORT,CfgVarListMode::TXT_FULL));
  const bool shortMode{ mode == CfgVarListMode::TXT_SHORT };

  StrView line_prefix(line_prefix_cstr);

  auto dumpGroup = [&os,line_prefix,&shortMode]( Optional<VarGroupId> optSelectedGroup, const StrView& groupTitle )
  {
    auto wrap_prefix_descr = line_prefix.to_string();
    constexpr auto description_title = StrView::make("    Description: ");
    wrap_prefix_descr += "                 ";
    WordWrapCfg wrapcfg;
    wrapcfg.colwidth = 80 + line_prefix.size();
    wrapcfg.prefix = wrap_prefix_descr;
    wrapcfg.initial_offset = line_prefix.size() + description_title.size();
    os << line_prefix << groupTitle << ":\n";
    bool first = true;
    auto prVarNameAndType = [&first,&line_prefix,&os] ( StrView name, StrView typedescr )
    {
      if ( first )
        first = false;
      if ( !first )
        os << '\n';
      os << line_prefix << "  "<<name<<":\n";
      os << line_prefix << "    Type: "<<typedescr<< "\n";
    };
    auto prVarDescription = [&line_prefix,&os,&description_title,&wrapcfg] ( StrView descr )
    {
      os<<line_prefix<<description_title;
      streamWrappedText(os,descr,wrapcfg );
    };

    if ( optSelectedGroup.has_value() ) {
      auto selectedGroup = optSelectedGroup.value();
      auto varlist_span = Span<const VarInfo>(&varlist[0],&varlist[0]+(sizeof(varlist)/sizeof(varlist[0])));
      for ( const auto& var : varlist_span ) {
        if ( var.groupId() != selectedGroup )
          continue;
        if ( shortMode ) {
          os << line_prefix << "  "<<var.nameSV()<<"\n";
          continue;
        }
        prVarNameAndType(var.nameSV(),var.valueTypeDescr());
        if ( var.hasUnits() ) {
          os << line_prefix <<"    Allowed input units: ";
          var.streamUnitsDescr(os);
          os<<'\n';
        }

        if ( var.hasDefaultValue() ) {
          os << line_prefix <<"    "<<"Default value: ";
          static constexpr auto sv_string = StrView::make("string");
          bool needs_quotation = ( sv_string == var.valueTypeDescr() );
          if ( needs_quotation )
            os << '"';
          var.streamDefaultValue(os);
          if ( needs_quotation )
            os << '"';
          if ( var.hasUnits() ) {
            nc_assert(var.actualUnitName()!=nullptr);
            os << " "<<var.actualUnitName();//space to make "infAa" look a bit less silly: "inf Aa"
          }
          os<<'\n';
        } else {
          os << line_prefix <<"    "<<"No default value.\n";
        }
        prVarDescription(var.description());
      }
      for ( auto& pv : getPseudoVars( selectedGroup ) ) {
        static constexpr auto sv_pseudovar = StrView::make("pseudo");
        if ( shortMode ) {
          os << line_prefix << "  "<<pv.name<<" (pseudo parameter)\n";
          continue;
        }
        prVarNameAndType(pv.name,sv_pseudovar);
        prVarDescription(pv.description);
      }
    } else {
      for ( auto& tv : getTopLvlVars() ) {
        static constexpr auto sv_specialvar = StrView::make("special");
        if ( shortMode ) {
          os << line_prefix << "  "<<tv.name<<"\n";
          continue;
        }
        prVarNameAndType(tv.name,sv_specialvar);
        if ( tv.allowedUnits.has_value() )
          os << line_prefix <<"    Allowed input units: " << tv.allowedUnits.value() << '\n';
        prVarDescription(tv.description);
      }
    }
  };

  const auto& grs = detail_printvargrouplist();
  for ( auto i : ncrange(grs.size()) ) {
    dumpGroup( grs[i].first,grs[i].second );
    if ( !shortMode && i+1 != grs.size() )
      os << line_prefix<<'\n';
  }

}

NC::Optional<NC::Cfg::VarId> NC::Cfg::varIdFromName( StrView name )
{
  static constexpr auto itB = &varlist[0];
  static constexpr auto itE = itB + (sizeof(varlist)/sizeof*varlist);
  auto it = std::lower_bound( itB, itE, name,
                              [](const VarInfo& a, const StrView& b) -> bool { return a.nameSV() < b; } );
  if ( it != itE && it->nameSV() == name )
    return static_cast<VarId>(std::distance(itB,it));
  return NullOpt;
}


std::string NC::Cfg::FactNameRequest::to_string() const {
  if (m_excluded.empty())
    return m_specific;
  //precalc size needed:
  auto i = m_specific.size();
  for ( auto& e : m_excluded )
    i += ( e.size() + (i?2:1) );
  std::string res;
  res.reserve( i );
  res += m_specific;
  for ( auto& e : m_excluded ) {
    res += ( res.empty() ? "!" : "@!" );
    res += e;
  }
  return res;
}

namespace NCRYSTAL_NAMESPACE {
  namespace Cfg {
    namespace {
      //compile-time check of sorted varlist:
      constexpr bool check_varlist_sorted( unsigned idx = 1 ) {
        return ( idx == sizeof(varlist)/sizeof(*varlist)
                 ? true
                 : ( varlist[idx-1].nameSV().constexprLessThan( varlist[idx].nameSV() ) ? check_varlist_sorted(idx+1) : false ) );
      }
      static_assert(check_varlist_sorted(),"");
    }
  }
}

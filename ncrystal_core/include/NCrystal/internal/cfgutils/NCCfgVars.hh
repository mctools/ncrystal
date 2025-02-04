#ifndef NCrystal_CfgVars_hh
#define NCrystal_CfgVars_hh

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
#include "NCrystal/internal/utils/NCAtomUtils.hh"

namespace NCRYSTAL_NAMESPACE {

  ///////////////////////////////////////////////////////////////////////////
  //                                                                       //
  // Definition of cfg-variables (for pseudo variables, see                //
  // NCCfgManip.cc).                                                       //
  //                                                                       //
  // NB: To add a new parameter xxx:                                       //
  //    1) Add a suitable vardef_xxx struct below.                         //
  //    2) Add it into the varlist array and VarId enum below.             //
  //    4) Add C++ interface methods for it in the NCCfgManip.hh file.     //
  //    5) Add it to suitable user-visible Cfg classes (e.g. NCMatCfg.hh,  //
  //       NCFactRequests.hh).                                             //
  ///////////////////////////////////////////////////////////////////////////

  namespace Cfg {

    struct vardef_temp final : public ValDbl<vardef_temp> {
      static constexpr auto name = "temp";
      static constexpr auto group = VarGroupId::Info;
      static constexpr auto description =
        "Temperature of material in Kelvin. The special value of -1.0"
        " implies 293.15K unless input data is only valid at a specific"
        " temperature, in which case that temperature is used instead.";
      static constexpr value_type default_value() { return -1.0; }
      using units = units_temperature;
      static double value_validate( double value )
      {
        if ( ! ( value == -1.0
                 || ( value >= Temperature::allowed_range.first
                      && value <= Temperature::allowed_range.second ) ) ) {
          //Fluff in the next lines is to force intel oneapi compiler to not
          //print negative zero with a minus sign.
          Temperature printval{value};
          if ( value == 0.0 && std::signbit(value) )
            printval = Temperature{0.0};
          NCRYSTAL_THROW2(BadInput,"Out of range temperature value "
                          <<printval
                          <<" provided for parameter \""<<name
                          <<"\" (valid temperatures must be in the range "
                          <<Temperature{Temperature::allowed_range.first}
                          <<" .. "
                          <<Temperature{Temperature::allowed_range.second}<<")"
                          );
        }
        return value;
      }
    };

    struct vardef_dcutoff final : public ValDbl<vardef_dcutoff> {
      static constexpr auto name = "dcutoff";
      static constexpr auto group = VarGroupId::Info;
      static constexpr auto description =
        "Crystal planes with d-spacing below this value will be ignored. The special value of"
        " 0 implies an automatic selection of this threshold. Note that for backwards compatibility -1 is treated as 0 (for now).";
      static constexpr value_type default_value() { return 0.0; }
      using units = units_length;
      static double value_validate( double value )
      {
        if ( value == -1.0 || value == 0.0 )
          return 0.0;//For backwards compatibility we map dcutoff=-1 to 0
        if ( !(value > 0.0) )
          NCRYSTAL_THROW2(BadInput,name<<" must be >=0.0");
        if ( !( value>=1e-3 && value<=1e5 ) )//nb: value==0.0 accepted above
          NCRYSTAL_THROW2(BadInput,name<<" must be 0 (for automatic selection), or in range [1e-3,1e5] (Aa)");
        return value;
      }
    };

    struct vardef_dcutoffup final : public ValDbl<vardef_dcutoffup> {
      static constexpr auto name = "dcutoffup";
      static constexpr auto group = VarGroupId::Info;
      static constexpr auto description = "Crystal planes with d-spacing above this value will be ignored.";
      static constexpr value_type default_value() { return kInfinity; }
      using units = units_length;
      static double value_validate( double value )
      {
        if ( !(value >= 0.0) )
          NCRYSTAL_THROW2(BadInput,name<<" must be >=0.0");
        return value;
      }
    };

    struct vardef_sccutoff final : public ValDbl<vardef_sccutoff> {
      static constexpr auto name = "sccutoff";
      static constexpr auto group = VarGroupId::ScatterExtra;
      static constexpr auto description =
        "Single-crystal modelling cutoff. Crystal planes with d-spacing below this"
        " value will be approximated as having infinite mosaicity (as in a powder)."
        " A value of 0 naturally disables this approximation entirely.";
      static constexpr value_type default_value() { return 0.4; }
      using units = units_length;
      static double value_validate( double value )
      {
        if ( ! (value >= 0.0) )
          NCRYSTAL_THROW2(BadInput,name<<" must be >=0.0");
        return value;
      }
    };

    struct vardef_mos final : public ValDbl<vardef_mos> {
      static constexpr auto name = "mos";
      static constexpr auto group = VarGroupId::ScatterExtra;
      static constexpr auto description =
        "Mosaic FWHM spread in mosaic single crystals."
        " When this parameter is set, the parameters dir1 and dir2 must also be provided."
        ;
      static constexpr NullOptType default_value() { return NullOpt; }//no default value!
      using units = units_angle;
      static double value_validate( double value )
      {
        if ( !( value > 0.0) || value > kPiHalf )
          NCRYSTAL_THROW2(BadInput,name<<" must be in range (0.0,pi/2]");
        return value;
      }
    };

    struct vardef_dir1 final : public ValOrientDir<vardef_dir1> {
      static constexpr auto name = "dir1";
      static constexpr auto group = VarGroupId::ScatterExtra;
      static constexpr auto description =
        "Primary orientation axis of a single crystal. This is specified by indicating the direction of given axis in both the crystal (c1,c2,c2) and lab frames (l1,l2,l3),"
        " using the format \"@crys:c1,c2,c3@lab:l1,l2,l3\". The direction in the crystal frame can alternatively"
        " be provided in HKL space (indicating the normal of a given HKL plane), by using \"@crys_hkl:\" instead of \"@crys:\": \"dir1=@crys_hkl:c1,c2,c3@lab:l1,l2,l3\"."
        " When this parameter is set, the parameters mos and dir2 must also be provided."
        ;
      static constexpr NullOptType default_value() { return NullOpt; }//no default value!
    };

    struct vardef_dir2 final : public ValOrientDir<vardef_dir2> {
      static constexpr auto name = "dir2";
      static constexpr auto group = VarGroupId::ScatterExtra;
      static constexpr auto description =
        "Secondary orientation axis of a single crystal. This is specified using the same syntax as for the dir1 parameter."
        " In general the opening angle between the dir1 and dir2 vectors must be nonzero and identical in the crystal and lab frames, but a discrepancy"
        " up to the value of the dirtol parameter is allowed. In any case, the components of the dir2 vectors parallel to the dir1 vectors are ignored."
        " When this parameter is set, the parameters mos and dir1 must also be provided."
        ;
      static constexpr NullOptType default_value() { return NullOpt; }//no default value!
    };

    struct vardef_dirtol final : public ValDbl<vardef_dirtol> {
      static constexpr auto name = "dirtol";
      static constexpr auto group = VarGroupId::ScatterExtra;
      static constexpr auto description =
        "Tolerance parameter for the secondary direction of the single crystal orientation"
        " (see the dir2 parameter description for more information)."
        " A value of 180deg can be used to easily set up a single crystal monochromator"
        " where one is only interested in the primary direction."
        " When this parameter is set, the parameters mos, dir1, and dir2 must also be provided."
        ;
      static constexpr value_type default_value() { return 1e-4; }
      using units = units_angle;
      static double value_validate( double value )
      {
        if ( !(value>0.0 && value<=kPi ) )
          NCRYSTAL_THROW2(BadInput,name<<" must be in range (0.0,pi]");
        return value;
      }
    };

    struct vardef_mosprec final : public ValDbl<vardef_mosprec> {
      static constexpr auto name = "mosprec";
      static constexpr auto group = VarGroupId::ScatterExtra;
      static constexpr auto description = "Approximate relative numerical precision in implementation of mosaic model in single crystals.";
      static constexpr value_type default_value() { return 1e-3; }
      using units = units_purenumberonly;
      static double value_validate( double value )
      {
        if ( !(value>=1e-7) || value>1e-1 )
          NCRYSTAL_THROW2(BadInput,name<<" must be in range [1e-7,1e-1]");
        return value;
      }
    };

    struct vardef_vdoslux final : public ValInt<vardef_vdoslux> {
      static constexpr auto name = "vdoslux";
      static constexpr auto group = VarGroupId::ScatterBase;
      static constexpr auto description =
        "Setting affecting \"luxury\" level when expanding phonon"
        " spectrums (VDOS) into scattering kernels."
        " This primarily impacts the granularity of the kernel and the upper neutron energy (Emax)"
        " beyond which free-gas extrapolation is used, with implication for memory usage"
        " and initialisation time. Allowed values are:"
        " 0 (Extremely crude, 100x50 grid, Emax=0.5eV, 0.1MB, 0.02s init),"
        " 1 (Crude, 200x100 grid, Emax=1eV, 0.5MB, 0.02s init),"
        " 2 (Decent, 400x200 grid, Emax=3eV, 2MB, 0.08s init),"
        " 3 (Good, 800x400 grid, Emax=5eV, 8MB, 0.2s init),"
        " 4 (Very good, 1600x800 grid, Emax=8eV, 30MB, 0.8s init),"
        " 5 (Overkill, 3200x1600 grid, Emax=12eV, 125MB, 5s init)."
        " Note that when no actual VDOS input curve is available and one is approximated from a Debye temperature,"
        " the vdoslux level actually used will be 3 less than the one specified in this parameter (but at least 0)."
        ;
      static constexpr value_type default_value() { return 3; }
      static value_type value_validate( value_type value )
      {
        if ( value < 0 || value > 5 )
          NCRYSTAL_THROW2(BadInput,name<<" must be an integral value from 0 to 5");
        return value;
      }
    };

    struct vardef_lcaxis final : public ValVector<vardef_lcaxis> {
      static constexpr auto name = "lcaxis";
      static constexpr auto group = VarGroupId::ScatterExtra;
      static constexpr auto description =
        "Symmetry axis of anisotropic layered"
        " crystals with a layout similar to pyrolytic graphite (PG). The axis"
        " must be provided in direct lattice coordinates using a format"
        " like \"0,0,1\". Specifying"
        " this parameter along with an orientation (see dir1 and dir2"
        " parameters) will result in the appropriate anisotropic single"
        " crystal scatter model being used for Bragg diffraction."
        ;
      static constexpr NullOptType default_value() { return NullOpt; }//no default value!
      static constexpr bool auto_normalise = false;
      static void extraChecks( const Vector& v )
      {
        double m2 = v.mag2();
        if ( !(m2>0.0) )
          NCRYSTAL_THROW2(BadInput,"Null vector provided for parameter \""<<name<<"\"");
        if ( std::isinf(m2) || std::isinf(v[0]) || std::isinf(v[1]) || std::isinf(v[2]) )
          NCRYSTAL_THROW2(BadInput, "Infinities or too large values specified in "<<name<<" vector");
      }
    };

    struct vardef_lcmode final : public ValInt<vardef_lcmode> {
      static constexpr auto name = "lcmode";
      static constexpr auto group = VarGroupId::ScatterExtra;
      static constexpr auto description =
        "Choose which modelling is used for layered crystals like PG"
        " (ignored unless the lcaxis, dir1, and dir2 parameters are set)."
        " The default value 0"
        " enables the recommended model, which is both fast and"
        " accurate. A positive value N triggers a very slow but simple"
        " reference model, in which N crystallite orientations"
        " are sampled internally (the model is accurate only when N"
        " is very high). A negative value -N triggers a different (and"
        " multi-thread unsafe!) model in which each crossSection call"
        " triggers a new selection of N randomly oriented"
        " crystallites."
        ;
      static constexpr value_type default_value() { return 0; }
      static value_type value_validate( value_type value )
      {
        constexpr value_type a = 4000000000ll;
        constexpr value_type b = -a;
        if ( value < b || value > a )
          NCRYSTAL_THROW2(BadInput,name<<" must be an integral value from "<<b<<" to "<<a);
        return value;
      }
    };

    struct vardef_incoh_elas final : public ValBool<vardef_incoh_elas> {
      static constexpr auto name = "incoh_elas";
      static constexpr auto group = VarGroupId::ScatterBase;
      static constexpr auto description =
        "If enabled, incoherent elastic scattering"
        " components will be included for solid materials."
        ;
      static constexpr value_type default_value() { return true; }
    };

    struct vardef_coh_elas final : public ValBool<vardef_coh_elas> {
      static constexpr auto name = "coh_elas";
      static constexpr auto group = VarGroupId::ScatterBase;
      static constexpr auto description =
        "If enabled, coherent elastic components will be"
        " included for solid materials. In the case of crystalline"
        " materials this is essentially Bragg diffraction."
        ;
      static constexpr value_type default_value() { return true; }
    };

    struct vardef_sans final : public ValBool<vardef_sans> {
      static constexpr auto name = "sans";
      static constexpr auto group = VarGroupId::ScatterBase;
      static constexpr auto description =
        "Control presence of SANS models.  Note that this parameter"
        " is primarily added to support future developments."
        ;
      static constexpr value_type default_value() { return true; }
    };

    struct vardef_inelas final : public ValStr<vardef_inelas> {
      static constexpr auto name = "inelas";
      static constexpr auto group = VarGroupId::ScatterBase;
      static constexpr auto description =
        "Influence choice of inelastic scattering models. The default value"
        " of \"auto\" leaves the choice to the code, and values of \"none\","
        " \"0\", \"false\", or \"sterile\", all disable inelastic scattering. The"
        " standard scatter plugin currently supports"
        " additional values: \"external\", \"dyninfo\", \"vdosdebye\", and"
        " \"freegas\", and internally the \"auto\" mode will simply"
        " select the first possible of those in the listed order"
        " (falling back to \"none\" when nothing is possible). Note"
        " that \"external\" is only currently supported by .nxs files."
        " The \"dyninfo\" mode will simply base modelling on"
        " whatever dynamic information is available for each element"
        " in the input data. The \"vdosdebye\" and \"freegas\" modes"
        " overrides this, and force those models for all elements if"
        " possible (thus \"inelas=freegas;elas=0\" can be used to"
        " force a pure free-gas scattering model). The \"external\""
        " mode implies usage of an externally provided cross-section"
        " curve with an isotropic-elastic scattering model."
        ;
      static constexpr value_type default_value() { return StrView::make("auto"); }
      static Variant<StrView,std::string> str2val( StrView sv )
      {
        if (sv.empty()||!sv.contains_only("abcdefghijklmnopqrstuvwxyz_0123456789"))
          NCRYSTAL_THROW2(BadInput,"invalid value specified for parameter "<<name<<": \""<<sv<<"\"");
        nc_assert(sv.has_value());
        static constexpr auto sv_0 = StrView::make("0");
        if (isOneOf(sv,"none","0","sterile","false"))
          return sv_0;
        return sv;
      }
    };

    class FactNameRequest {
    public:
      //Book-keeping class, tracking requests for a specific named factory
      //and/or exclusion of a list of named factories.
      bool hasSpecificRequest() const noexcept { return !m_specific.empty(); }
      const std::string& specificRequest() const noexcept { return m_specific; }
      bool excludes( StrView fn ) const noexcept;
      std::string to_string() const;
      //For constructing new objects (throws BadInput errors in case of syntax issues):
      FactNameRequest withAdditionalExclude( StrView factname ) const;//just "notthis"
      FactNameRequest withNoSpecificRequest() const;//clears m_specific
      class Parser;
    private:
      FactNameRequest() = default;
      friend class Parser;
      //FactNameRequest( StrView );
      std::string m_specific;
      SmallVector_IC<std::string,2> m_excluded;//sorted list of excluded
    };

    class FactNameRequest::Parser {
      //Parses strings (primarily for scatfactory/absnfactory cfg parameters)
      //into the factory name itself plus a list of excluded factories.
      //
      //Factory names can be excluded by adding them with a "!" in front of
      //their name, and multiple entries can be added by separating them with an
      //"@" sign. However, at most one non-excluded entry can appear.
    public:
      static FactNameRequest doParse( StrView sv )
      {
        FactNameRequest res;
        auto checkValidFactoryName = []( StrView svn )
        {
          nc_assert(svn.trimmed() == svn);
          bool ok = !svn.empty();
          for ( auto& e : svn )
            if ( !isAlphaNumeric(e) && e!='_' && e!='-' )
              ok=false;
          if (!ok)
            NCRYSTAL_THROW2(BadInput,"Not a valid factory name: \""<<svn<<"\"");
        };

        //e.g. "myfact", or "myfact,!notthis,!notthiseither"
        for ( auto& e : sv.splitTrimmedNoEmpty('@') ) {
          if ( e.startswith('!') ) {
            auto exlname = e.substr(1).trimmed();
            checkValidFactoryName(exlname);
            if ( !res.excludes(exlname) )
              res.m_excluded.push_back(exlname.to_string());
          } else {
            checkValidFactoryName(e);
            if ( ! res.m_specific.empty() )
              NCRYSTAL_THROW2(BadInput,"Contains more than one (non-negated) entry (\""
                              <<res.m_specific<<"\" and \""<<e<<"\".");
            res.m_specific.assign(e.data(),e.size());
          }
        }
        if ( !res.m_specific.empty() && res.excludes(res.m_specific) )
          NCRYSTAL_THROW2(BadInput,"The factory \""<<res.m_specific<<"\" is both specified"
                          " as being simultaneously required and excluded.");
        return res;
      }
    };

    struct vardef_infofactory final : public ValStr<vardef_infofactory> {
      static constexpr auto name = "infofactory";
      static constexpr auto group = VarGroupId::Info;
#ifdef ncrystal_detail_xxxfact_descr_helper
#  undef ncrystal_detail_xxxfact_descr_helper
#endif
#define ncrystal_detail_xxxfact_descr_helper(x)                         \
      "This parameter can be used by experts to bypass the usual"      \
      " factory selection logic for " x " objects."                     \
      " A factory can be selected by providing its name, or"            \
      " excluded by prefixing the name with \"!\". Multiple"            \
      " entries must be separated by an \"@\" sign (obviously at most"  \
      " one non-excluded entry can appear)."
      static constexpr auto description = ncrystal_detail_xxxfact_descr_helper("material Info");
      static constexpr value_type default_value() { return StrView::make(""); }
      static Variant<StrView,std::string> str2val( StrView sv ) {
        try {
          return FactNameRequest::Parser::doParse(sv).to_string();
        } catch (Error::BadInput&err) {
          NCRYSTAL_THROW2(BadInput,"Syntax error in "<<name<<" parameter. Error is: "<<err.what());
        }
      }
    };

    struct vardef_scatfactory final : public ValStr<vardef_scatfactory> {
      static constexpr auto name = "scatfactory";
      static constexpr auto group = VarGroupId::ScatterBase;
      static constexpr auto description = ncrystal_detail_xxxfact_descr_helper("Scatter");
      static constexpr value_type default_value() { return StrView::make(""); }
      static Variant<StrView,std::string> str2val( StrView sv ) {
        try {
          return FactNameRequest::Parser::doParse(sv).to_string();
        } catch (Error::BadInput&err) {
          NCRYSTAL_THROW2(BadInput,"Syntax error in "<<name<<" parameter. Error is: "<<err.what());
        }
      }
    };

    struct vardef_absnfactory final : public ValStr<vardef_absnfactory> {
      static constexpr auto name = "absnfactory";
      static constexpr auto group = VarGroupId::Absorption;
      static constexpr auto description = ncrystal_detail_xxxfact_descr_helper("Absorption");
      static constexpr value_type default_value() { return StrView::make(""); }
      static Variant<StrView,std::string> str2val( StrView sv ) {
        try {
          return FactNameRequest::Parser::doParse(sv).to_string();
        } catch (Error::BadInput&err) {
          NCRYSTAL_THROW2(BadInput,"Syntax error in "<<name<<" parameter. Error is: "<<err.what());
        }
      }
    };

    struct vardef_ucnmode final : public ValStr<vardef_ucnmode> {
      static constexpr auto name = "ucnmode";
      static constexpr auto group = VarGroupId::ScatterExtra;
      static constexpr auto description =
        "Modify how UCN (ultra cold neutron) production is handled in inelastic"
        " models. The value \"refine\" simply improves the modelling by replacing the usual"
        " scattering kernel treatment near the kinematic endpoint, where the neutron ends"
        " with less than 300neV, with a different model. The values \"only\" and \"remove\""
        " performs the same split of the modelling, but then leaves out either all non-UCN"
        " or all UCN processes, respectively, from the inelastic cross sections. Finally, the"
        " threshold value of 300neV can be modified by appending the desired value to the"
        " first keyword, separated by a \":\" character. The default unit is eV, but meV and"
        " neV are supported as well, so \"ucnmode=refine:200neV\", \"ucnmode=remove:2e-7eV\","
        " \"ucnmode=remove:2e-7\", and \"ucnmode=only:0.0002meV\" all specify the same threshold. In"
        " addition to simply refining the UCN model, the primary intended purpose of the"
        " ucnmode parameter is to allow one to split out the UCN process from the rest, in"
        " order to perform biased Monte Carlo simulations of UCN production in moderators.";
      static constexpr value_type default_value() { return StrView::make(""); }
      static bool isStdKeyword( const StrView& ssv)
      {
        static constexpr auto sv_refine = StrView::make("refine");
        static constexpr auto sv_remove = StrView::make("remove");
        static constexpr auto sv_only = StrView::make("only");
        return isOneOf(ssv,sv_refine,sv_remove,sv_only);
      };

      static Optional<UCNMode> decode_value( StrView sv )
      {
        static constexpr auto sv_refine = StrView::make("refine");
        static constexpr auto sv_only = StrView::make("only");
        //assumes "sv" comes from str2val below.
        if ( sv.empty() )
          return NullOpt;
        auto decodeMode = []( StrView svmode )
        {
          if ( svmode == sv_refine )
            return UCNMode::Mode::Refine;
          if ( svmode == sv_only )
            return UCNMode::Mode::Only;
          nc_assert( svmode == StrView::make("remove") );
          return UCNMode::Mode::Remove;
        };
        if ( !sv.contains(':') ) {
          UCNMode res;
          res.mode = decodeMode(sv);
          return res;
        }
        auto parts = sv.splitTrimmed<2>(':');
        nc_assert( parts.size() == 2 && isStdKeyword(parts.at(0)) );
        StrView thrstr = parts.at(1);
        // Optional<double> thrval;
        // double unit=1.0;
        auto decodeWithUnit = [thrstr]( StrView unitname, StrView unitfpprefix, double unitvalue ) -> Optional<double> {
          if ( !thrstr.endswith(unitname) )
            return NullOpt;
          auto valstr = thrstr.substr(0,thrstr.size()-unitname.size());
          if ( !unitfpprefix.empty() && !valstr.contains_any("eE") ) {
            //adding stuff like e-9 to the string before conversion to double
            //seems to avoid introducing imprecision, unlike multiplying with
            //1e-9 afterwards!
            auto tmp = valstr.to_string() + unitfpprefix.to_string();
            auto val = StrView(tmp).toDbl();
            if ( val.has_value() )
              return val;
          }
          auto val2 = valstr.toDbl();
          if ( val2.has_value() )
            val2.value() *= unitvalue;
          return val2;
        };
        auto val = decodeWithUnit("neV","e-9",1e-9);
        if ( !val.has_value() )
          val = decodeWithUnit("meV","e-3",1e-3);
        if ( !val.has_value() )
          val = decodeWithUnit("eV","",1.0);//check AFTER neV/meV since they end with "eV".
        if ( !val.has_value() )
          val = thrstr.toDbl();
        nc_assert( val.has_value() );//because we encoded it outselves!
        UCNMode res;
        res.mode = decodeMode(parts.at(0));
        res.threshold = NeutronEnergy{ val.value() };
        return res;
      }

      static Variant<StrView,std::string> str2val( StrView sv ) {
        try {
          if ( sv.empty() )
            return sv;
          if ( isStdKeyword(sv) )
            return sv;
          auto parts = sv.splitTrimmed<2>(':');
          if ( parts.size() != 2 || !isStdKeyword(parts.at(0)) )
            NCRYSTAL_THROW2(BadInput,"\""<<sv<<"\" is not in a valid format");
          StrView thr = parts.at(1);
          double unit = 1.0;
          auto orig_thr_with_units = thr.to_string();
          if ( thr.endswith("neV") ) {
            unit = 1e-9;
            thr = thr.substr(0,thr.size()-3).rtrimmed();
            orig_thr_with_units = thr.to_string() + "neV";
          } else if ( thr.endswith("meV") ) {
            unit = 1e-3;
            thr=thr.substr(0,thr.size()-3).rtrimmed();
            orig_thr_with_units = thr.to_string() + "meV";
          } else if ( thr.endswith("eV") ) {
            thr=thr.substr(0,thr.size()-2).rtrimmed();
            orig_thr_with_units = thr.to_string();
          }
          auto thr_val = thr.toDbl();
          if ( !thr_val.has_value() )
            NCRYSTAL_THROW2(BadInput,"Invalid number: "<<thr);
          double val = thr_val.value() * unit;
          if ( !(val>0.0) || val > 1e3 )
            NCRYSTAL_THROW2(BadInput,"UCN threshold out of range: "<<sv);
          //encode in best unit (same logic as in NCTypes.cc):
          std::ostringstream ss, ssthr;
          ss << parts.at(0) << ':';
          if ( val >= 1e-9 && val < 1000e-9 ) {
            ssthr << fmt(val * 1e9) << "neV";
          } else if ( val >= 1e-3 && val < 1.0 ) {
            ssthr << fmt(val * 1e3) << "meV";
          } else {
            ssthr << fmt(val);
          }
          auto proposed_thr_with_units = ssthr.str();
          if ( proposed_thr_with_units.size() < orig_thr_with_units.size() )
            ss << proposed_thr_with_units;
          else
            ss << orig_thr_with_units;
          return ss.str();
        } catch (Error::BadInput&err) {
          NCRYSTAL_THROW2(BadInput,"Syntax error in "<<name<<" parameter. Error is: "<<err.what());
        }
      }
    };

    struct vardef_atomdb final : public ValStr<vardef_atomdb> {
      static constexpr auto name = "atomdb";
      static constexpr auto group = VarGroupId::Info;
      static constexpr auto description =
        "Modify atomic definitions if supported"
        " (in practice this is unlikely to be supported by anything"
        " except NCMAT data). The string must follow a syntax"
        " identical to that used in @ATOMDB sections of NCMAT file"
        " (cf. https://github.com/mctools/ncrystal/wiki/NCMAT-format),"
        " with a few exceptions explained here:"
        " First of all, colons (':') are interpreted as whitespace characters,"
        " which might occasionally be useful (e.g. on the command line). Next,"
        " '@' characters play the role of line separators. Finally,"
        " when used with an NCMAT file that already includes an"
        " internal @ATOMDB section, the effect will essentially be to"
        " combine the two sections by appending the atomdb lines from"
        " this cfg parameter to the lines already present in the input"
        " data. The exception is the case where the cfg parameter"
        " contains an initial line with the single word \"nodefaults\""
        " the effect of which will always be the same as if it was"
        " placed on the very first line in the @ATOMDB section"
        " (i.e. NCrystal's internal database of elements and isotopes"
        " will be ignored)."
        ;

      static constexpr value_type default_value() { return StrView::make(""); }
      static Variant<StrView,std::string> str2val( StrView sv )
      {
        //Split lines on @, replace ':' with space, normalise space in each line and throw away the space there.
        std::string result;
        for ( auto& line_sv : sv.splitTrimmedNoEmpty('@') ) {
          auto line = line_sv.to_string();
          strreplace(line,":"," ");//support ':' as ' '
          auto parts = StrView(line).split();
          if ( parts.empty() )
            continue;
          auto part_joined = joinstr(parts,StrView::make(":"));
          try {
            validateAtomDBLine( split2(part_joined,0,':') );
          } catch (Error::BadInput&e) {
            NCRYSTAL_THROW2(BadInput,"Invalid entry in "<<name<<" cfg parameter in the line: \""
                            <<part_joined<<"\". Error is: "<<e.what());
          }
          if ( part_joined == "nodefaults" && !result.empty() )
            NCRYSTAL_THROW2(BadInput,"Invalid entry in "<<name<<" cfg parameter (\"nodefaults\" must be the first line).");
          if ( !result.empty() )
            result += '@';
          result.append(part_joined.data(),part_joined.size());
        }
        return result;
      }
    };

    static constexpr VarInfo varlist[] = {
      make_varinfo<vardef_absnfactory>(),
      make_varinfo<vardef_atomdb>(),
      make_varinfo<vardef_coh_elas>(),
      make_varinfo<vardef_dcutoff>(),
      make_varinfo<vardef_dcutoffup>(),
      make_varinfo<vardef_dir1>(),
      make_varinfo<vardef_dir2>(),
      make_varinfo<vardef_dirtol>(),
      make_varinfo<vardef_incoh_elas>(),
      make_varinfo<vardef_inelas>(),
      make_varinfo<vardef_infofactory>(),
      make_varinfo<vardef_lcaxis>(),
      make_varinfo<vardef_lcmode>(),
      make_varinfo<vardef_mos>(),
      make_varinfo<vardef_mosprec>(),
      make_varinfo<vardef_sans>(),
      make_varinfo<vardef_scatfactory>(),
      make_varinfo<vardef_sccutoff>(),
      make_varinfo<vardef_temp>(),
      make_varinfo<vardef_ucnmode>(),
      make_varinfo<vardef_vdoslux>()
    };

    inline constexpr unsigned constexpr_varName2Idx( const char * name )
    {
      return constexpr_name2Idx( varlist, name );
    }

    enum class detail::VarId : std::uint32_t {
      temp    = constexpr_varName2Idx("temp"),
      sccutoff = constexpr_varName2Idx("sccutoff"),
      dcutoff = constexpr_varName2Idx("dcutoff"),
      dcutoffup = constexpr_varName2Idx("dcutoffup"),
      dirtol = constexpr_varName2Idx("dirtol"),
      mosprec = constexpr_varName2Idx("mosprec"),
      vdoslux = constexpr_varName2Idx("vdoslux"),
      lcmode = constexpr_varName2Idx("lcmode"),
      lcaxis = constexpr_varName2Idx("lcaxis"),
      ucnmode = constexpr_varName2Idx("ucnmode"),
      mos = constexpr_varName2Idx("mos"),
      dir1 = constexpr_varName2Idx("dir1"),
      dir2 = constexpr_varName2Idx("dir2"),
      sans = constexpr_varName2Idx("sans"),
      incoh_elas = constexpr_varName2Idx("incoh_elas"),
      inelas = constexpr_varName2Idx("inelas"),
      coh_elas = constexpr_varName2Idx("coh_elas"),
      atomdb = constexpr_varName2Idx("atomdb"),
      infofactory = constexpr_varName2Idx("infofactory"),
      scatfactory = constexpr_varName2Idx("scatfactory"),
      absnfactory = constexpr_varName2Idx("absnfactory")
    };

    constexpr VarId constexpr_varIdFromName( const char * );
    constexpr VarGroupId varGroup( VarId ) noexcept;
    constexpr const char * varName( VarId ) noexcept;
    constexpr const VarInfo& varInfo( VarId ) noexcept;
    Optional<VarId> varIdFromName( StrView name );

    enum class CfgVarListMode { TXT_SHORT, TXT_FULL, JSON };
    void dumpCfgVarList( std::ostream&,
                         CfgVarListMode,
                         const char * line_prefix = "" );
  }
}


////////////////////////////
// Inline implementations //
////////////////////////////

namespace NCRYSTAL_NAMESPACE {
  namespace Cfg {

    inline constexpr VarId constexpr_varIdFromName( const char * name )
    {
      return static_cast<VarId>(constexpr_varName2Idx(name));
    }

    inline constexpr VarGroupId varGroup( VarId varid ) noexcept
    {
      return varlist[enumAsInt(varid)].groupId();
    }

    inline constexpr const char * varName( VarId varid ) noexcept
    {
      return varlist[enumAsInt(varid)].name();
    }

    inline constexpr const VarInfo& varInfo( VarId varid ) noexcept {
      return varlist[enumAsInt(varid)];
    }

    inline bool FactNameRequest::excludes( StrView fn ) const noexcept
    {
      //tiny, just do linear search:
      for ( auto& e : m_excluded )
        if ( fn == e )
          return true;
      return false;
    }

    inline FactNameRequest FactNameRequest::withAdditionalExclude( StrView factname ) const
    {
      if ( excludes( factname ) )
        return *this;
      FactNameRequest res( *this );
      res.m_excluded.emplace_back( factname.to_string() );
      return res;
    }

    inline FactNameRequest FactNameRequest::withNoSpecificRequest() const
    {
      FactNameRequest res;
      res.m_excluded = m_excluded;
      return res;
    }
  }
}

#endif

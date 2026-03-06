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

#include "NCrystal/internal/minimc/NCMMC_EngineOpts.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCMMC_ParseCfg.hh"
namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

//Fixme tallyref (and possibly other parameters) is not unit tested.
//fixme: shorter "short descriptions"?

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace {
      namespace TallyFlagsStrDB {
        //NB: Keep in alphabetical order!!  NB: If modifying, also check if we
        //need to modify the description of the "tally" variable below.
        using F = TallyFlags::Flags;
        constexpr static const TallyFlags::value_type vals[]
        = { F::de, F::e, F::highres, F::l, F::lowres, F::mu,
            F::nobreakdown, F::nscat, F::nscat_uw, F::q, F::theta, F::w  };
        constexpr static const char* strs[]
        = { "de", "e","highres", "l", "lowres","mu",
            "nobreakdown", "nscat", "nscat_uw", "q", "theta", "w"  };
        constexpr static const char* strs_descr[]
        = { "energy loss", "energy","more bins", "wavelength", "less bins",
            "cosine scattering angle", "disable breakdown histograms",
            "number of scatterings", "unweighted number of scatterings",
            "scattering vector", "scattering angle", "weight" };
        constexpr static int n = sizeof(vals)/sizeof(*vals);
        static_assert( n == sizeof(strs)/sizeof(*strs), "" );
        static_assert( n == sizeof(vals)/sizeof(*vals), "" );
        static_assert( n == sizeof(strs_descr)/sizeof(*strs_descr), "" );
      }

      //Default values. NOTE: These should be kept synchronized with default
      //values encoded on the structs in NCMMC_EngineOpts.hh! The order of the
      //variables is not important:
      static constexpr auto enginecfg_defvals_str
      = "ignoremiss=0"
        ";nthreads=auto"
        ";absorption=1"
        ";seed=0"
        ";nscatlimit=none"
        ";roulette=0.1,0.01,2"
        ";tallyref=truth"
        ";tally=theta"
        ";tallybins="
        ";tallybreakdown=1"
        ;

      static constexpr auto enginecfg_nscatlimmax = 32000;
      //All the variable names should be listed here (not required to be in any
      //particular order):
      static constexpr auto enginecfg_varnames_str
      = "seed;roulette"
        ";ignoremiss;nthreads;absorption;nscatlimit"
        ";tally;tallybins;tallybreakdown;tallyref";
    }
  }
}

void NCMMC::engineOptsDocsToJSON( std::ostream& os )
{
  //NOTE: Keep synchronised with parseEngineOpts and other functions.

  //Most parameters are documented through static entries in the following
  //array, but a few special keywords like <SHORT> (short_descr), <TALLYLIST>,
  //etc. are expanded dynamically further down. The order of entries in the docs
  //array is not important.

  const char* docs[][4] = {
    {
      "ignoremiss",//name
      "0#1",//example values (split by #)
      "Whether to exclude neutrons missing the geometry from tallies.",//short descr
      "<SHORT> Set to 1 to exclude them and 0 to include them.",//long descr
    },
    {
      "nthreads", "auto#1#2#20",
      "Number of threads to use for simulation.",
      "<SHORT> Either accepts an integer value >=0 or the special"
      " keyword \"auto\"."
    },
    {
      "seed","1000#123456789",
      "Seed value for the random number stream used for the simulation.",
      "<SHORT> Accepts any 64 non-negative integral value that fits into a"
      " 64bit unsigned integer. Note that even with the same seed, results are"
      " not expected to be reproducible unless also using nthreads=1."
    },
    {
      "tallyref", "truth#src",
      "Definition of incoming quantities in tallies.",
      "<SHORT> This affects the calculation of tally quantities"
      " like \"q\", \"theta\", or \"mu\" which needs parameters from"
      " the incoming neutron state.  The default value of \"truth\""
      " implies that the actual original quantity of each individual"
      " neutron as it came from the source is used. A value of"
      " \"src\" implies that the reference value will be taken as"
      " the mean value provided by the source. In the special case"
      " of energy specified as a wavelength (via the \"wl\" source"
      " parameter), the mean wavelength rather than energy will be"
      " used where that makes a difference."
    },
    {
      "absorption", "0#1",
      "Whether to include absorption effects.",
      "<SHORT> If set to 1 (the default) they will be included implicitly"
      " through attenuation of neutrons weights. If set to 0 it will be as"
      " if the absorption cross section was always 0 barn."
    },
    {
      "nscatlimit", "none#0#1#2#10#999",
      "Maximum number of scatterings to allow for each neutron.",
      "<SHORT> If set to none (the default), no limit will by explicitly"
      " applied. If set to any other non-negative integral value, the"
      " scattering cross section of a neutron will essentially become 0"
      " just after a neutron has undergone that many scatterings. A typical"
      " usage might be to set this value to 1, in order to exclude multiple"
      " scattering effects."
    },
    {
      "roulette", "0.5,1e-4,5",
      "Roulette parameters for the neutron termination strategy.",
      "<SHORT> At each simulation step inside the geometry, neutrons are"
      " split into a transmitted (and tallied) part as well as a part"
      " undergoing a forced collision. Weights are handled appropriately,"
      " decreasing in magnitude after each step. To maximise utility of"
      " simulation time, a roulette choice is made for neutrons having"
      " undergone at least N scatterings and whose weights are below a"
      " threshold value W. The choice allow the neutron to survive with"
      " a probability of P, and if it survives its weight will be"
      " increased by a factor of 1/P. This allows the simulation to avoid"
      " using most time on insignificant neutrons, while still providing the"
      " correct tally distributions in the limit of large statistics. The"
      " roulette parameter thus needs three comma-separated values,"
      " encoded as \"P,W,N\"."
    },
    {
      "tally", "q#mu,q,e,l#q",
      "Quantities to tally (available: <TALLYLISTSHORT>).",
      "Comma separated list of quantities to tally."
      " The quantities available for tally are: <TALLYLISTLONG>."
    },
    {
      "tallybins",
      "#q:1000:0.0:20.0#q:1000:0.0:20.0,l:200:3.0:5.0#+#-#q:1000:0.0,+",
      "Modify binning of tally histograms.",
      "<SHORT> This is supplied as a comma separated list of"
      " \"<tallyname>,<nbins>,<lower edge>,<upper edge>\". If a tally"
      " is enabled via the tally keyword, it will use the binning"
      " specified in this list, if any. If a binning is not specified,"
      " an automatic choice is made for the binning. Special entries"
      " \"+\" or \"-\" can be used to generally increase or"
      " decrease binning in such cases."
    },
    {
      "tallybreakdown", "0#1",
      "Whether to also tally per-component.",
      "If set to 1 (the default), tallies will also include per-component"
      " breakdowns, allowing to distinguish contributions depending on"
      " number of scatterings (0, 1, 2+) and whether the"
      " scatterings were elastic or not."
    },
  };

  std::string short_descr_tallylist;
  std::string long_descr_tallylist;
  {
    std::ostringstream sss, ssl;
    bool first = true;
    for ( auto i : ncrange(TallyFlagsStrDB::n) ) {
      if ( !( TallyFlags::Flags::ALLHISTS & TallyFlagsStrDB::vals[i] ) )
        continue;
      if ( first ) {
        first = false;
      } else if ( i ) {
        sss << ',';
        ssl << ", ";
      }
      sss << TallyFlagsStrDB::strs[i];
      ssl << TallyFlagsStrDB::strs[i] << " ("
          <<TallyFlagsStrDB::strs_descr[i]<<')';
    }
    short_descr_tallylist = std::move(sss).str();
    long_descr_tallylist = std::move(ssl).str();
  }

  constexpr auto nvars = sizeof(docs)/sizeof(*docs);
  auto defvals = parseMMCCfg::tokenize( enginecfg_defvals_str ).tokens;
  auto ref_varnames = StrView(enginecfg_varnames_str).splitTrimmedNoEmpty(';');
  nc_assert_always( ref_varnames.size() == nvars );
  os << '{';
  for ( auto i : ncrange(nvars) ) {
    auto docval = [i,&docs](int j) { return StrView(docs[i][j]).trimmed(); };
    auto v_name = docval(0);
    auto v_examplevals = docval(1);
    auto v_descr_short = docval(2);
    nc_assert_always(v_descr_short.size()<70);
    //Expand special words in the descriptions:
    auto v_descr_long = docval(3);
    std::string tmpbuf;
    if ( v_descr_long.contains("<SHORT>") ) {
      tmpbuf = v_descr_long.to_string();
      strreplace(tmpbuf,"<SHORT>",v_descr_short.to_string());
      v_descr_long = tmpbuf;
    }
    std::string tmpbuf2, tmpbuf3;
    if ( v_name == "tally" ) {
      nc_assert( v_descr_short.contains("<TALLYLISTSHORT>") );
      nc_assert( v_descr_long.contains("<TALLYLISTLONG>") );
      tmpbuf2 = v_descr_short.to_string();
      strreplace(tmpbuf2,"<TALLYLISTSHORT>",short_descr_tallylist);
      v_descr_short = tmpbuf2;
      tmpbuf3 = v_descr_long.to_string();
      strreplace(tmpbuf3,"<TALLYLISTLONG>",long_descr_tallylist);
      v_descr_long = tmpbuf3;
    }
    //Find default value if any:
    Optional<StrView> defval;
    for ( auto& e : defvals ) {
      if ( e.first == v_name ) {
        defval = e.second;
        break;
      }
    }
#ifndef NDEBUG
    {
      int dbg_in_ref = 0;
      for (auto& e : ref_varnames )
        if ( e == v_name )
          ++dbg_in_ref;
      nc_assert(dbg_in_ref==1);
    }
#endif
    if ( v_descr_long.empty() )
      v_descr_long = v_descr_short;

    //Example values:
    auto exvals = v_examplevals.splitTrimmed('#');

    //Output json:
    if ( i )
      os << ',';
    streamJSON(os,v_name);
    os << ':';
    streamJSONDictEntry( os, "name", v_name, JSONDictPos::FIRST );
    streamJSONDictEntry( os, "default_value", defval );
    streamJSONDictEntry( os, "example_values", exvals );
    streamJSONDictEntry( os, "descr_short", v_descr_short );
    streamJSONDictEntry( os, "descr_long", v_descr_long, JSONDictPos::LAST );
  }
  os << '}';

}

NCMMC::EngineOpts NCMMC::parseEngineOpts( StrView raw_eoptsstr )
{
  //NOTE: Keep synchronised with engineOptsDocsToJSON and other functions.
  namespace PMC = parseMMCCfg;

  auto tokeninfo = PMC::tokenize( raw_eoptsstr );
  auto& tokens = tokeninfo.tokens;

  StrView engine_name =  tokeninfo.mainName;
  if ( !engine_name.has_value() )
    engine_name = "std";//Default to "std" engine

  //For now we only have the one engine (fixme: add an "analogue" engine?)
  if (engine_name != "std" )
    NCRYSTAL_THROW2(BadInput,"Invalid MiniMC engine \""<<engine_name<<"\"");

  NCMMC::EngineOpts res;
  //Defaults (with asserts that they are compatible with the defaults on the
  //struct): The default values must be kept synchronised with the stream
  //operator (which omits values at default values) and the default values on
  //the EngineOpt struct itself.
  const char * defaults = enginecfg_defvals_str;
  static_assert( EngineOpts::IgnoreMiss::Default ==
                 EngineOpts::IgnoreMiss::NO, "" );
  static_assert( EngineOpts::IncludeAbsorption::Default ==
                 EngineOpts::IncludeAbsorption::YES, "" );
  using TR = EngineOpts::TallyReference;
  static_assert( TR::Default == TR::Truth, "");
  nc_assert( res.nthreads.indicatesAutoDetect() );

  //Apply defaults and error in case of unknown parameters:
  PMC::applyDefaults( tokens, defaults );
  PMC::checkNoUnknown(tokens,enginecfg_varnames_str,"engine");
  {
    auto sv_tr = PMC::getValue_str( tokens, "tallyref" );
    if ( sv_tr == "truth" ) {
      res.tallyRef = TR::Truth;
    } else if ( sv_tr == "src" ) {
      res.tallyRef = TR::Source;
    } else {
      NCRYSTAL_THROW(BadInput,"tallyref must be \"truth\" or \"src\"");
    }
  }
  {
    auto sv_roulette = PMC::getValue_str( tokens, "roulette" );
    auto parts = sv_roulette.splitTrimmed<3>(',');
    Optional<double> val_psurvive;
    Optional<double> val_thr_w;
    Optional<std::int32_t> val_thr_nscat;
    if ( parts.size() == 3 ) {
      val_psurvive = parts.at(0).toDbl();
      val_thr_w = parts.at(1).toDbl();
      val_thr_nscat = parts.at(2).toInt32();
    }
    if ( !(val_thr_w.value_or(0.5) > 0.0 ) )
      val_thr_w = NullOpt;
    if ( ncisnanorinf(val_thr_w.value_or(0.5)) )
      val_thr_w = NullOpt;
    if ( !(val_psurvive.value_or(0.5) > 0.0 ) )
      val_psurvive = NullOpt;
    if ( !(val_psurvive.value_or(0.5) < 1.0 ) )
      val_psurvive = NullOpt;
    if ( ncisnanorinf(val_psurvive.value_or(0.5)) )
      val_psurvive = NullOpt;
    if ( val_thr_nscat.value_or(2) > enginecfg_nscatlimmax )
      val_thr_nscat = NullOpt;
    if ( val_thr_nscat.value_or(2) < 0 )
      val_thr_nscat = 0;
    if ( !val_psurvive.has_value()
         || !val_thr_w.has_value()
         || !val_thr_nscat.has_value() )
      NCRYSTAL_THROW2(BadInput,
                      "Invalid roulette options: \""<<sv_roulette
                      <<"\" Please provide valid values in the form"
                      " \"psurvival,weight_threshold,nscat_threshold\"");
    res.roulette.survival_probability = val_psurvive.value();
    res.roulette.weight_threshold = val_thr_w.value();
    static_assert( std::numeric_limits<int>::max()
                   >= std::numeric_limits<std::int32_t>::max(), "" );
    res.roulette.nscat_threshold = static_cast<int>(val_thr_nscat.value());
  }

  {
    auto sv_seed = PMC::getValue_str( tokens, "seed" );
    auto seedval = sv_seed.toUInt64();
    if ( !seedval.has_value() )
      NCRYSTAL_THROW2(BadInput,"Invalid seed \""<<sv_seed<<"\".");
    res.seed = seedval.value();
  }

  res.ignoreMiss = ( PMC::getValue_bool(tokens,"ignoremiss")
                       ? EngineOpts::IgnoreMiss::YES
                       : EngineOpts::IgnoreMiss::NO );
  res.includeAbsorption = ( PMC::getValue_bool(tokens,"absorption")
                            ? EngineOpts::IncludeAbsorption::YES
                            : EngineOpts::IncludeAbsorption::NO );

  //nthreads must be "auto" or an unsigned integer
  auto nthreadsval_str = PMC::getValue_str(tokens,"nthreads");
  if ( nthreadsval_str == "auto" ) {
    //do nothing, res.nthreads already defaults to "auto".
  } else {
    auto v = PMC::getValue_sizet( tokens, "nthreads" );
    //Values >=9999 also indicates auto
    if ( v < 9999 )
      res.nthreads = ThreadCount{ static_cast<std::uint32_t>(v) };
  }

  if ( PMC::getValue_str( tokens, "nscatlimit" ) != "none" ) {
    auto val = PMC::getValue_sizet(tokens,"nscatlimit");
    if ( val > enginecfg_nscatlimmax )
      NCRYSTAL_THROW2(BadInput,
                      "nscatlimit can be at most "<<enginecfg_nscatlimmax);
    res.nScatLimit = val;
  }

  using F = TallyFlags::Flags;
  if ( PMC::hasValue(tokens,"tally") ) {
    auto tallyflags_strlist =
      PMC::getValue_str(tokens,"tally")
      .splitTrimmedNoEmpty<TallyFlags::strlist_type::nsmall>(',');
    auto tf = TallyFlags(tallyflags_strlist);
    constexpr auto tf_nonhist = F::ALL & (~F::ALLHISTS);
    auto tf_forbidden = tf.getValue() & tf_nonhist;
    static_assert( tf_nonhist==( F::nobreakdown | F::lowres | F::highres), "" );
    if ( tf_forbidden & F::nobreakdown )
      NCRYSTAL_THROW2(BadInput,
                      "Invalid tally \""
                      <<TallyFlags::singleFlagToString(F::nobreakdown)
                      <<"\" (use \"tallybreakdown=0\" to disable breakdowns)");
    if ( tf_forbidden & F::lowres )
      NCRYSTAL_THROW2(BadInput,
                      "Invalid tally \""
                      <<TallyFlags::singleFlagToString(F::lowres)
                      <<"\" (use \"tallybins=-\" to decrease binnings)");
    if ( tf_forbidden & F::highres )
      NCRYSTAL_THROW2(BadInput,
                      "Invalid tally \""
                      <<TallyFlags::singleFlagToString(F::highres)
                      <<"\" (use \"tallybins=+\" to increase binnings)");
    nc_assert_always( !tf_forbidden );
    res.tallyFlags = tf;
  }

  //default is to include breakdown hists:
  static_assert( !(F::DEFAULT & F::nobreakdown),"" );
  if ( !PMC::getValue_bool(tokens,"tallybreakdown") )
    res.tallyFlags.add(F::nobreakdown);

  {
    auto strlist = PMC::getValue_str_allowempty(tokens,"tallybins")
      .splitTrimmedNoEmpty<TallyFlags::strlist_type::nsmall>(',');
    for ( auto& bstr : strlist ) {
      //Special support for "+" and "-" to enable highres/lowres TallyFlags:
      if ( bstr == '+' || bstr == '-' ) {
        res.tallyFlags.add(bstr=='+'?F::highres:F::lowres);
        continue;
      }
      TallyFlags::value_type v = 0;
      Optional<std::int32_t> nbins;
      Optional<double> xmin;
      Optional<double> xmax;
      auto parts = bstr.splitTrimmedNoEmpty<4>(':');
      if ( parts.size() == 4 ) {
        v = TallyFlags::lookup(parts.at(0));
        nbins = parts.at(1).toInt32();
        xmin = parts.at(2).toDbl();
        xmax = parts.at(3).toDbl();
      }
      if ( !v || !TallyFlags::isSingleFlag(v) || !( v & F::ALLHISTS ) )
        NCRYSTAL_THROW2(BadInput,"Invalid enginecfg tallybins tally name \""
                        <<parts.at(0) <<"\".");
      if ( !nbins.has_value() || !xmin.has_value() || !xmax.has_value() )
        NCRYSTAL_THROW2(BadInput,"Invalid enginecfg tallybins entry \""
                        <<bstr <<"\" (should have the form"
                        " \"tallyname:nbins:xmin:xmax\").");
      Hists::Binning b( nbins.value(), xmin.value(), xmax.value() );
      b.validate();
      res.tallyBinnings.add( v, b );
    }
  }

  return res;
}

std::ostream& NCMMC::operator<<( std::ostream& os, const EngineOpts& eopts )
{
  os << "EngineOpts("<<engineOptsToString(eopts)<<')';
  return os;
}

std::string NCMMC::engineOptsToString( const EngineOpts& eopts )
{
  using EO = EngineOpts;
  std::ostringstream ss;
  bool is_empty(true);
  auto delim = [&ss,&is_empty] () -> std::ostringstream&
  {
    if (!is_empty)
      ss<<';';
    is_empty=false;
    return ss;
  };

  if ( eopts.ignoreMiss != EO::IgnoreMiss::Default )
    delim()<<"ignoremiss="
           << ( eopts.ignoreMiss == EO::IgnoreMiss::YES ? '1' : '0' );

  if ( eopts.includeAbsorption != EO::IncludeAbsorption::Default )
    delim()<<"absorption="
           << ( eopts.includeAbsorption == EO::IncludeAbsorption::YES ? '1' : '0' );

  using TR = EngineOpts::TallyReference;
  static_assert( TR::Default == TR::Truth, "");
  if ( eopts.tallyRef != TR::Truth ) {
    nc_assert( eopts.tallyRef == TR::Source );
    delim()<<"tallyref=src";
  }

  if ( eopts.seed != 0 )
    delim()<<"seed="<<eopts.seed;

  if ( ! ( eopts.roulette == RouletteOptions{} ) )
    delim()<<"roulette="<<fmt(eopts.roulette.survival_probability)
           <<','<<fmt(eopts.roulette.weight_threshold)
           <<','<<fmt(eopts.roulette.nscat_threshold);

  if ( eopts.nScatLimit.has_value() )
    delim()<<"nscatlimit="<<eopts.nScatLimit.value();

  if ( !eopts.nthreads.indicatesAutoDetect() )
    delim()<<"nthreads="<<eopts.nthreads.get();

  using F = TallyFlags::Flags;
  const auto tflags_hists = eopts.tallyFlags.getValue() & F::ALLHISTS;

  if ( tflags_hists != ( F::DEFAULT & F::ALLHISTS ) ) {
    delim()<<"tally=";
    bool firsttl = true;
    for ( auto& e : TallyFlags(tflags_hists).toStringList() ) {
      if (!firsttl)
        ss<<',';
      firsttl=false;
      ss<<e;
    }
  }

  static_assert( !(F::DEFAULT & F::nobreakdown),"" );
  if ( eopts.tallyFlags.getValue() & F::nobreakdown )
    delim()<<"tallybreakdown=0";

  {
    auto& tbins = eopts.tallyBinnings;
    bool tbins_lowres = eopts.tallyFlags.getValue() & F::lowres;
    bool tbins_highres = eopts.tallyFlags.getValue() & F::highres;
    if ( tbins.flagsAffected().getValue() || tbins_lowres || tbins_highres ) {
      delim() << "tallybins=";
      tbins.toString(ss);
      if ( tbins_highres ) {
        if ( tbins.flagsAffected().getValue() )
          ss << ',';
        ss << '+';
      }
      if ( tbins_lowres ) {
        if ( tbins.flagsAffected().getValue() || tbins_highres )
          ss << ',';
        ss << '-';
      }
    }
  }

  return ss.str();
}

void NCMMC::engineOptsToJSON(std::ostream& os, const EngineOpts& eopts)
{
  using EO = EngineOpts;
  //We could include engine name if we have more at some point:
  streamJSONDictEntry( os, "name", "std", JSONDictPos::FIRST );
  streamJSONDictEntry( os, "cfgstr", engineOptsToString(eopts)  );
  os << ",\"decoded\":";
  streamJSONDictEntry( os, "seed", eopts.seed, JSONDictPos::FIRST );
  {
    os << ",\"roulette\":";
    streamJSONDictEntry( os, "psurvival", eopts.roulette.survival_probability,
                         JSONDictPos::FIRST );
    streamJSONDictEntry( os, "weight_threshold",
                         eopts.roulette.weight_threshold );
    streamJSONDictEntry( os, "nscat_threshold", eopts.roulette.nscat_threshold,
                         JSONDictPos::LAST );
  }
  streamJSONDictEntry( os, "ignoremiss",
                       bool( eopts.ignoreMiss == EO::IgnoreMiss::YES ) );
  streamJSONDictEntry( os, "absorption",
                       bool( eopts.includeAbsorption
                             == EO::IncludeAbsorption::YES ) );

  using TR = EngineOpts::TallyReference;
  nc_assert( eopts.tallyRef == TR::Source || eopts.tallyRef == TR::Truth );
  streamJSONDictEntry( os, "tallyref",
                       ( eopts.tallyRef == TR::Truth ? "truth" : "src" ) );

  if ( eopts.nScatLimit.has_value() ) {
    streamJSONDictEntry( os, "nscatlimit", eopts.nScatLimit.value() );
  } else {
    streamJSONDictEntry( os, "nscatlimit", json_null_t{} );
  }

  using F = TallyFlags::Flags;
  const auto tflags_hists = eopts.tallyFlags.getValue() & F::ALLHISTS;
  streamJSONDictEntry( os, "tally",
                       TallyFlags(tflags_hists).toStringList() );
  streamJSONDictEntry( os, "tallybreakdown",
                       !(eopts.tallyFlags.getValue() & F::nobreakdown) );
  {
    os << ",\"tallybins\":{\"overrides\":";
    eopts.tallyBinnings.toJSON(os);
    os << ",\"autobin_flags\":[";
    const auto tf = eopts.tallyFlags.getValue();
    if ( tf & F::highres )
      os << "\"+\"";
    if ( tf & F::lowres )
      os << (tf&F::highres?",":"") << "\"-\"";
    os << "]}";
  }

  if ( eopts.nthreads.indicatesAutoDetect() ) {
    streamJSONDictEntry( os, "nthreads", "auto", JSONDictPos::LAST );
  } else {
    streamJSONDictEntry( os, "nthreads", eopts.nthreads.get(),
                         JSONDictPos::LAST );
  }
  os << '}';
}

bool NCMMC::TallyFlags::isSingleFlag( value_type v )
{
  for ( int i = 0; i < TallyFlagsStrDB::n; ++i ) {
    if ( v == TallyFlagsStrDB::vals[i] )
      return true;
  }
  return false;
}

const char * NCMMC::TallyFlags::singleFlagToString( value_type v )
{
  for ( int i = 0; i < TallyFlagsStrDB::n; ++i ) {
    if ( v == TallyFlagsStrDB::vals[i] )
      return TallyFlagsStrDB::strs[i];
  }
  NCRYSTAL_THROW2( BadInput, "Not a single flag (must have"
                   " exactly one bit set): \""<<v<<'"');
  return "";
}

NCMMC::TallyFlags::strlist_type
NCMMC::TallyFlags::toStringList() const
{
  strlist_type res;
  for ( int i = 0; i < TallyFlagsStrDB::n; ++i ) {
    if ( has(TallyFlagsStrDB::vals[i]) )
      res.emplace_back(TallyFlagsStrDB::strs[i]);
  }
  return res;
}

NCMMC::TallyFlags::value_type
NCMMC::TallyFlags::lookup( StrView sv )
{
  sv = sv.trimmed();

  //collections:
  constexpr StrView sv_default = StrView::make("default");
  constexpr StrView sv_all = StrView::make("all");
  constexpr StrView sv_allhists = StrView::make("allhists");

  if ( sv == sv_default )
    return Flags::DEFAULT;
  if ( sv == sv_all )
    return Flags::ALL;
  if ( sv == sv_allhists )
    return Flags::ALLHISTS;

  //flags:
  if ( !sv.empty() ) {
    char c = sv[0];
    for ( int i = 0; i < TallyFlagsStrDB::n; ++i ) {
      if ( c == TallyFlagsStrDB::strs[i][0] && sv == TallyFlagsStrDB::strs[i] )
        return TallyFlagsStrDB::vals[i];
    }
  }

  NCRYSTAL_THROW2( BadInput, "Not a valid tally name: \""<<sv<<'"');
  return Flags::NONE;
}

NCMMC::TallyFlags::TallyFlags( const strlist_type& sl )
  : m_value(Flags::NONE)
{
  static_assert( std::is_same<decltype(m_value),value_type>::value, "");
  for ( auto& sv : sl ) {
    StrView svt = sv.trimmed();
    if ( !svt.empty() ) {
      add(svt);
    }
  }
}

void NCMMC::TallyBinningOverrides::toString( std::ostream& os ) const
{
  if (!m_all)
    return;
  bool first = true;
  //Output in sorted (by names) order:
  SmallVector<std::pair<StrView,std::size_t>,decltype(m_db)::nsmall> nb;
  for ( auto i : ncrange(m_db.size()) )
    nb.emplace_back( TallyFlags::singleFlagToString( m_db[i].first ), i );
  std::sort(nb.begin(),nb.end());
  for ( auto& e : nb ) {
    if ( first )
      first = false;
    else
      os << ',';
    const auto& binning = m_db[e.second].second;
    os << e.first
       << ':' << binning.nbins
       << ':' << fmt(binning.xmin)
       << ':' << fmt(binning.xmax);
  }
}

void NCMMC::TallyBinningOverrides::toJSON( std::ostream& os ) const
{
  using F = TallyFlags::Flags;
  os << '{';
  bool first = true;
  const Binning dummy( 100, 0.0, 1.0 );
  for ( auto i : ncrange(TallyFlagsStrDB::n) ) {
    const auto f = TallyFlagsStrDB::vals[i];
    if ( !( F::ALLHISTS & f ) )
      continue;
    if ( first )
      first = false;
    else
      os << ',';
    streamJSON(os,TallyFlagsStrDB::strs[i]);
    os<<':';
    if ( !( m_all & f ) ) {
      streamJSON(os,json_null_t{});
    } else {
      auto& bins = lookup( f, dummy );
      os<<'[';
      streamJSON(os,bins.nbins);
      os<<',';
      streamJSON(os,bins.xmin);
      os<<',';
      streamJSON(os,bins.xmax);
      os<<']';
    }
  }
  os << '}';
}

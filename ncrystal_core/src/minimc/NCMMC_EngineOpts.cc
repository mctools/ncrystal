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

#include "NCrystal/internal/minimc/NCMMC_EngineOpts.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCMMC_ParseCfg.hh"
namespace NC = NCrystal;
namespace NCMMC = NCrystal::MiniMC;

NCMMC::EngineOpts NCMMC::parseEngineOpts( StrView raw_srcstr )
{
  namespace PMC = parseMMCCfg;
  auto tokens = PMC::tokenize( raw_srcstr );
  StrView src_name =  PMC::mainName( tokens );
  if ( !src_name.has_value() )
    src_name = "std";//Default to "std" engine

  //For now we only have the one engine (fixme: add an "analogue" engine?)
  if (src_name != "std" )
    NCRYSTAL_THROW2(BadInput,"Invalid MiniMC engine \""<<src_name<<"\"");

  NCMMC::EngineOpts res;
  //Defaults (with asserts that they are compatible with the defaults on the
  //struct): The default values must be kept synchronised with the stream
  //operator (which omits values at default values) and the default values on
  //the EngineOpt struct itself.
  const char * defaults
    = "ignoremiss=0;nthreads=auto;absorption=1;seed=0";
  static_assert( EngineOpts::IgnoreMiss::Default ==
                 EngineOpts::IgnoreMiss::NO, "" );
  static_assert( EngineOpts::IncludeAbsorption::Default ==
                 EngineOpts::IncludeAbsorption::YES, "" );
  nc_assert( res.nthreads.indicatesAutoDetect() );

  //Apply defaults and error in case of unknown parameters:
  PMC::applyDefaults( tokens, defaults );
  PMC::checkNoUnknown(tokens,
                      //
                      "seed;roulette"
                      ";ignoremiss;nthreads;absorption;nscatlimit"
                      ";beamdirx;beamdiry;beamdirz;beamenergy"
                      ";tally;tallybins",
                      //
                      "engine");

  int nbeamdir = ( (PMC::hasValue(tokens,"beamdirx")?1:0)
                   + (PMC::hasValue(tokens,"beamdiry")?1:0)
                   + (PMC::hasValue(tokens,"beamdirz")?1:0) );
  if ( nbeamdir ) {
    if ( nbeamdir<3)
      NCRYSTAL_THROW(BadInput,"Must set all or none of the parameters:"
                     " \"beamdirx\", \"beamdiry\", and \"beamdirz\".");
    Vector bd( PMC::getValue_dbl( tokens, "beamdirx" ),
               PMC::getValue_dbl( tokens, "beamdiry" ),
               PMC::getValue_dbl( tokens, "beamdirz" ) );
    if ( bd.mag2()<1e-12 )
      NCRYSTAL_THROW(BadInput,"Provided beamdir vector is too"
                     " close to a null vector.");
    if ( bd.mag2()>1e99 )
      NCRYSTAL_THROW(BadInput,"Provided beamdir vector has too "
                     "large components.");
    //Normalise and set:
    res.tallyBeamDir = bd.unit().as<NeutronDirection>();
  }

  if ( PMC::hasValue(tokens,"beamenergy") ) {
    double be = PMC::getValue_dbl( tokens, "beamenergy" );
    if ( !(be>0.0) || !std::isfinite(be) )
      NCRYSTAL_THROW(BadInput,
                     "Provided beamenergy must be positive and finite");
    res.tallyBeamEnergy = NeutronEnergy{ be };
  }


  if ( PMC::hasValue( tokens, "roulette" ) ) {
      auto sv_roulette = PMC::getValue_str( tokens, "roulette" );
      auto parts = sv_roulette.splitTrimmedNoEmpty<3>(',');
      Optional<double> val_psurvive;
      Optional<double> val_thr_w;
      Optional<int> val_thr_nscat;
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
      if ( val_thr_nscat.value_or(2) > 30000 )
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

  if ( PMC::hasValue(tokens,"nscatlimit")
       && PMC::getValue_str( tokens, "nscatlimit" ) != "none" ) {
    auto val = PMC::getValue_sizet(tokens,"nscatlimit");
    constexpr auto nscatlimmax = 32000;
    if ( val > nscatlimmax )
      NCRYSTAL_THROW2(BadInput,
                      "nscatlimit can be at most "<<nscatlimmax);
    res.nScatLimit = val;
  }

  if ( PMC::hasValue(tokens,"tally") ) {
    auto tallyflags_strlist =
      PMC::getValue_str(tokens,"tally")
      .splitTrimmedNoEmpty<TallyFlags::strlist_type::nsmall>(',');
    res.tallyFlags = TallyFlags(tallyflags_strlist);
  }

  if ( PMC::hasValue(tokens,"tallybins") ) {
    auto strlist = PMC::getValue_str(tokens,"tallybins")
      .splitTrimmedNoEmpty<TallyFlags::strlist_type::nsmall>(',');
    for ( auto& bstr : strlist ) {
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
      if ( !v || !TallyFlags::isSingleFlag(v)
           || !(v&TallyFlags::Flags::ALLHISTS) )
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
  if ( eopts.ignoreMiss != EO::IgnoreMiss::Default ) {
    is_empty = false;
    ss<<"ignoremiss="
      << ( eopts.ignoreMiss == EO::IgnoreMiss::YES ? '1' : '0' );
  }
  if ( eopts.includeAbsorption != EO::IncludeAbsorption::Default ) {
    if (!is_empty)
      ss << ';';
    is_empty = false;
    ss<<"absorption="
      << ( eopts.includeAbsorption == EO::IncludeAbsorption::YES ? '1' : '0' );
  }
  if ( !eopts.nthreads.indicatesAutoDetect() ) {
    if (!is_empty)
      ss << ';';
    is_empty = false;
    ss<<"nthreads="<<eopts.nthreads.get();
  }
  if ( eopts.tallyBeamDir.has_value() ) {
    ss<<";beamdirx="<<fmt(eopts.tallyBeamDir.value()[0])
      <<";beamdiry="<<fmt(eopts.tallyBeamDir.value()[1])
      <<";beamdirz="<<fmt(eopts.tallyBeamDir.value()[2]);
  }
  if ( eopts.seed != 0 )
    ss<<";seed="<<eopts.seed;
  if ( ! ( eopts.roulette == RouletteOptions{} ) ) {
    ss<<";roulette="<<fmt(eopts.roulette.survival_probability)
      <<','<<fmt(eopts.roulette.weight_threshold)
      <<','<<fmt(eopts.roulette.nscat_threshold);
  }
  if ( eopts.tallyBeamEnergy.has_value() )
    ss<<";beamenergy="<<fmt(eopts.tallyBeamEnergy.value().dbl());;
  if ( eopts.nScatLimit.has_value() )
    ss<<";nscatlimit="<<eopts.nScatLimit.value();
  if ( !eopts.nthreads.indicatesAutoDetect() ) {
    if (!is_empty)
      ss << ';';
    is_empty = false;
    ss<<"nthreads="<<eopts.nthreads.get();
  }
  if ( eopts.tallyFlags.getValue() != TallyFlags().getValue() ) {
    if (!is_empty)
      ss << ';';
    is_empty = false;
    ss<<"tally=";
    bool firsttl = true;
    for ( auto& e : eopts.tallyFlags.toStringList() ) {
      if (!firsttl)
        ss<<',';
      firsttl=false;
      ss<<e;
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
  if ( eopts.tallyBeamDir.has_value() ) {
    VectD vtbd = { eopts.tallyBeamDir.value()[0],
                   eopts.tallyBeamDir.value()[1],
                   eopts.tallyBeamDir.value()[2] };
    //fixme: also make srcBeamDir + energy avail in json
    streamJSONDictEntry( os, "beamdir", vtbd );
  } else {
    streamJSONDictEntry( os, "beamdir", json_null_t{} );
  }
  if ( eopts.tallyBeamEnergy.has_value() ) {
    streamJSONDictEntry( os, "beamenergy",
                         eopts.tallyBeamEnergy.value().dbl() );
  } else {
    streamJSONDictEntry( os, "beamenergy", json_null_t{} );
  }
  if ( eopts.nScatLimit.has_value() ) {
    streamJSONDictEntry( os, "nscatlimit", eopts.nScatLimit.value() );
  } else {
    streamJSONDictEntry( os, "nscatlimit", json_null_t{} );
  }
  streamJSONDictEntry( os, "tally",
                       eopts.tallyFlags.toStringList() );
  if ( eopts.nthreads.indicatesAutoDetect() ) {
    streamJSONDictEntry( os, "nthreads", "auto", JSONDictPos::LAST );
  } else {
    streamJSONDictEntry( os, "nthreads", eopts.nthreads.get(),
                         JSONDictPos::LAST );
  }
  os << '}';
}

namespace NCRYSTAL_NAMESPACE {
  namespace MiniMC {
    namespace {
      struct TallyFlagsStrDB {
        //NB: Keep in alphabetical order!!
        using F = TallyFlags::Flags;
        constexpr static const TallyFlags::value_type vals[]
        = { F::cosmu, F::e, F::highres, F::l, F::lowres, F::mu,
            F::nobreakdown, F::nscat, F::nscat_uw, F::q, F::w  };
        constexpr static const char* strs[]
        = { "cosmu", "e","highres", "l", "lowres", "mu",
            "nobreakdown", "nscat", "nscat_uw", "q", "w"  };
        constexpr static int n = sizeof(vals)/sizeof(*vals);
        static_assert( n == sizeof(strs)/sizeof(*strs), "" );
        static_assert( n == sizeof(vals)/sizeof(*vals), "" );
      };
    }
  }
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

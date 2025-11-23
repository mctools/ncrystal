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
  if (src_name != "std" )
    NCRYSTAL_THROW2(BadInput,"Invalid MiniMC engine \""<<src_name<<"\"");

  NCMMC::EngineOpts res;

  //Defaults (with asserts that they are compatible with the defaults on the
  //struct): The default values must be kept synchronised with the stream
  //operator (which omits values at default values) and the default values on
  //the EngineOpt struct itself.
  const char * defaults = "ignoremiss=0;tallybreakdown=1;nthreads=auto";
  static_assert( EngineOpts::IgnoreMiss::Default ==
                 EngineOpts::IgnoreMiss::NO, "" );
  static_assert( EngineOpts::TallyBreakdown::Default ==
                 EngineOpts::TallyBreakdown::YES, "" );
  nc_assert( res.nthreads.indicatesAutoDetect() );

  //Apply defaults and error in case of unknown parameters:
  PMC::applyDefaults( tokens, defaults );
  PMC::checkNoUnknown(tokens,"ignoremiss;tallybreakdown;nthreads","engine");

  res.ignoreMiss = ( PMC::getValue_bool(tokens,"ignoremiss")
                       ? EngineOpts::IgnoreMiss::YES
                       : EngineOpts::IgnoreMiss::NO );
  res.tallyBreakdown = ( PMC::getValue_bool(tokens,"tallybreakdown")
                         ? EngineOpts::TallyBreakdown::YES
                         : EngineOpts::TallyBreakdown::NO );
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
  if ( eopts.tallyBreakdown != EO::TallyBreakdown::Default ) {
    if (!is_empty)
      ss << ';';
    is_empty = false;
    ss<<"tallybreakdown="
      << ( eopts.tallyBreakdown == EO::TallyBreakdown::YES ? '1' : '0' );
  }
  if ( !eopts.nthreads.indicatesAutoDetect() ) {
    if (!is_empty)
      ss << ';';
    is_empty = false;
    ss<<"nthreads="<<eopts.nthreads.get();
  }
  return ss.str();
}

void NCMMC::engineOptsToJSON(std::ostream& os, const EngineOpts& eopts)
{
  using EO = EngineOpts;
  //We could include engine name if we have more at some point (but note that if
  //makes it a tiny bit more difficult to construct a eopts-string from a json
  //dict)
  //streamJSONDictEntry( os, "engine", "std", JSONDictPos::FIRST );
  streamJSONDictEntry( os, "ignoremiss",
                       bool( eopts.ignoreMiss == EO::IgnoreMiss::YES ),
                       JSONDictPos::FIRST );
  streamJSONDictEntry( os, "tallybreakdown",
                       bool( eopts.tallyBreakdown == EO::TallyBreakdown::YES ),
                       JSONDictPos::OTHER );
  if ( eopts.nthreads.indicatesAutoDetect() ) {
    streamJSONDictEntry( os, "nthreads", "auto", JSONDictPos::LAST );
  } else {
    streamJSONDictEntry( os, "nthreads", eopts.nthreads.get(),
                         JSONDictPos::LAST );
  }
}

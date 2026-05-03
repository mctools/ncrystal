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

#include "NCrystal/internal/query/NCQuery.hh"
#include "NCrystal/internal/minimc/NCMMC_Query.hh"
#include "NCrystal/internal/cfgutils/NCCfgTypes.hh"
#include "NCSABQuery.hh"
#include "NCUtilQuery.hh"

namespace NC = NCrystal;

void NC::JSONQuery( std::ostream& os, const Query& query )
{
  if ( query.empty() )
    NCRYSTAL_THROW2(BadInput,
                    "Invalid JSON query (query is empty)");

  for ( auto& e : query ) {
    //Enforce lack of leading '-' in query strings to make it easier to accept
    //query strings from the command line, and in general not allowing more
    //chars than absolutely necessary to accept NCrystal cfg-strings.
    auto badchar
      = findForbiddenChar( e, Cfg::forbidden_chars_multiphase );
    if ( badchar.has_value() )
      NCRYSTAL_THROW2(BadInput,"Forbidden character "<<badchar.value()
                      <<" in query! Problem found in: \""<<e<<'"');
    if ( e.ltrimmed().startswith('-') )
      NCRYSTAL_THROW2(BadInput,"Query items can not start with a \"-\""
                      " character. Problem found in: \"" << e << '"' );
  }

  const auto& key = query.front();
  constexpr auto sv_version = StrView::make("version");
  constexpr auto sv_list = StrView::make("list");
  constexpr auto sv_mmc = StrView::make("mmc");
  constexpr auto sv_util = StrView::make("util");
  constexpr auto sv_sab = StrView::make("sab");

  if ( key == sv_mmc ) {
    MiniMC::Query::JSONQuery( os, query );
  } else if ( key == sv_list ) {
    if ( query.size() != 1 )
      NCRYSTAL_THROW2(BadInput, "Invalid JSON query (no arguments"
                      " should follow key \"list\")");
    os << "[\"version\", \"util\", \"mmc\", \"sab\"]";
  } else if ( key == sv_version ) {
    queryimpl_version( os, query );
  } else if ( key == sv_util ) {
    queryimpl_util( os, query );
  } else if ( key == sv_sab ) {
    SABUtils::JSONQuery( os, query );
  } else {
    NCRYSTAL_THROW2(BadInput, "Invalid JSON query key: \""<<key<<'"');
  }
}

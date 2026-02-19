#ifndef NCrystal_UtilQuery_hh
#define NCrystal_UtilQuery_hh

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

namespace NCRYSTAL_NAMESPACE {
  namespace {
    void queryimpl_version( std::ostream& os, const Query& query )
    {
      if ( query.size() != 1 || query.front() != "version" )
        NCRYSTAL_THROW2(BadInput, "Invalid JSON query (no arguments"
                        " should follow key \"version\")");
      int vmajor = NCRYSTAL_VERSION_MAJOR;
      int vminor = NCRYSTAL_VERSION_MINOR;
      int vpatch = NCRYSTAL_VERSION_PATCH;
      int vint = NCRYSTAL_VERSION;
      if ( ncgetenv_bool("FIX_QUERY_VERSION_FOR_TESTS") ) {
        vmajor = 123;
        vminor = 45;
        vpatch = 6;
        vint = 1000000*vmajor+1000*vminor+vpatch;
      }
      std::array<int,3> vlist = { vmajor, vminor, vpatch };
      streamJSONDictEntry( os, "major", vmajor, JSONDictPos::FIRST );
      streamJSONDictEntry( os, "minor", vminor );
      streamJSONDictEntry( os, "patch", vpatch );
      streamJSONDictEntry( os, "int", vint );
      streamJSONDictEntry( os, "list", vlist );
      std::ostringstream vstr;
      vstr << vmajor << '.' << vminor << '.' << vpatch;
      streamJSONDictEntry( os, "str", vstr.str(), JSONDictPos::LAST );
    }

    void queryimpl_util( std::ostream& os, const Query& query )
    {
      auto invalid = [&query](const char * reason){
        std::ostringstream ss;
        ss << "Invalid util query: ";
        streamJSON( ss, query );//trick: use streamJSON for easy format.
        if ( reason )
          ss<< " ("<<reason<<')';
        NCRYSTAL_THROW( BadInput, ss.str() );
      };
      constexpr auto sv_util = StrView::make("util");
      if ( query.size() < 2 || query.front() != sv_util ) {
        invalid(nullptr);
        return;
      }
      //shift off the "util" and key entries:
      const auto& key = query.at(1).trimmed();
      auto arg = [&query]( std::size_t i ) { return query.at(i+2); };
      //auto argstr = [&arg]( std::size_t i ) { return arg(i).to_string(); };
      const std::size_t nargs = static_cast<std::size_t>(query.size()-2);

      constexpr auto sv_list = StrView::make("list");
      constexpr auto sv_wl2ekin = StrView::make("wl2ekin");
      constexpr auto sv_ekin2wl = StrView::make("ekin2wl");
      if ( key == sv_list ) {
        if ( nargs != 0 )
          invalid("no arguments should come after: [\"util\",\"list\"]");
        os<<"[\"wl2ekin\", \"ekin2wl\"]";
      } else if ( isOneOf(key,sv_wl2ekin,sv_ekin2wl) ) {
        double val = ( nargs == 1
                       ? arg(0).toDbl().value_or(-1.0)
                       : -1.0 );
        if ( !(val>=0.0) )
          invalid( key==sv_wl2ekin
                   ? "correct usage: [\"util\",\"wl2ekin\",VALUE_ANGSTROM]"
                   : "correct usage: [\"util\",\"ekin2wl\",VALUE_EV]" );
        streamJSON( os, ( key==sv_wl2ekin
                          ? wl2ekin( val )
                          : ekin2wl( val ) ) );
      } else {
        invalid(nullptr);
      }

    }
  }
}

#endif

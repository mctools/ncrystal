#ifndef NCrystal_BrowseQuery_hh
#define NCrystal_BrowseQuery_hh

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


#include "NCrystal/internal/utils/NCStrView.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace BrowseQuery {

    //Implementation of JSON queries for browsing available data:
    //
    //["util","browsedb"]: dict of TextData factory names and their
    //  number of browsable entries.
    //["util","browsedb",FACTNAME,(I,N,)("cheap")]: list with a dict per
    //  entry of the factory (optionally only the I'th of N contiguous
    //  chunks). Unless "cheap", each entry is loaded with createInfo
    //  (in parallel if factory threads are enabled), and physics
    //  properties (or a load error) are included.
    //["util","browsefactories"]: names of all factories, by type.
    //
    //The args are the query items following the "browsedb" key.
    void browseDB( std::ostream&, const std::vector<StrView>& args );
    void browseFactories( std::ostream& );

  }

}

#endif

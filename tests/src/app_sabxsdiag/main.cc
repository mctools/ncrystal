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

//SABXSDIAG-TEMPORARY: dedicated, easily revertible diagnostic commit. Dumps
//full-precision intermediate SAB grid/table state for the one known
//cross-platform-unstable case (Li_from_Li2O.ncmat;vdoslux=2004;knllux=4;
//temp=10, matching tests/scripts/sabxs.py), so a real CI run on the two
//actually-diverging legs can be diffed against a passing leg. Not meant to
//stay in the codebase: revert this whole commit once the data has been
//captured.

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCTestUtils/NCTestFindData.hh"
#include <cstdlib>
#include <iomanip>

namespace NC = NCrystal;

int main()
{
#if !defined(_WIN32)
  setenv("NCRYSTAL_SABXS_DIAG","1",1);
#endif

  auto path = nctest::find_test_data("Li_from_Li2O.ncmat");
  NC::MatCfg cfg( path + ";vdoslux=2004;knllux=4;temp=10" );
  auto sc = NC::FactImpl::createScatter(cfg);
  NC::CachePtr cp;
  for ( auto e : { 0.02654279, 0.04641589 } ) {
    auto xs = sc->crossSectionIsotropic( cp, NC::NeutronEnergy{e} ).get();
    NCRYSTAL_MSG( "SABXSDIAG xsect(" << std::setprecision(17) << e
                 << ")=" << std::setprecision(17) << xs );
  }
  return 0;
}

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

////////////////////////////////////////////////////////////////////////////////
// Self-timed benchmark of VDOS->SAB construction (expandVDOSToGnFcts +        //
// FastConvolve, the hot path chased across many commits fixing a              //
// Windows-only crash) for a small, fixed set of materials. Deliberately has   //
// no test.log reference: timings vary run to run (CI runner load, platform,  //
// compiler), so this is registered as a plain pass/fail (exit code) ctest     //
// test rather than one whose *output* is diffed against a fixed reference --  //
// what matters is the timing lines printed to stdout, meant to be read from  //
// the ctest log capture (e.g. in CI), not compared automatically.            //
//                                                                            //
// Each material's SAB data is (re-)extracted from scratch nreps times        //
// (useCache=false), timed individually with std::chrono::steady_clock, and   //
// reported as min/median/max -- min is the most representative of steady-    //
// state cost (least affected by transient CI-runner noise), median guards    //
// against a single lucky outlier, max flags unusually high jitter.           //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/factories/NCMatCfg.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/dyninfoutils/NCDynInfoUtils.hh"
#include <iostream>
#include <chrono>
#include <algorithm>
#include <cstdio>

namespace NC = NCrystal;

namespace {

  double toMS( std::chrono::steady_clock::duration d )
  {
    return std::chrono::duration<double,std::milli>(d).count();
  }

  //One material's worth of DI_ScatKnl entries, each timed nreps times:
  void benchMaterial( const std::string& cfgstr, unsigned nreps )
  {
    NC::MatCfg matcfg( cfgstr );
    const NC::VDOS::VDOSLux vdoslux( matcfg.get_vdoslux() );
    auto info = NC::FactImpl::createInfo( cfgstr );

    for ( auto& di : info->getDynamicInfoList() ) {
      auto di_scatknl = dynamic_cast<const NC::DI_ScatKnl*>( di.get() );
      if ( !di_scatknl )
        continue;//not a phonon/VDOS-based component (e.g. free gas, sterile)

      std::vector<double> times_ms;
      times_ms.reserve( nreps );
      for ( unsigned i = 0; i < nreps; ++i ) {
        const auto t0 = std::chrono::steady_clock::now();
        auto sab = NC::extractSABDataFromDynInfo( di_scatknl, vdoslux,
                                                  false/*useCache*/ );
        const auto t1 = std::chrono::steady_clock::now();
        nc_assert_always( sab != nullptr );
        times_ms.push_back( toMS( t1 - t0 ) );
      }
      std::sort( times_ms.begin(), times_ms.end() );
      const double tmin = times_ms.front();
      const double tmax = times_ms.back();
      const double tmedian = times_ms[ times_ms.size() / 2 ];
      std::printf( "PERFVDOS: material=\"%s\" element=\"%s\" nreps=%u"
                   " min_ms=%.3f median_ms=%.3f max_ms=%.3f\n",
                   cfgstr.c_str(), di->atomData().elementName().c_str(),
                   nreps, tmin, tmedian, tmax );
    }
  }

}

int main()
{
  //Deliberately small, fixed set: enough phonon orders to exercise the
  //VDOS Gn expansion/FastConvolve hot path meaningfully (Al_sg225 alone
  //reaches phonon order ~168 at the default vdoslux), while keeping total
  //runtime modest even on a slower CI runner or unoptimised build:
  const std::vector<std::string> materials = {
    "Al_sg225.ncmat",
    "Li2O_sg225_LithiumOxide.ncmat",
  };
  constexpr unsigned nreps = 5;

  for ( auto& cfgstr : materials )
    benchMaterial( cfgstr, nreps );

  return 0;
}

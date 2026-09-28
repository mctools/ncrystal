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
// Self-timed benchmark of the *simulation-time* SAB cross-section/scattering  //
// hot path, as opposed to app_perfvdos which only times one-time VDOS->SAB    //
// kernel *construction* at material load. Uses a vdoslux>=2000 ("next-gen")   //
// material, whose crossSectionIsotropic/sampleScatterIsotropic calls run     //
// through NCSABScatterNG.cc -> NCSABProcessor.cc -> NCSABCellInteg.cc's       //
// StdLogLinCellIntegrator, exercising exactly the on-demand SAB cell-        //
// integration cost paid throughout an actual simulation, not just at setup.  //
// Sibling of the developer-only, ad hoc devel/simplebuild/static/NCDev/      //
// app_benchsabxs/app_benchsabsample (simplebuild-only, not CI/Windows-       //
// buildable); this app instead follows app_perfvdos's own                    //
// CMake-buildable/self-timed/no-reference-log convention, so it can run in   //
// benchmark.yml across platforms including Windows.                         //
//                                                                            //
// The energy is cycled log-uniformly across a wide range for every call,    //
// rather than reusing one fixed energy throughout: a repeated energy always //
// lands in the same grid cell, which (after the first call) reuses a cached //
// cell instead of re-running the actual cell-integration math -- see        //
// app_benchsabxs/app_benchsabsample's own comments for the same point.      //
////////////////////////////////////////////////////////////////////////////////

#include "NCrystal/factories/NCFactImpl.hh"
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

  double cycledEkin( double ekin_lo, double logratio,
                     unsigned i, unsigned nsample )
  {
    return ekin_lo * std::exp( logratio * ( double(i) / nsample ) );
  }

  void reportTimes( const char* mode, const std::string& cfgstr,
                    unsigned nsample, unsigned nreps,
                    std::vector<double>& times_ms )
  {
    std::sort( times_ms.begin(), times_ms.end() );
    std::printf( "PERFSAB: mode=\"%s\" material=\"%s\" nsample=%u nreps=%u"
                 " min_ms=%.3f median_ms=%.3f max_ms=%.3f\n",
                 mode, cfgstr.c_str(), nsample, nreps,
                 times_ms.front(), times_ms[ times_ms.size() / 2 ],
                 times_ms.back() );
  }

  void benchXS( const std::string& cfgstr, double ekin_lo, double ekin_hi,
               unsigned nsample, unsigned nreps )
  {
    std::printf( "PERFSAB: mode=\"xs\" material=\"%s\" starting %u reps\n",
                cfgstr.c_str(), nreps );
    auto scatter = NC::FactImpl::createScatter( cfgstr );
    const double logratio = std::log( ekin_hi / ekin_lo );
    std::vector<double> times_ms;
    times_ms.reserve( nreps );
    double xssum = 0.0;//accumulate to prevent the loop being optimised away
    for ( unsigned r = 0; r < nreps; ++r ) {
      NC::CachePtr cache;
      const auto t0 = std::chrono::steady_clock::now();
      for ( unsigned i = 0; i < nsample; ++i ) {
        const NC::NeutronEnergy ekin{ NC::DoValidate,
          cycledEkin(ekin_lo,logratio,i,nsample) };
        xssum += scatter->crossSectionIsotropic( cache, ekin ).dbl();
      }
      const auto t1 = std::chrono::steady_clock::now();
      times_ms.push_back( toMS( t1 - t0 ) );
    }
    (void)xssum;
    reportTimes( "xs", cfgstr, nsample, nreps, times_ms );
  }

  void benchSample( const std::string& cfgstr, double ekin_lo, double ekin_hi,
                    unsigned nsample, unsigned nreps )
  {
    std::printf( "PERFSAB: mode=\"sample\" material=\"%s\" starting %u reps\n",
                cfgstr.c_str(), nreps );
    auto scatter = NC::FactImpl::createScatter( cfgstr );
    auto rng = NC::getRNG();
    const double logratio = std::log( ekin_hi / ekin_lo );
    std::vector<double> times_ms;
    times_ms.reserve( nreps );
    for ( unsigned r = 0; r < nreps; ++r ) {
      NC::CachePtr cache;
      const auto t0 = std::chrono::steady_clock::now();
      for ( unsigned i = 0; i < nsample; ++i ) {
        const NC::NeutronEnergy ekin{ NC::DoValidate,
          cycledEkin(ekin_lo,logratio,i,nsample) };
        (void)scatter->sampleScatterIsotropic( cache, rng, ekin );
      }
      const auto t1 = std::chrono::steady_clock::now();
      times_ms.push_back( toMS( t1 - t0 ) );
    }
    reportTimes( "sample", cfgstr, nsample, nreps, times_ms );
  }

}

int main()
{
  //Unbuffer stdout -- see app_perfvdos's own top comment for why this
  //matters on Windows specifically (lost output on a hard crash/hang):
  std::setvbuf(stdout, nullptr, _IONBF, 0);

  //vdoslux=2000 selects the "next-gen" SAB expansion, whose cross-section/
  //sampling path runs through NCSABCellInteg.cc's StdLogLinCellIntegrator:
  const std::string cfgstr = "Al_sg225.ncmat;vdoslux=2000";
  constexpr double ekin_lo = 1e-5;//eV
  constexpr double ekin_hi = 5.0;//eV
  constexpr unsigned nreps = 5;

  benchXS( cfgstr, ekin_lo, ekin_hi, 1000000, nreps );
  benchSample( cfgstr, ekin_lo, ekin_hi, 400000, nreps );

  std::printf( "BENCH: all done\n" );
  return 0;
}

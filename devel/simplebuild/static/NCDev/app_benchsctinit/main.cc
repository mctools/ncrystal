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

//Benchmark of SCTXSProvider construction cost (i.e. the fixed-order
//quadrature + PCHIP lookup-table initialisation), intended for profiling
//with e.g. valgrind/callgrind. Usage:
//
//   sb_ncdev_benchsctinit [amu r nctor]
//
//With no arguments a small default matrix is run.

#include "NCrystal/internal/phys_utils/NCSCTUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include <chrono>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>

namespace NC = NCrystal;

namespace {
  unsigned envNPts()
  {
    //The production node count is tied to knllux via
    //SABCfg::Cfg::sct_table_npts; this dev-only override supports
    //accuracy/cost studies:
    return static_cast<unsigned>(
      NC::ncgetenv_int("SCTLUTNPTS", 90 ) );
  }

  double benchCtor( double amu, double r, unsigned nctor )
  {
    const NC::Temperature T{293.15};
    const NC::Temperature Teff{ T.dbl() * r };
    const NC::AtomMass mass{amu};
    NC::StableSum dump;
    auto t0 = std::chrono::steady_clock::now();
    for ( unsigned i = 0; i < nctor; ++i ) {
      NC::SCTXSProvider p( T, Teff, mass, NC::SigmaBound{1.0},
                           envNPts() );
      dump.add( p.crossSection( NC::NeutronEnergy{1.0} ).dbl() );
    }
    auto dt = std::chrono::duration<double,std::milli>(
      std::chrono::steady_clock::now() - t0 ).count();
    std::printf("BENCHSCTINIT amu=%g r=%g nctor=%u : %.3f ms/ctor"
                " (checksum %g)\n", amu, r, nctor, dt / nctor,
                dump.sum() );
    if ( NC::ncgetenv_bool("SCTLUT_PROBE") ) {
      //Accuracy probe: densely compare the production lookup against
      //accurate direct quadrature:
      NC::SCTXSProvider p( T, Teff, mass, NC::SigmaBound{1.0},
                           envNPts() );
      const double A = mass.relativeToNeutronMass();
      const double csat = 35.0 * ( A > 1.0 ? A : 1.0/A );
      double worst = 0.0, worst_at = 0.0;
      constexpr unsigned nprobe = 300;
      for ( unsigned i = 0; i < nprobe; ++i ) {
        const double c = 1.1e-3 * std::pow( 0.999*csat/1.1e-3,
                                            i/(nprobe-1.0) );
        const double v1 = p.evalUpscatterCorr( c );
        const double v2 = p.evalUpscatterCorrIntegral( c, 1e-10 );
        if ( v2 > 1e-10 ) {
          const double e = std::fabs( v1/v2 - 1.0 );
          if ( e > worst ) { worst = e; worst_at = c; }
        }
      }
      std::printf("  PROBE worst reldiff %.3g at c=%.4g (c/csat=%.3f)\n",
                  worst, worst_at, worst_at/csat );
    }
    return dt;
  }
}

namespace {
  void tsweepProbe( double amu, double e_ev )
  {
    //Smoothness probe: sweep temperature smoothly (with a smooth
    //synthetic Teff(T) model with a zero-point-like floor) and look
    //for kinks/jumps in the resulting cross section curve via scaled
    //second differences:
    const NC::AtomMass mass{amu};
    const unsigned n = static_cast<unsigned>(
      NC::ncgetenv_int("SCTLUT_TSWEEP_N", 2000 ) );
    constexpr double t0 = 150.0, t1 = 900.0;
    constexpr double teff_floor = 1100.0;
    std::vector<double> xs;
    xs.reserve(n);
    for ( unsigned i = 0; i < n; ++i ) {
      const double T = t0 + (t1-t0) * i / (n-1.0);
      const double teff = std::sqrt( T*T + teff_floor*teff_floor );
      NC::SCTXSProvider p( NC::Temperature{T}, NC::Temperature{teff},
                           mass, NC::SigmaBound{1.0}, envNPts() );
      xs.push_back( p.crossSection( NC::NeutronEnergy{e_ev} ).dbl() );
    }
    //Roughness metric: second differences normalised to the local
    //value. For a smooth curve this is ~(dT)^2*f''/f ~ 1e-8; kinks
    //and jumps stand out far above that:
    double worst = 0.0, worst_T = 0.0;
    NC::StableSum sum;
    std::vector<double> d2s;
    d2s.reserve(n);
    for ( unsigned i = 1; i+1 < n; ++i ) {
      const double d2 = std::fabs( xs[i+1] - 2.0*xs[i] + xs[i-1] )
        / xs[i];
      d2s.push_back( d2 );
      sum.add( d2 );
      if ( d2 > worst ) {
        worst = d2;
        worst_T = t0 + (t1-t0)*i/(n-1.0);
      }
    }
    std::sort( d2s.begin(), d2s.end() );
    const double median = d2s[ d2s.size()/2 ];
    std::printf("  TSWEEP amu=%g E=%geV n=%u : worst |d2 xs|/xs %.3g"
                " at T=%.1fK (median %.3g, worst/median %.1f)\n",
                amu, e_ev, n, worst, worst_T, median,
                worst / ( median > 0.0 ? median : 1e-300 ) );
  }
}

int main( int argc, char** argv )
{
  if ( argc == 4 ) {
    benchCtor( std::atof(argv[1]), std::atof(argv[2]),
               unsigned(std::atol(argv[3])) );
    return 0;
  }
  if ( argc != 1 ) {
    std::printf("Please provide args: [amu r nctor] (or none for"
                " default matrix)\n");
    return 1;
  }
  for ( double amu : { 1.008, 27.0, 207.2 } )
    for ( double r : { 1.05, 4.12, 59.3 } )
      benchCtor( amu, r, 50 );
  if ( NC::ncgetenv_bool("SCTLUT_TSWEEP") ) {
    for ( double amu : { 1.008, 207.2 } )
      for ( double e_ev : { 1.5, 0.6 } )
        tsweepProbe( amu, e_ev );
  }
  return 0;
}

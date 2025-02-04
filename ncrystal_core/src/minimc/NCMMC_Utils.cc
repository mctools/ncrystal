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

#include "NCrystal/internal/minimc/NCMMC_Utils.hh"
#include "NCrystal/internal/extd_utils/NCABIUtils.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"

namespace NC = NCrystal;
namespace NCMMCU = NCrystal::MiniMC::Utils;

void NCMMCU::calcProbTransm( NumberDensity nd, std::size_t N,
                             const double * ncrestrict xs_or_nullptr,
                             const double * ncrestrict dist,
                             double * ncrestrict out )
{
  if ( !xs_or_nullptr ) {
    for ( auto i : ncrange( N ) )
      out[i] = 1.0;
    return;
  }
  for ( auto i : ncrange( N ) )
    out[i] = macroXS( nd, CrossSect{xs_or_nullptr[i]} );
  for ( auto i : ncrange( N ) )
    out[i] *= dist[i];
  for ( auto i : ncrange( N ) )
    out[i] = -out[i];
  for ( auto i : ncrange( N ) )
    out[i] = std::exp( out[i] );
}

void NCMMCU::propagate( NeutronBasket& b,
                        const double* ncrestrict dists )
{
  for ( auto i : ncrange( b.size() ) )
    b.x[i] += dists[i] * b.ux[i];
  for ( auto i : ncrange( b.size() ) )
    b.y[i] += dists[i] * b.uy[i];
  for ( auto i : ncrange( b.size() ) )
    b.z[i] += dists[i] * b.uz[i];
}

void NCMMCU::propagateAndAttenuate( NeutronBasket& b,
                                    NumberDensity nd,
                                    const double* ncrestrict dists,
                                    const double* ncrestrict xsvals )
{
  propagate( b, dists );
  if ( xsvals ) {
    for ( auto i : ncrange( b.size() ) )
      b.w[i] *= std::exp( -macroXS( nd, CrossSect{ xsvals[i] } ) * dists[i] );
  }
}

void NCMMCU::sampleRandDists( RNG& rng, NumberDensity nd,
                              const double * ncrestrict dists,
                              const double * ncrestrict xsvals,
                              std::size_t N,
                              double * ncrestrict tgt )
{
  NewABI::generateMany( rng, N, tgt );
  //TODO: Better vectorisation possible?
  for ( auto i : ncrange(N) ) {
    RandExpIntervalSampler rs( 0.0, dists[i],
                               macroXS( nd, CrossSect{ xsvals[i] } ) );
    tgt[i] = rs.sample( tgt[i] );
  }
}

void NCMMCU::scatterGivenMu( RNG& rng,
                             NeutronBasket& b,
                             double * ncrestrict mu_vals )
{
  for ( auto i : ncrange(b.size()) ) {
    //TODO: better vectorisation possible?
    auto newdir
      = randNeutronDirectionGivenScatterMu( rng,
                                            CosineScatAngle{ mu_vals[i] },
                                            NeutronDirection{ b.ux[i],
                                                              b.uy[i],
                                                              b.uz[i] } );
    b.ux[i] = newdir[0];
    b.uy[i] = newdir[1];
    b.uz[i] = newdir[2];
  }
}

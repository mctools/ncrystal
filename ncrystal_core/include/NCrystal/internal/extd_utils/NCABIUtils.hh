#ifndef NCrystal_ABIUtils_hh
#define NCrystal_ABIUtils_hh

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

#include "NCrystal/core/NCDefs.hh"
#include "NCrystal/interfaces/NCProcImpl.hh"

namespace NCRYSTAL_NAMESPACE {

  // ABI-breaking changes to the public NCrystal ABI should be rare and
  // well-tested before being introduced, and should induce a major version
  // number bump of the NCrystal release. To help with this procedure, we should
  // ideally introduce ABI-breaking changes (in particular new virtual methods
  // on interface base classes) protected with NCRYSTAL_ALLOW_ABI_BREAKAGE, and
  // use functions in the present header file, to these new ABI-breaking methods
  // or their ABI-preserving fall-back implementations, as appropriate for a
  // given NCrystal build. This allows us to take advantage of the new features
  // in development builds, but to hopefully be able to postpone the actual ABI
  // breakage to a point when it is most suitable.

  namespace NewABI {
    inline void generateMany( RNG& rng, std::size_t n, double* tgt )
    {
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
      rng.generateMany(n,tgt);
#else
      for ( std::size_t i = 0; i < n; ++i )
        *tgt++ = rng.generate();
#endif
    }
  }

  namespace ProcImpl {

    namespace NewABI {
      inline bool isPureElasticScatter( const Process& p )
      {
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
        return p.isPureElasticScatter();
#else
        auto proc_comp = dynamic_cast<const ProcImpl::ProcComposition*>(&p);
        if ( proc_comp ) {
          for ( auto& comp : proc_comp->components() ) {
            if ( !isPureElasticScatter( comp.process ) )
              return false;
          }
          return true;
        } else {
          std::string name(p.name());
          return ( name == "PowderBragg"
                   || name == "LCBragg"
                   || name == "SCBragg"
                   || name == "ElIncScatter"
                   || name == "NullScatter" );
        }
#endif
      }

      //If you know you need both the cross section and to sample a scatterings,
      //it is usually more efficient to request both at once:
      inline std::pair<CrossSect,ScatterOutcome>
      evalXSAndSampleScatter( const Process& p,
                              CachePtr& cp, RNG& rng,
                              NeutronEnergy ekin, const NeutronDirection& ndir )
      {
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
        return p.evalXSAndSampleScatter(cp,rng,ekin,ndir);
#else
        auto xs = p.crossSection( cp, ekin, ndir );
        return { xs, p.sampleScatter( cp, rng, ekin, ndir ) };
#endif
      }

      inline std::pair<CrossSect,ScatterOutcomeIsotropic>
      evalXSAndSampleScatterIsotropic(const Process& p, CachePtr& cp, RNG& rng, NeutronEnergy ekin )
      {
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
        return p.evalXSAndSampleScatterIsotropic( cp, rng, ekin );
#else
        auto xs = p.crossSectionIsotropic( cp, ekin );
        return { xs, p.sampleScatterIsotropic( cp, rng, ekin ) };
#endif
      }

      inline void evalManyXS( const Process& p, CachePtr& cp, const double* ekin,
                              const double* ux, const double* uy, const double* uz,
                              std::size_t N,
                              double* out_xs )
      {
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
        return p.evalManyXS( cp, ekin, ux, uy, uz, N, out_xs );
#else
        for ( auto i: ncrange(N) ) {
          (void)i;
          *out_xs++ = p.crossSection( cp,
                                      NeutronEnergy{*ekin++},
                                      NeutronDirection{ *ux++, *uy++, *uz++ } ).dbl();
        }
#endif
      }

      inline void evalManyXSIsotropic( const Process& p, CachePtr& cp, const double* ekin, std::size_t N,
                                       double* out_xs )
      {
#ifdef NCRYSTAL_ALLOW_ABI_BREAKAGE
        return p.evalManyXSIsotropic( cp, ekin, N, out_xs );
#else
        for ( auto i : ncrange(N) ) {
          (void)i;
          *out_xs++ = p.crossSectionIsotropic( cp, NeutronEnergy{*ekin++} ).dbl();
        }
#endif
      }
    }
  }
}

#endif

#ifndef NCrystal_SABProcessor_hh
#define NCrystal_SABProcessor_hh

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

#include "NCrystal/interfaces/NCSABData.hh"

namespace NCRYSTAL_NAMESPACE {

  namespace SABUtils {

    class SABProcessor final : private MoveOnly {

      // Class which processes an S(alpha,beta) table and provides both cross
      // sections and sampling capabilities. In addition to an understanding of
      // which cells are touched or covered by kinematic bounds of various
      // energies (provided by SABSurveyor), this crucially depends on
      // establishing a grid of energy values. At each energy point in this
      // grid, the S(alpha,beta) table is integrated in order to provide both
      // data needed for cross-sections, as well as sufficient per-cell integral
      // information for an efficient sampling of scattering events later.

      // fixme: more details here or elsewhere, including about egrid, algs,
      // etc.

    public:

      ~SABProcessor();

      enum class SampleSupport { YES, NO };
      enum class StoreExtraDiagnostics { YES, NO };
      SABProcessor( shared_obj<const SABData>,
                    std::shared_ptr<const VectD> egrid,
                    SampleSupport = SampleSupport::YES,
                    StoreExtraDiagnostics = StoreExtraDiagnostics::NO );

      //integral of S(alpha,beta) within phasespace of a given neutron energy:
      double phaseSpaceIntegral( NeutronEnergy ) const;

      //Sample a scattering event in (dE,mu) or (alpha,beta) space, with or
      //without diagnostics:
      ScatterOutcomeIsotropic sampleScatter( RNG&, NeutronEnergy ) const;

      struct AlphaBetaOutcome {
        double alpha, beta;
        unsigned ntries;
      };
      AlphaBetaOutcome sampleScatterAlphaBeta( RNG&, NeutronEnergy ) const;

      struct ScatOutcomeDiag {
        ScatterOutcomeIsotropic ekinmu;
        AlphaBetaOutcome alphabeta;
      };
      ScatOutcomeDiag sampleScatterDiag( RNG&, NeutronEnergy ) const;


      bool hasSampleSupport() const;

      //fixme: also a streamJSONSummary?
      void toJSON( std::ostream& ) const;

      //Maximal energy, after which the SAB does not reliably cover neutron
      //physics (this is the last point in the internal energy grid):
      NeutronEnergy getEMax() const;

      //For reference, access the energy grid (in units of kT):
      const VectD& getEDivKTGrid() const;

      //Move-only:
      SABProcessor( const SABProcessor& ) = delete;
      SABProcessor& operator=( const SABProcessor& ) = delete;
      SABProcessor( SABProcessor&& ) noexcept;
      SABProcessor& operator=( SABProcessor&& ) noexcept;

    private:
      void * m_impl = nullptr;
    };
  }
}

#endif

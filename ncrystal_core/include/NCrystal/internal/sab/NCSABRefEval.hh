#ifndef NCrystal_SABRefEval_hh
#define NCrystal_SABRefEval_hh

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

#include "NCrystal/internal/sab/NCSABEval.hh"
#include "NCrystal/interfaces/NCSABData.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/internal/utils/NCStrView.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {

    //Slow but extremely trustworthy reference evaluation of cross sections and
    //samplings of SABData at a given neutron energy point.

    template< InterpolationScheme alphaInterpType_ = InterpolationScheme::LOGLIN,
              SABInterpolationOrder interpOrder_ = SABInterpolationOrder::ALPHA_FIRST>
    class SABRefEval final : public NoCopyMove {
    public:
      using celleval_t = SABCellEval<alphaInterpType_,interpOrder_>;
      using sabeval_t = SABEval<alphaInterpType_,interpOrder_>;
      using cell_index_t = PackedCellIndex::index_t;
    private:
      sabeval_t m_sabeval;
      VectD m_cell_integrals_commul;
      std::vector<cell_index_t> m_cell_idx;
      NeutronEnergy m_ekin;
      Temperature m_temp;
      double m_kT;
      double m_ekin_div_kT = -1.0;
      double m_acceptanceRate = -1.0;
      double m_xsect_val = -1.0;
    public:
      SABRefEval( shared_obj<const SABData> sabdata, NeutronEnergy eval )
        : m_sabeval( sabdata ),
          m_ekin( eval ),
          m_temp(sabdata->temperature()),
          m_kT(sabdata->temperature().kT()),
          m_ekin_div_kT( eval.get() / m_kT )
      {
        nc_assert_always(sabdata!=nullptr);
        std::vector<std::pair<double,cell_index_t>> cell_contribs;
        StableSum touched_cells_total_integral;
        for ( auto cellidx : ncrange( m_sabeval.nRawCellIdxMax() ) ) {
          auto cell = m_sabeval.getCell( cellidx );
          double contrib = cell.integralWithinKinematicBounds( m_ekin_div_kT );
          if (contrib) {
            cell_contribs.emplace_back(contrib,cellidx);
            touched_cells_total_integral.add( cell.integral() );
          }
        }
        std::sort( cell_contribs.begin(), cell_contribs.end() );
        StableSum contrib_commul;
        for ( auto cc : cell_contribs ) {
          contrib_commul.add( cc.first );
          m_cell_integrals_commul.push_back(contrib_commul.sum());
          m_cell_idx.push_back(cc.second);
        }

        m_acceptanceRate = contrib_commul.sum()/touched_cells_total_integral.sum();
        const double xs_scale_fact
          = 0.25 * sabdata->boundXS().get() * sabdata->temperature().kT();
        m_xsect_val = xs_scale_fact * contrib_commul.sum() / m_ekin.get();
      }

      double sabIntegral() const { return m_cell_integrals_commul.back(); };
      CrossSect xsect() const { return CrossSect{ m_xsect_val }; }
      double samplingAcceptanceRate() const { return m_acceptanceRate; }

      struct SampleResult {
        PairDD alphaBeta;
        DeltaEMuVal_t dEMu;
      };

      SampleResult sampleScatter( RNG& rng ) const
      {
        SampleResult result;
        auto cell_contrib_nbr = pickRandIdxByWeight( rng, m_cell_integrals_commul );
        auto cell = m_sabeval.getCell( m_cell_idx.at(cell_contrib_nbr) );
        while ( true ) {
          result.alphaBeta = cell.sample(rng);
          if ( result.alphaBeta.second < -m_ekin_div_kT )
            continue;
          auto alims = getAlphaLimits( m_ekin_div_kT, result.alphaBeta.second );
          if ( valueInInterval( alims.first, alims.second, result.alphaBeta.first ) )
            break;
        }
        if ( muIsotropicAtBeta( result.alphaBeta.second, m_ekin_div_kT ) ) {
          result.dEMu.deltaE = ncmax(-m_ekin.get(),result.alphaBeta.second*m_kT);
          result.dEMu.mu = 2.0*rng()-1.0;
        } else {
          result.dEMu = convertAlphaBetaToDeltaEMu( result.alphaBeta.first,
                                                    result.alphaBeta.second,
                                                    m_ekin, m_kT );
        }
        return result;
      }

      void toJSON( std::ostream& os,
                   std::uint64_t nsample = 0,
                   RNG* rng = nullptr ) const
      {
        if ( nsample && !rng )
          NCRYSTAL_THROW(BadInput,"SABRefEval::toJSON needs RNG if nsample>0");
        os << "{\"neutron_energy\":";
        streamJSON(os,m_ekin.get());
        os << ",\"temperature\":";
        streamJSON(os,m_temp.get());
        os << ",\"kT\":";
        streamJSON(os,m_kT);
        os << ",\"sab_phasespace_integral\":";
        streamJSON(os,sabIntegral());
        os << ",\"cross_section\":";
        streamJSON(os,xsect().get());
        os << ",\"sampling_acceptance_rate\":";
        streamJSON(os,m_acceptanceRate);
        os << ",\"samples_format\":[\"alpha\",\"beta\",\"deltaE\",\"mu\"]";
        os << ",\"samples\":[";
        for ( std::uint64_t i = 0; i < nsample; ++i ) {
          auto res = this->sampleScatter( *rng );
          if ( i )
            os << ',';
          os << "[";
          streamJSON(os,res.alphaBeta.first);
          os << ',';
          streamJSON(os,res.alphaBeta.second);
          os << ',';
          streamJSON(os,res.dEMu.deltaE);
          os << ',';
          streamJSON(os,res.dEMu.mu);
          os << ']';
        }
        os << "]}";
      }

    };

  }
}

#endif

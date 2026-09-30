#ifndef NCrystal_SABAnalyser_hh
#define NCrystal_SABAnalyser_hh

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

////////////////////////////////////////////////////////////////////////////////
// Estimation of the effective temperature (Teff, as needed by the SCT       //
// kernel-extension model) and the mean-squared displacement directly from   //
// a tabulated S(alpha,beta) kernel, via exact moment sum rules of the       //
// incoherent/Gaussian approximation (valid at every alpha, so robust        //
// median-combined per-row estimators result):                               //
//                                                                           //
//   Teff: Wick's second-moment rule on the asymmetric kernel,               //
//         Var_beta(alpha) = 2*alpha*(m_n/M)*(Teff/T)                        //
//         (NCrystal's alpha convention is alpha=hbar^2*q^2/(2*m_n*kT),      //
//         cf. convertAlphaBetaToDeltaEMu, i.e. A times the ENDF alpha).     //
//   msd:  the Debye-Waller zeroth-moment rule for kernels with the          //
//         elastic line excluded, N(alpha) = 1 - exp(-c*alpha) with          //
//         c = (2*m_n*kT/hbar^2)*msd (exactly NCVDOSExpand's alpha2x).       //
//         Kernels without a Debye-Waller deficit (e.g. liquids, where       //
//         N(alpha)~=1) correctly yield "not inferable".                     //
//                                                                           //
// Guards per row: kernels store tiny floor-marker placeholder values in     //
// far wings (masked away), the +-coverage_nsigma wings must be covered by   //
// real data, and the recoil ridge must be resolved by the grid. Estimates   //
// are returned raw together with robustness metrics -- acceptance policy    //
// (spread thresholds etc.) deliberately belongs to the consumers, so the    //
// tuning phase can explore cuts without rebuilding. Validated against       //
// VDOSEval-derived truth to <1% (Teff) and <5% (msd) on VDOS-expanded       //
// kernels across A=1..56 and T=20..800K.                                    //
////////////////////////////////////////////////////////////////////////////////

namespace NCRYSTAL_NAMESPACE {

  namespace SABAnalyser {

    struct Options {
      //NB: the per-row validity gates (floor-marker detection, wing
      //coverage, ridge resolution, msd m0 window, ...) are hardcoded
      //constants in the .cc: their values follow from error-propagation
      //or convergence arguments rather than tuning, as confirmed by a
      //sensitivity scan (each varied x2 up/down over 263 ENDF8-, VDOS-
      //and tests/data-derived kernels: accepted results flat at <<0.1%).
      //The single genuinely data-derived value is center_tol:
      //
      //Per-row self-validation via the first-moment (recoil) sum rule:
      //only rows whose measured mean matches the expected recoil center
      //-alpha*m_n/M to within this relative tolerance contribute to
      //Teff. This certifies the row's sum-rule integrity: stored (e.g.
      //ENDF-converted) kernels break the delicate detailed-balance
      //cancellation at small alpha (finite precision/grids), where the
      //estimators otherwise produce garbage amplified by 1/alpha, while
      //the large-alpha deep-recoil rows -- exactly the SCT regime --
      //plateau at the true value. Genuinely non-Gaussian kernels (e.g.
      //cryogenic ortho-H) have no conforming rows at all and correctly
      //yield "not inferable":
      double center_tol = 0.05;
      bool collect_diagnostics = false;
    };

    //Why a given row did (not) contribute (values are stable, used as
    //integers in the JSON query output):
    enum class RowStatus : std::uint8_t { Ok = 0,
                                          AlphaTooLow = 1,
                                          NoData = 2,
                                          WingsNotCovered = 3,
                                          Unresolved = 4,
                                          M0OutOfWindow = 5,
                                          CenterOffRecoil = 6 };

    struct Diagnostics {
      //Per tabulated alpha row (all vectors aligned with alpha; the
      //moment entries are -1 where not computable), for inspection and
      //tuning scripts:
      VectD alpha;
      VectD m0;       //zeroth beta-moments N(alpha) (floor-masked data)
      VectD mean;     //first-moment row means
      VectD variance; //central second moments
      VectD teff_row; //per-row Teff estimate in K (-1 where inadmissible)
      VectD msd_row;  //per-row msd estimate in Aa^2 (-1 where inadmissible)
      std::vector<RowStatus> teff_row_status;
      std::vector<RowStatus> msd_row_status;
      double floor_value = -1.0;
    };

    struct Result {
      //Raw estimates: absent when no admissible rows. NB: deliberately
      //no acceptance policy applied here; consumers should threshold on
      //the relspread metrics (IQR/median over the per-row estimates)
      //and validate the values. Recommended acceptance for automatic
      //Teff usage, tuned on 271 ENDF8-converted kernels with recorded
      //ENDF effective temperatures: teff_nrows>=20 and
      //teff_relspread<=0.05, giving worst-case 0.8% and p95 0.34%
      //there (kernels failing the policy keep the free-gas fallback):
      Optional<Temperature> teff;
      Optional<double> msd;//[Aa^2]
      double teff_relspread = -1.0;
      double msd_relspread = -1.0;
      unsigned teff_nrows = 0;
      unsigned msd_nrows = 0;
      //Consistency diagnostic: median over admissible rows of the row
      //mean divided by the expected recoil center -alpha*m_n/M (should
      //be ~1 for kernels obeying the Gaussian-approximation sum rules):
      Optional<double> recoil_center_ratio;
      //Set iff Options::collect_diagnostics:
      std::shared_ptr<const Diagnostics> diagnostics;
    };

    Result analyse( const SABData&, const Options& = Options() );

  }
}

#endif

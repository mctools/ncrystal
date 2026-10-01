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

#include "NCrystal/internal/extd_utils/NCSABAnalyser.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/core/NCMem.hh"
#include <mutex>
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace SABAnalyser {
    namespace {

      //Per-row validity gates. These are NOT tuning parameters: each
      //value follows from an error-propagation or convergence argument
      //(and the floor pair from how kernel files encode placeholder
      //floors), and a sensitivity scan varying each x2 up/down (orders
      //of magnitude for the floors) over 263 ENDF8-, VDOS- and
      //tests/data-derived kernels moved accepted results by <<0.1%
      //(see nocheck/ana_scan.py on the tk_volatile branch):
      constexpr double opt_alpha_min = 0.1;       //1/alpha blowup guard
      constexpr double opt_coverage_nsigma = 4.0; //+-4sigma trunc. bias ~1e-4
      constexpr unsigned opt_resolution_minpts = 8;//points per +-2sigma
      constexpr double opt_floor_relmin = 10.0;   //S<=relmin*min(S>0) or...
      constexpr double opt_floor_relmax = 1e-30;  //...S<=relmax*max(S)
      constexpr double opt_msd_m0_low = 0.05;     //1/N and 1/(1-N) error
      constexpr double opt_msd_m0_high = 0.95;    //  amplification guards
      constexpr unsigned opt_n_sigma_iter = 2;    //converged (1 vs 2: <0.1%)
      //msd rows: require the trapezoidal m0 to be converged, tested by
      //halving the beta resolution (error ratio 4:1), relative to the
      //deficit 1-m0 which sets the msd error scale. Catches coarse
      //kernels whose sharp small-alpha quasi-elastic peak is
      //under-sampled (bias invisible to the spread metric):
      constexpr double opt_msd_m0_convergence = 0.15;

      //Median and relative IQR spread of a scratch vector (modified):
      double calcMedian( VectD& v )
      {
        nc_assert( !v.empty() );
        std::sort( v.begin(), v.end() );
        const std::size_t n = v.size();
        return ( n % 2 == 1
                 ? vectAt( v, n/2 )
                 : 0.5 * ( vectAt( v, n/2 - 1 ) + vectAt( v, n/2 ) ) );
      }

      double calcRelIQR( VectD& v, double median )
      {
        //v must already be sorted; interpolation-free quartiles are
        //fine for a robustness metric:
        nc_assert( !v.empty() && median != 0.0 );
        const std::size_t n = v.size();
        const double q1 = vectAt( v, n/4 );
        const double q3 = vectAt( v, (3*n)/4 );
        return ( q3 - q1 ) / median;
      }

      //Zeroth beta-moment of a row using only every second retained
      //point, for the msd trapezoid-convergence self-test (an
      //under-resolved quasi-elastic peak converges badly):
      double calcRowM0Stride2( const VectD& beta, const VectD& sab,
                               std::size_t nalpha, std::size_t ia,
                               double floor )
      {
        StableSum m0;
        double bprev(0.0), wprev(0.0);
        unsigned nkept(0);
        bool have_prev(false), take(true);
        const std::size_t nb = beta.size();
        for ( std::size_t j = 0; j < nb; ++j ) {
          const double w = vectAt( sab, j*nalpha + ia );
          if ( !( w > floor ) )
            continue;
          if ( !take ) {
            take = true;
            continue;
          }
          take = false;
          ++nkept;
          const double b = vectAt( beta, j );
          if ( have_prev )
            m0.add( 0.5*( b - bprev )*( w + wprev ) );
          bprev = b; wprev = w; have_prev = true;
        }
        return nkept >= 2 ? m0.sum() : -1.0;
      }

      //Row moments over the floor-masked polyline (trapezoidal):
      struct RowMoments {
        double m0 = -1.0, mean = -1.0, var = -1.0;
        unsigned npts = 0;
        double b_first = 0.0, b_last = 0.0;
        unsigned npts_core = 0;//points within +-2 sigma of mean
      };

      RowMoments calcRowMoments( const VectD& beta, const VectD& sab,
                                 std::size_t nalpha, std::size_t ia,
                                 double floor, double sigma_guess )
      {
        RowMoments res;
        //First pass: m0 and m1 over the masked polyline:
        StableSum m0, m1;
        double bprev(0.0), wprev(0.0);
        bool have_prev = false;
        const std::size_t nb = beta.size();
        for ( std::size_t j = 0; j < nb; ++j ) {
          const double w = vectAt( sab, j*nalpha + ia );
          if ( !( w > floor ) )
            continue;
          const double b = vectAt( beta, j );
          if ( !res.npts ) {
            res.b_first = b;
          } else {
            const double halfdb = 0.5*( b - bprev );
            m0.add( halfdb*( w + wprev ) );
            m1.add( halfdb*( w*b + wprev*bprev ) );
          }
          res.b_last = b;
          bprev = b;
          wprev = w;
          have_prev = true;
          ++res.npts;
        }
        (void)have_prev;
        if ( res.npts < 10 || !( m0.sum() > 0.0 ) )
          return res;
        res.m0 = m0.sum();
        res.mean = m1.sum() / res.m0;
        //Second pass: central second moment + core-resolution count:
        StableSum m2;
        have_prev = false;
        for ( std::size_t j = 0; j < nb; ++j ) {
          const double w = vectAt( sab, j*nalpha + ia );
          if ( !( w > floor ) )
            continue;
          const double b = vectAt( beta, j );
          const double d = b - res.mean;
          if ( ncabs( d ) < 2.0*sigma_guess )
            ++res.npts_core;
          if ( have_prev ) {
            const double dprev = bprev - res.mean;
            m2.add( 0.5*( b - bprev )*( w*d*d + wprev*dprev*dprev ) );
          }
          bprev = b;
          wprev = w;
          have_prev = true;
        }
        res.var = m2.sum() / res.m0;
        return res;
      }

    }
  }
}

NC::SABAnalyser::Result
NC::SABAnalyser::analyse( const SABData& sab, const Options& opt )
{
#ifndef NDEBUG
  nc_assert_always( opt.center_tol > 0.0 && opt.center_tol < 1.0 );
#endif
  const auto& agrid = sab.alphaGrid();
  const auto& bgrid = sab.betaGrid();
  const auto& S = sab.sab();
  const std::size_t na = agrid.size();
  nc_assert_always( S.size() == na * bgrid.size() );
  const double T = sab.temperature().dbl();
  const double kT = sab.temperature().kT();
  const double massratio
    = sab.elementMassAMU().dbl() / const_neutron_mass_amu;

  Result res;

  //Global floor-marker level (placeholder values in far wings are not
  //data):
  double smin_pos = kInfinity, smax = 0.0;
  for ( auto v : S ) {
    if ( v > 0.0 )
      smin_pos = ncmin( smin_pos, v );
    smax = ncmax( smax, v );
  }
  if ( !( smax > 0.0 ) )
    return res;//empty/degenerate kernel
  const double floor = ncmax( smin_pos * opt_floor_relmin,
                              smax * opt_floor_relmax );

  std::shared_ptr<Diagnostics> diag;
  if ( opt.collect_diagnostics ) {
    diag = std::make_shared<Diagnostics>();
    diag->alpha = agrid;
    diag->m0.resize( na, -1.0 );
    diag->mean.resize( na, -1.0 );
    diag->variance.resize( na, -1.0 );
    diag->teff_row.resize( na, -1.0 );
    diag->msd_row.resize( na, -1.0 );
    diag->teff_row_status.resize( na, RowStatus::NoData );
    diag->msd_row_status.resize( na, RowStatus::NoData );
    diag->floor_value = floor;
  }

  //msd constant: 2W = alpha * alpha2x_factor * kT * msd (exactly as in
  //NCVDOSExpand.cc):
  constexpr double alpha2x_factor
    = 2.0*const_neutron_mass_evc2/(constant_hbar*constant_hbar);

  VectD teff_rows, msd_rows, center_rows;
  teff_rows.reserve( na );
  msd_rows.reserve( na );
  center_rows.reserve( na );

  //Teff needs a sigma guess for the coverage/resolution guards, refined
  //over n_sigma_iter passes (msd results are identical each pass, so
  //only the final pass records them):
  double r_guess = 2.0;
  for ( unsigned ipass = 0; ipass < opt_n_sigma_iter; ++ipass ) {
    const bool final_pass = ( ipass + 1 == opt_n_sigma_iter );
    teff_rows.clear();
    msd_rows.clear();
    center_rows.clear();
    for ( std::size_t ia = 0; ia < na; ++ia ) {
      const double a = vectAt( agrid, ia );
      const double sigma
        = std::sqrt( 2.0 * a * ncmax( r_guess, 1e-3 ) / massratio );
      auto rm = calcRowMoments( bgrid, S, na, ia, floor, sigma );
      auto setstat = [&diag,ia,final_pass]( RowStatus st, bool for_teff )
      {
        if ( diag != nullptr && final_pass ) {
          if ( for_teff )
            vectAt( diag->teff_row_status, ia ) = st;
          else
            vectAt( diag->msd_row_status, ia ) = st;
        }
      };
      if ( rm.m0 < 0.0 ) {
        setstat( RowStatus::NoData, true );
        setstat( RowStatus::NoData, false );
        continue;
      }
      if ( diag != nullptr && final_pass ) {
        vectAt( diag->m0, ia ) = rm.m0;
        vectAt( diag->mean, ia ) = rm.mean;
        vectAt( diag->variance, ia ) = rm.var;
      }
      //msd row (independent of sigma iteration):
      if ( final_pass ) {
        if ( !( a > opt_alpha_min ) ) {
          setstat( RowStatus::AlphaTooLow, false );
        } else if ( !( rm.m0 > opt_msd_m0_low
                       && rm.m0 < opt_msd_m0_high ) ) {
          setstat( RowStatus::M0OutOfWindow, false );
        } else if ( ncabs( calcRowM0Stride2( bgrid, S, na, ia, floor )
                           - rm.m0 )
                    > opt_msd_m0_convergence * ( 1.0 - rm.m0 ) ) {
          setstat( RowStatus::Unresolved, false );
        } else {
          const double c = -std::log1p( -rm.m0 ) / a;
          const double msd_row = c / ( alpha2x_factor * kT );
          msd_rows.push_back( msd_row );
          setstat( RowStatus::Ok, false );
          if ( diag != nullptr )
            vectAt( diag->msd_row, ia ) = msd_row;
        }
      }
      //Teff row:
      if ( !( a > opt_alpha_min ) ) {
        setstat( RowStatus::AlphaTooLow, true );
        continue;
      }
      if ( rm.b_first > rm.mean - opt_coverage_nsigma*sigma
           || rm.b_last < rm.mean + opt_coverage_nsigma*sigma ) {
        setstat( RowStatus::WingsNotCovered, true );
        continue;
      }
      if ( rm.npts_core < opt_resolution_minpts ) {
        setstat( RowStatus::Unresolved, true );
        continue;
      }
      const double c_row = rm.mean * massratio / ( -a );
      if ( !( ncabs( c_row - 1.0 ) <= opt.center_tol ) ) {
        setstat( RowStatus::CenterOffRecoil, true );
        continue;
      }
      const double teff_row = massratio * T * rm.var / ( 2.0 * a );
      teff_rows.push_back( teff_row );
      center_rows.push_back( c_row );
      setstat( RowStatus::Ok, true );
      if ( diag != nullptr && final_pass )
        vectAt( diag->teff_row, ia ) = teff_row;
    }
    if ( teff_rows.empty() )
      break;//no admissible rows; another pass will not change that
    VectD tmp( teff_rows );
    r_guess = calcMedian( tmp ) / T;
  }

  res.teff_nrows = static_cast<unsigned>( teff_rows.size() );
  res.msd_nrows = static_cast<unsigned>( msd_rows.size() );
  if ( !teff_rows.empty() ) {
    const double med = calcMedian( teff_rows );//NB: sorts
    res.teff = Temperature{ med };
    res.teff_relspread = calcRelIQR( teff_rows, med );
    res.recoil_center_ratio = calcMedian( center_rows );
  }
  if ( !msd_rows.empty() ) {
    const double med = calcMedian( msd_rows );//NB: sorts
    res.msd = med;
    res.msd_relspread = calcRelIQR( msd_rows, med );
  }
  res.diagnostics = std::move( diag );
  return res;
}

namespace NCRYSTAL_NAMESPACE {
  namespace SABAnalyser {
    namespace {

      //Acceptance policies for the turn-key estimates, per quantity
      //since the two consumers have opposite failure asymmetries:
      //
      //Teff: refusal is NOT neutral -- the consumer then extends with
      //free gas at T, an error of (Teff_true-T), usually far larger
      //than any credible estimate's error (e.g. H in 100K solid
      //ethanol: fallback ~12x off vs a refused 3% estimate). Hence a
      //loose policy: over 501 pooled truth-records (ENDF8- and
      //JENDL5-converted, stdlib and tests/data kernels incl.
      //temperature extremes), nrows>=10 && relspread<=0.1 never
      //accepted an estimate worse than the fallback, worst error 5.3%
      //(p95 well under 1%); garbage only appears below the row floor
      //(140% error at nrows=2).
      //
      //msd: refusal IS the status quo (no elastic component), while a
      //wrong value adds wrong physics, so strictness is correct. The
      //0.03 spread cut refuses the self-flagged tail, e.g. U in UO2
      //at 5K (spread 0.035, error 3%), for a worst accepted error of
      //0.77% (the msd values remain provisional pending ENDF-grade
      //ground truth):
      constexpr unsigned policy_teff_min_nrows = 10;
      constexpr double policy_teff_max_relspread = 0.1;
      constexpr unsigned policy_msd_min_nrows = 20;
      constexpr double policy_msd_max_relspread = 0.03;
      constexpr double policy_msd_sanity_max = 10.0;//[Aa^2]

      struct TeffMSDCacheDB {
        std::mutex mtx;
        //Both quantities are strictly positive, so -1.0 encodes
        //absence:
        struct Entry { double teff, msd; };
        std::map<std::uint64_t,Entry> map;
        bool cleanup_registered = false;
      };

      TeffMSDCacheDB& getTeffMSDCacheDB()
      {
        static TeffMSDCacheDB db;
        return db;
      }

      TeffMSD decodeEntry( const TeffMSDCacheDB::Entry& e )
      {
        TeffMSD res;
        if ( e.teff >= 0.0 )
          res.effectiveTemperature = Temperature{ e.teff };
        if ( e.msd >= 0.0 )
          res.msd = e.msd;
        return res;
      }

    }
  }
}

NC::SABAnalyser::TeffMSD
NC::SABAnalyser::estimateTeffMSD( const SABData& sab )
{
  auto& db = getTeffMSDCacheDB();
  const std::uint64_t uid = sab.getUniqueID().value;
  {
    NCRYSTAL_LOCK_GUARD(db.mtx);
    if ( !db.cleanup_registered ) {
      db.cleanup_registered = true;
      registerCacheCleanupFunction( []()
      {
        auto& thedb = getTeffMSDCacheDB();
        NCRYSTAL_LOCK_GUARD(thedb.mtx);
        thedb.map.clear();
      });
    }
    auto it = db.map.find( uid );
    if ( it != db.map.end() )
      return decodeEntry( it->second );
  }
  //Not cached; analyse outside the lock (a racing duplicate
  //computation is benign):
  auto r = analyse( sab );
  TeffMSD res;
  const double T = sab.temperature().dbl();
  if ( r.teff.has_value()
       && r.teff_nrows >= policy_teff_min_nrows
       && r.teff_relspread <= policy_teff_max_relspread ) {
    double tv = r.teff.value().dbl();
    //Same semantics as sct_checkedTeff, but discarding rather than
    //throwing:
    if ( std::isfinite( tv ) && tv >= 0.999*T ) {
      tv = ncmax( tv, T );
      if ( tv >= Temperature::allowed_range.first
           && tv <= Temperature::allowed_range.second )
        res.effectiveTemperature = Temperature{ tv };
    }
  }
  if ( r.msd.has_value()
       && r.msd_nrows >= policy_msd_min_nrows
       && r.msd_relspread <= policy_msd_max_relspread
       && std::isfinite( r.msd.value() )
       && r.msd.value() > 0.0
       && r.msd.value() < policy_msd_sanity_max )
    res.msd = r.msd.value();
  {
    NCRYSTAL_LOCK_GUARD(db.mtx);
    //Reset if somehow reaching an unreasonable size, on principle:
    if ( db.map.size() >= 5000 )
      db.map.clear();
    db.map[uid]
      = TeffMSDCacheDB::Entry{ ( res.effectiveTemperature.has_value()
                                 ? res.effectiveTemperature.value().dbl()
                                 : -1.0 ),
                               res.msd.has_value() ? res.msd.value()
                                                   : -1.0 };
  }
  return res;
}

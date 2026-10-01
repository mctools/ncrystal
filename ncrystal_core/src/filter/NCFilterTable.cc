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

#include "NCrystal/internal/filter/NCFilterTable.hh"
#include "NCrystal/factories/NCFactImpl.hh"
#include "NCrystal/interfaces/NCInfo.hh"
#include "NCrystal/interfaces/NCProcImpl.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCString.hh"

//FIXME: Get back to reproducibility for the filter tables: the selected table
//points depend on floating point details, which might differ between platforms
//and compilers (e.g. due to FMA contractions or auto-vectorisation).

namespace NC = NCrystal;
namespace NCF = NCrystal::Filter;

namespace NCRYSTAL_NAMESPACE {
  namespace Filter {
    namespace {

      //Relative distance from a discontinuity, at which the cross sections
      //just below and above it are evaluated. Discontinuities closer than
      //twice this are treated as a single one.
      constexpr double edge_eps = 1e-10;

      //The table is built with a slightly smaller tolerance than requested,
      //so the final verification at random points (which are not part of the
      //dense sample) has a margin:
      constexpr double tol_build_factor = 0.98;

      //Neighbouring dense sample points closer than this (relative) are not
      //refined further. If the tolerance is not met between such points, the
      //cross section has a discontinuity (or an extremely sharp feature) which
      //was not known in advance, and an exception is thrown:
      constexpr double min_rel_spacing = 1e-9;

      //The tolerance is relative to max(xs,abs_floor), with xs in barn per
      //atom (an absolute floor is needed for cross sections vanishing at
      //wavelength 0):
      constexpr double xs_abs_floor = 1e-12;

      //Max number of refinement iterations (each iteration at least halves
      //the spacing of the points where the tolerance is not met):
      constexpr unsigned max_refine_iterations = 100;

      //Leaf physics processes whose cross sections are known to have no
      //features which the table creation could miss: they are smooth, apart
      //from the Bragg edges (PowderBragg) and the boundaries of their energy
      //domains, which are both included explicitly. Any other process (e.g. a
      //process with absorption resonances, or the UCN processes) gives a
      //warning, and must be studied before being added here.
      bool isSupportedProcess( const std::string& name )
      {
        static const std::set<std::string> supported
          = { "NullScatter", "NullAbsorption", "FreeGas", "ElIncScatter",
              "SABScatter", "PowderBragg", "AbsOOV" };
        return supported.count( name ) > 0;
      }

      //Collect the leaf processes (of ProcComposition trees), and the
      //wavelengths corresponding to the boundaries of their energy domains:
      void inspectProcess( const ProcImpl::Process& proc,
                           std::set<std::string>& names, VectD& domain_wl )
      {
        auto comp = dynamic_cast<const ProcImpl::ProcComposition*>( &proc );
        if ( comp ) {
          for ( auto& c : comp->components() )
            inspectProcess( *c.process, names, domain_wl );
          return;
        }
        names.insert( proc.name() );
        const auto domain = proc.domain();
        if ( domain.isNull() )
          return;
        for ( auto e : { domain.elow.dbl(), domain.ehigh.dbl() } )
          if ( e > 0.0 && std::isfinite( e ) )
            domain_wl.push_back( ekin2wl( e ) );
      }

      void collectBraggEdges( const Info& info, VectD& edges )
      {
        if ( info.isMultiPhase() ) {
          for ( auto& ph : info.getPhases() )
            collectBraggEdges( *ph.second, edges );
          return;
        }
        if ( !info.hasHKLInfo() )
          return;
        for ( auto& hkl : info.hklList() )
          edges.push_back( 2.0 * hkl.dspacing );
      }

      //Sorted discontinuities within (0,wlmax), with clusters of nearby ones
      //represented by their largest member (for Bragg edges, the cross section
      //drops to the value above the whole cluster there).
      VectD mergeDiscontinuities( VectD d, double wlmax )
      {
        std::sort( d.begin(), d.end() );
        VectD res;
        res.reserve( d.size() );
        for ( auto v : d ) {
          if ( !( v > 0.0 && v < wlmax*(1.0-4*edge_eps) ) )
            continue;
          if ( !res.empty() && v <= res.back() * ( 1.0 + 2*edge_eps ) )
            res.back() = v;//merge with previous
          else
            res.push_back( v );
        }
        return res;
      }

      //Simple deterministic random numbers in (0,1), for the final
      //verification (splitmix64):
      class VerificationRNG {
        std::uint64_t m_state = 0x9e3779b97f4a7c15ull;
      public:
        double generate()
        {
          std::uint64_t z = ( m_state += 0x9e3779b97f4a7c15ull );
          z = ( z ^ ( z >> 30 ) ) * 0xbf58476d1ce4e5b9ull;
          z = ( z ^ ( z >> 27 ) ) * 0x94d049bb133111ebull;
          z = z ^ ( z >> 31 );
          return ( static_cast<double>( z >> 11 ) + 0.5 ) * ( 1.0 / 9007199254740992.0 );
        }
      };

      //Deviation of the chord from point i to point j, at point k:
      inline bool withinTol( const VectD& x, const VectD& y,
                             std::size_t i, std::size_t j, std::size_t k,
                             double tol, double abs_floor )
      {
        //NB: x[i]<x[j] is guaranteed by the callers when there are interior points.
        const double t = ( x[k] - x[i] ) / ( x[j] - x[i] );
        const double yc = y[i] + t * ( y[j] - y[i] );
        return ncabs( yc - y[k] ) <= tol * ncmax( ncabs( y[k] ), abs_floor );
      }

      //Check if chord i->j reproduces all interior points within tol:
      bool chordOK( const VectD& x, const VectD& y,
                    std::size_t i, std::size_t j, double tol, double abs_floor )
      {
        nc_assert( j > i );
        if ( j == i + 1 )
          return true;
        if ( !( x[j] > x[i] ) )
          return false;//vertical chord with interior points (never happens
                       //for valid input, where only pairs are identical).
        for ( std::size_t k = i + 1; k < j; ++k )
          if ( !withinTol( x, y, i, j, k, tol, abs_floor ) )
            return false;
        return true;
      }

      std::vector<std::size_t> reduceGreedy( const VectD& x, const VectD& y,
                                             double tol, double abs_floor )
      {
        //From the current point, find the most distant point which can be
        //reached with a chord satisfying the tolerance. An exponential search
        //followed by a binary search is used, which assumes that feasibility
        //is (approximately) monotonic in the end point. Every selected chord
        //is explicitly verified, so the tolerance is always respected at all
        //the sample points.
        const std::size_t n = x.size();
        std::vector<std::size_t> res;
        res.push_back( 0 );
        std::size_t i = 0;
        while ( i + 1 < n ) {
          std::size_t good = i + 1;
          std::size_t bad = n;//first known infeasible end point (n: none)
          std::size_t step = 2;
          while ( true ) {
            std::size_t j = std::min<std::size_t>( i + step, n - 1 );
            if ( j <= good )
              break;
            if ( chordOK( x, y, i, j, tol, abs_floor ) ) {
              good = j;
              if ( j == n - 1 )
                break;
              step *= 2;
            } else {
              bad = j;
              break;
            }
          }
          if ( bad < n ) {
            while ( bad - good > 1 ) {
              std::size_t mid = good + ( bad - good ) / 2;
              if ( chordOK( x, y, i, mid, tol, abs_floor ) )
                good = mid;
              else
                bad = mid;
            }
          }
          res.push_back( good );
          i = good;
        }
        return res;
      }

      std::vector<std::size_t> reduceDP( const VectD& x, const VectD& y,
                                         double tol, double abs_floor )
      {
        //Douglas-Peucker style: recursively split at the point with the
        //largest relative deviation from the chord.
        const std::size_t n = x.size();
        std::vector<char> keep( n, 0 );
        keep.front() = keep.back() = 1;
        std::vector<std::pair<std::size_t,std::size_t>> stack;
        stack.emplace_back( 0, n - 1 );
        while ( !stack.empty() ) {
          auto ij = stack.back();
          stack.pop_back();
          const std::size_t i = ij.first;
          const std::size_t j = ij.second;
          if ( j <= i + 1 )
            continue;
          std::size_t kworst = i + 1;
          double worst = -1.0;
          if ( !( x[j] > x[i] ) ) {
            worst = kInfinity;
          } else {
            for ( std::size_t k = i + 1; k < j; ++k ) {
              const double t = ( x[k] - x[i] ) / ( x[j] - x[i] );
              const double yc = y[i] + t * ( y[j] - y[i] );
              const double d = ncabs( yc - y[k] );
              const double yref = ncmax( ncabs( y[k] ), abs_floor );
              const double rel = ( yref > 0.0
                                   ? d / yref
                                   : ( d > 0.0 ? kInfinity : 0.0 ) );
              if ( rel > worst ) {
                worst = rel;
                kworst = k;
              }
            }
          }
          if ( worst > tol ) {
            keep[kworst] = 1;
            stack.emplace_back( i, kworst );
            stack.emplace_back( kworst, j );
          }
        }
        std::vector<std::size_t> res;
        for ( std::size_t k = 0; k < n; ++k )
          if ( keep[k] )
            res.push_back( k );
        return res;
      }

    }
  }
}

std::vector<std::size_t> NCF::reducePoints( const VectD& x, const VectD& y,
                                            double tol, ReductionAlgo algo,
                                            double abs_floor )
{
  nc_assert_always( x.size() == y.size() );
  nc_assert_always( x.size() >= 2 );
  return ( algo == ReductionAlgo::Greedy
           ? reduceGreedy( x, y, tol, abs_floor )
           : reduceDP( x, y, tol, abs_floor ) );
}

unsigned NCF::autoNDense( double tol, double wlmax )
{
  //10000 points per decade for tol=1e-3. The segments of the table get
  //shorter for smaller tolerances (roughly as sqrt(tol) in smooth regions), so
  //the dense sample must get denser:
  nc_assert_always( tol > 0.0 && wlmax > denseWavelengthMin() );
  const double ndecades = std::log10( wlmax / denseWavelengthMin() );
  const double n = ndecades * 10000.0 * std::sqrt( 1e-3 / tol );
  return static_cast<unsigned>( ncclamp( n, 1000.0, 1e7 ) );
}

double NCF::evalTable( const double* wl, const double* xs, std::size_t n,
                       double wavelength )
{
  if ( n == 0 )
    return 0.0;
  if ( n == 1 || !( wavelength > wl[0] ) )
    return xs[0];
  if ( wavelength >= wl[n-1] ) {
    //Linear extrapolation of the last segment, clamped at 0:
    const double dx = wl[n-1] - wl[n-2];
    if ( !( dx > 0.0 ) )
      return xs[n-1];
    const double slope = ( xs[n-1] - xs[n-2] ) / dx;
    if ( slope == 0.0 )
      return xs[n-1];
    const double v = xs[n-1] + ( wavelength - wl[n-1] ) * slope;
    return v > 0.0 ? v : 0.0;
  }
  //Last point with wl[i] <= wavelength (the second of a pair at a
  //discontinuity):
  const std::size_t i = static_cast<std::size_t>( std::upper_bound( wl, wl + n, wavelength ) - wl ) - 1;
  nc_assert( i + 1 < n && wl[i+1] > wl[i] );
  const double t = ( wavelength - wl[i] ) / ( wl[i+1] - wl[i] );
  return xs[i] + t * ( xs[i+1] - xs[i] );
}

NCF::Table NCF::createTable( const std::function<double(double)>& xs_of_wl,
                             VectD discontinuities_in,
                             const TableParams& params )
{
  if ( !( params.wlmax > 10.0 * denseWavelengthMin() && std::isfinite( params.wlmax ) ) )
    NCRYSTAL_THROW2( BadInput, "Filter table: invalid wlmax " << params.wlmax
                     << " (must be larger than " << 10.0 * denseWavelengthMin()
                     << " Aa)" );
  if ( !( params.tol > 0.0 && params.tol < 1.0 ) )
    NCRYSTAL_THROW2( BadInput, "Filter table: invalid tolerance "
                     << params.tol << " (must be in the range (0,1))" );
  if ( params.ndense == 1 )
    NCRYSTAL_THROW( BadInput, "Filter table: ndense must be at least 2" );
  const unsigned ndense = ( params.ndense
                            ? params.ndense
                            : autoNDense( params.tol, params.wlmax ) );
  const double tol = params.tol;
  const double tol_build = tol_build_factor * tol;

  auto evalXS = [&xs_of_wl]( double wl )
  {
    const double v = xs_of_wl( wl );
    if ( !( std::isfinite( v ) && v >= 0.0 ) )
      NCRYSTAL_THROW2( CalcError, "Filter table: invalid cross section " << v
                       << " at wavelength " << wl << " Aa" );
    return v;
  };

  const VectD discontinuities = ( params.discontinuities
                                  ? mergeDiscontinuities( std::move( discontinuities_in ),
                                                          params.wlmax )
                                  : VectD() );

  //Dense sample: wavelength 0, log-spaced points, and points just below and
  //above each discontinuity. The x-values are the wavelengths, where a
  //discontinuity contributes two identical x-values (with the cross sections
  //below and above it).
  VectD wl_eval;//wavelengths at which to evaluate
  VectD x;//the corresponding x values of the dense sample
  wl_eval.reserve( 1 + ndense + 2 * discontinuities.size() );
  x.reserve( 1 + ndense + 2 * discontinuities.size() );
  wl_eval.push_back( 0.0 );//NB: the limit is calculated below
  x.push_back( 0.0 );
  {
    const VectD grid = geomspace( denseWavelengthMin(), params.wlmax, ndense );
    auto itD = discontinuities.begin();
    for ( auto wl : grid ) {
      for ( ; itD != discontinuities.end() && *itD <= wl*(1.0+2*edge_eps); ++itD ) {
        const double e = *itD;
        if ( x.size() > 1 && x.back() >= e*(1.0-2*edge_eps) ) {
          //previous grid point too close to the discontinuity, remove it:
          x.pop_back();
          wl_eval.pop_back();
        }
        wl_eval.push_back( e * ( 1.0 - edge_eps ) );
        x.push_back( e );
        wl_eval.push_back( e * ( 1.0 + edge_eps ) );
        x.push_back( e );
      }
      if ( x.size() > 1 && wl <= x.back()*(1.0+2*edge_eps) )
        continue;//grid point too close to the previous discontinuity
      wl_eval.push_back( wl );
      x.push_back( wl );
    }
    nc_assert_always( itD == discontinuities.end() );
  }

  VectD y;
  y.reserve( wl_eval.size() );
  y.push_back( 0.0 );
  for ( std::size_t i = 1; i < wl_eval.size(); ++i )
    y.push_back( evalXS( wl_eval[i] ) );
  wl_eval.clear();
  const double abs_floor = xs_abs_floor;

  //The limit for wavelength -> 0, from linear extrapolation to 0 from very
  //short wavelengths (which is exact for any cross section of the form a+b*wl
  //there), at three scales. When the limit exists, the differences between
  //the estimates at consecutive scales decrease rapidly (e.g. by a factor of
  //100 for a term c*wl^2), while they do not for cross sections without a
  //limit (e.g. with 1/wl, log(wl) or sqrt(wl) behaviour):
  {
    auto extrapTo0 = [&evalXS]( double h ) { return 2.0 * evalXS( h ) - evalXS( 2.0 * h ); };
    const double lim_a = extrapTo0( 1e-7 );
    const double lim_b = extrapTo0( 1e-8 );
    const double lim_c = extrapTo0( 1e-9 );
    const double d_ab = ncabs( lim_a - lim_b );
    const double d_bc = ncabs( lim_b - lim_c );
    const double lim_tol = tol * ncmax( ncabs( lim_c ), abs_floor );
    //Either the estimates agree, or they converge (geometrically, with a
    //remaining error of at most d_bc*0.05/0.95 in lim_c):
    const bool converges = ( d_bc <= 1e-3 * lim_tol
                             || ( d_bc <= 0.05 * d_ab && 0.06 * d_bc <= lim_tol ) );
    if ( !converges )
      NCRYSTAL_THROW2( CalcError, "Filter table: the cross section does not"
                       " converge for wavelength -> 0 (estimated limits "
                       << lim_a << ", " << lim_b << " and " << lim_c << " barn)" );
    if ( lim_c < -abs_floor )
      NCRYSTAL_THROW2( CalcError, "Filter table: invalid limit " << lim_c
                       << " of the cross section for wavelength -> 0" );
    y.front() = ( lim_c > 0.0 ? lim_c : 0.0 );
  }

  //Reduce the dense sample, and verify the result at the midpoints between
  //all neighbouring sample points. Midpoints where the tolerance is not met
  //are added to the sample, and the procedure is repeated.
  std::vector<std::size_t> idx;
  std::size_t nrefine = 0;
  while ( true ) {
    idx = reducePoints( x, y, tol_build, params.algo, abs_floor );
    VectD add_x, add_y;
    std::size_t seg = 0;//table segment is from x[idx[seg]] to x[idx[seg+1]]
    for ( std::size_t k = 0; k + 1 < x.size(); ++k ) {
      while ( idx[seg+1] <= k )
        ++seg;
      if ( !( x[k+1] > x[k] ) )
        continue;//known discontinuity
      const std::size_t i = idx[seg];
      const std::size_t j = idx[seg+1];
      const double xm = 0.5 * ( x[k] + x[k+1] );
      const double ym = evalXS( xm );
      const double t = ( xm - x[i] ) / ( x[j] - x[i] );
      const double yc = y[i] + t * ( y[j] - y[i] );
      if ( ncabs( yc - ym ) <= tol_build * ncmax( ym, abs_floor ) )
        continue;
      if ( !( x[k+1] - x[k] > min_rel_spacing * x[k+1] ) )
        NCRYSTAL_THROW2( CalcError, "Filter table: the cross section has an"
                         " unexpected discontinuity (or an extremely sharp"
                         " feature) at wavelength " << xm << " Aa" );
      add_x.push_back( xm );
      add_y.push_back( ym );
    }
    if ( add_x.empty() )
      break;
    if ( ++nrefine > max_refine_iterations )
      NCRYSTAL_THROW( CalcError, "Filter table: the refinement of the table"
                      " did not converge" );
    //Merge the new points into the (sorted) sample:
    VectD new_x, new_y;
    new_x.reserve( x.size() + add_x.size() );
    new_y.reserve( x.size() + add_x.size() );
    std::size_t ia = 0;
    for ( std::size_t k = 0; k < x.size(); ++k ) {
      for ( ; ia < add_x.size() && add_x[ia] < x[k]; ++ia ) {
        new_x.push_back( add_x[ia] );
        new_y.push_back( add_y[ia] );
      }
      new_x.push_back( x[k] );
      new_y.push_back( y[k] );
    }
    nc_assert_always( ia == add_x.size() );
    x.swap( new_x );
    y.swap( new_y );
  }

  Table res;
  res.wl.reserve( idx.size() );
  res.xs.reserve( idx.size() );
  for ( auto i : idx ) {
    res.wl.push_back( x[i] );
    res.xs.push_back( y[i] );
  }
  nc_assert_always( res.wl.size() >= 2 && res.wl.front() == 0.0
                    && res.wl.back() == params.wlmax
                    && res.wl[res.wl.size()-2] < params.wlmax );

  //Final verification at random points in each segment (in addition to the
  //sample points, which are verified by construction):
  {
    VerificationRNG rng;
    constexpr unsigned npts_per_segment = 4;
    for ( std::size_t i = 0; i + 1 < res.wl.size(); ++i ) {
      const double x0 = res.wl[i];
      const double x1 = res.wl[i+1];
      if ( !( x1 > x0 ) )
        continue;
      for ( unsigned k = 0; k < npts_per_segment; ++k ) {
        const double xv = x0 + ( 0.01 + 0.98 * rng.generate() ) * ( x1 - x0 );
        const double yv = evalXS( xv );
        const double yt = evalTable( res.wl.data(), res.xs.data(), res.wl.size(), xv );
        if ( !( ncabs( yt - yv ) <= tol * ncmax( yv, abs_floor ) ) )
          NCRYSTAL_THROW2( CalcError, "Filter table: verification failed at"
                           " wavelength " << xv << " Aa (table: " << yt
                           << " barn, exact: " << yv << " barn)" );
      }
    }
  }

  res.numberDensity = 0.0;
  res.nDiscontinuities = discontinuities.size();
  res.nDense = x.size();
  res.nRefine = nrefine;
  return res;
}

NCF::Table NCF::createTable( const MatCfg& cfg, const TableParams& params )
{
  auto info = FactImpl::createInfo( cfg );
  auto scat = FactImpl::createScatter( cfg );
  auto absn = FactImpl::createAbsorption( cfg );
  if ( scat->isOriented() || absn->isOriented() )
    NCRYSTAL_THROW( BadInput, "Filter table: only isotropic materials are"
                    " supported (the cfg-string must not specify a crystal"
                    " orientation)" );

  std::set<std::string> names;
  VectD discontinuities;
  inspectProcess( *scat, names, discontinuities );
  inspectProcess( *absn, names, discontinuities );
  for ( auto& n : names )
    if ( !isSupportedProcess( n ) )
      NCRYSTAL_WARN( "Filter table: the physics process \"" << n << "\" has"
                     " not been validated for filter tables (features of its"
                     " cross section might be missed)" );
  collectBraggEdges( *info, discontinuities );

  CachePtr cp_scat, cp_abs;
  auto xs_of_wl = [&scat,&absn,&cp_scat,&cp_abs]( double wl )
  {
    const NeutronEnergy ekin{ wl2ekin( wl ) };
    return ( scat->crossSectionIsotropic( cp_scat, ekin ).dbl()
             + absn->crossSectionIsotropic( cp_abs, ekin ).dbl() );
  };

  Table res = createTable( xs_of_wl, std::move( discontinuities ), params );
  res.numberDensity = info->getNumberDensity().dbl();
  res.processes.assign( names.begin(), names.end() );
  return res;
}

void NCF::JSONQuery( std::ostream& os, const Query& query )
{
  auto invalid = [&query](const std::string& reason){
    std::ostringstream ss;
    ss << "Invalid filtertable query: ";
    streamJSON( ss, query );
    ss << " (" << reason << ")";
    NCRYSTAL_THROW( BadInput, ss.str() );
  };
  if ( query.size() < 2 || query.front() != "filtertable" )
    invalid("usage: [\"filtertable\",CFGSTR,OPTIONS...] with options like"
            " \"tol=1e-3\", \"wlmax=500\", \"ndense=0\" (auto),"
            " \"discontinuities=1\" and \"algo=greedy\" (or \"algo=dp\")");
  const std::string cfgstr = query.at(1).to_string();

  TableParams params;
  for ( std::size_t i = 2; i < query.size(); ++i ) {
    auto opt = query.at(i).trimmed();
    auto parts = opt.split<2>('=');
    if ( parts.size() != 2 )
      invalid("options must have the form KEY=VALUE");
    auto key = parts.at(0).trimmed();
    auto val = parts.at(1).trimmed();
    auto getdbl = [&]() {
      auto v = val.toDbl();
      if ( !v.has_value() )
        invalid("invalid value for option "+key.to_string());
      return v.value();
    };
    auto getuint = [&]() {
      auto v = val.toInt();
      if ( !v.has_value() || v.value() < 0 || v.value() > 1000000000 )
        invalid("invalid value for option "+key.to_string());
      return static_cast<unsigned>( v.value() );
    };
    if ( key == "tol" ) {
      params.tol = getdbl();
    } else if ( key == "wlmax" ) {
      params.wlmax = getdbl();
    } else if ( key == "ndense" ) {
      params.ndense = getuint();
    } else if ( key == "discontinuities" ) {
      params.discontinuities = ( getuint() != 0 );
    } else if ( key == "algo" ) {
      if ( val == "greedy" )
        params.algo = ReductionAlgo::Greedy;
      else if ( val == "dp" )
        params.algo = ReductionAlgo::DouglasPeucker;
      else
        invalid("algo must be \"greedy\" or \"dp\"");
    } else {
      invalid("unknown option "+key.to_string());
    }
  }

  auto table = createTable( MatCfg( cfgstr ), params );

  streamJSONDictEntry( os, "cfgstr", cfgstr, JSONDictPos::FIRST );
  streamJSONDictEntry( os, "tol", params.tol );
  streamJSONDictEntry( os, "wlmax", params.wlmax );
  streamJSONDictEntry( os, "ndense", ( params.ndense
                                       ? params.ndense
                                       : autoNDense( params.tol, params.wlmax ) ) );
  streamJSONDictEntry( os, "discontinuities", params.discontinuities );
  streamJSONDictEntry( os, "algo", ( params.algo == ReductionAlgo::Greedy
                                     ? "greedy" : "dp" ) );
  streamJSONDictEntry( os, "processes", table.processes );
  streamJSONDictEntry( os, "numberdensity", table.numberDensity );
  streamJSONDictEntry( os, "ndiscontinuities", table.nDiscontinuities );
  streamJSONDictEntry( os, "ndense_total", table.nDense );
  streamJSONDictEntry( os, "nrefine", table.nRefine );
  streamJSONDictEntry( os, "npts", table.wl.size() );
  streamJSONDictEntry( os, "wl", table.wl );
  streamJSONDictEntry( os, "xs", table.xs, JSONDictPos::LAST );
}

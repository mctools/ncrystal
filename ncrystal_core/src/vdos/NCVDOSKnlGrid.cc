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

#include "NCrystal/internal/vdos/NCVDOSKnlGrid.hh"
#include "NCrystal/internal/vdos/NCVDOSUtils.hh"
#include "NCrystal/internal/utils/NCSpan.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCString.hh"

//TODO: reduce usage of nc_assert_always in this file once the model has been
//      used in production for a while.

namespace NC=NCrystal;

namespace NCRYSTAL_NAMESPACE {

  namespace VDOS {

    namespace {

      // evalPWLSum: Evaluates the weighted sum of piecewise-linear functions on
      // the supplied grid.
      //
      // Each function is zero outside its own [x0, xmax()] interval.
      // Values inside an interval are obtained by linear interpolation.
      // `fs` contains the input functions and `ws` their corresponding weights.
      // The returned vector has one value for each point in `grid`.

      struct PWLFct {
        double x0; //x{i=0}
        double binWidth;//distance between x{i} and x{i+1}
        Span<const double> f;//values of f at the x{i} points. The size of the
                             //span encodes the number of points.
        double x1() const { return x0 + (f.size()-1)*binWidth; }
        VectD dataHolder;//optional, so can hold its data if needed.
        Span<double> f_mutable()
        {
          //only possible when we hold our data
          nc_assert( dataHolder.size() == f.size() );
          return dataHolder;
        }
      };

      PWLFct pwlNarrowToPos(const PWLFct& p, double tol = 1e-3 )
      {
        nc_assert( p.binWidth > 0.0 );
        const double t = tol * p.binWidth;//threshold
        std::size_t i = ( p.x0 < t ? static_cast<std::size_t>
                          (std::ceil((t - p.x0) / p.binWidth)) : 0u );
        i = std::min<std::size_t>(i, p.f.size());
        PWLFct res;
        res.x0 = p.x0 + i * p.binWidth;
        res.binWidth = p.binWidth;
        res.dataHolder.assign(p.f.begin() + i, p.f.end());
        res.f = res.dataHolder;
        nc_assert_always( res.x0 >= t );
        nc_assert_always( res.f.size() >= 2 );
        return res;
      }

      VectD evalPWLSum(Span<const PWLFct> fs,
                       Span<const double> grid,
                       Span<const double> ws = {})
      {
        VectD out(grid.size(), 0.0);

        const double* ncrestrict gridPtr = grid.data();
        double* ncrestrict outPtr = out.data();

        const bool weighted = !ws.empty();

        for (std::size_t n = 0; n < fs.size(); ++n) {
          const PWLFct& p = fs[n];
          const auto& f = p.f;

          const std::size_t lastBin = f.size() - 1;
          const double invbw = 1.0 / p.binWidth;
          const double weight = weighted ? ws[n] : 1.0;
          const double xmax = p.x1();

          std::size_t g = static_cast<std::size_t>
            (std::lower_bound(grid.begin(), grid.end(), p.x0) - grid.begin());

          if (g == grid.size() || gridPtr[g] > xmax)
            continue;

          for (std::size_t i = 0; i < lastBin && g < grid.size(); ++i) {
            const double xLeft =
              p.x0 + static_cast<double>(i) * p.binWidth;
            const double xRight = xLeft + p.binWidth;
            const double y0 = vectAt(f, i);
            const double y1 = vectAt(f, i + 1);
            const double slope = (y1 - y0) * invbw;
            const double intercept = y0 - slope * xLeft;
            std::size_t end = g;
            while (end < grid.size() && gridPtr[end] < xRight)
              ++end;

            // This loop is deliberately simple so the compiler can
            // auto-vectorize it.
            for (std::size_t k = g; k < end; ++k)
              outPtr[k] += weight * (intercept + slope * gridPtr[k]);
            g = end;
          }

          // Handle the final sample at x == xmax.
          if (g < grid.size() && gridPtr[g] <= xmax) {
            const double y = vectAt(f, lastBin);
            std::size_t end = g;
            while (end < grid.size() && gridPtr[end] <= xmax)
              ++end;
            for (std::size_t k = g; k < end; ++k)
              outPtr[k] += weight * y;
          }
        }

        return out;
      }

      // Merge two finite, sorted, unique grids.
      // Endpoint values always take priority over averaging.
      VectD mergeGridsWithTol( const VectD& a, const VectD& b,
                               double rtol = 0.01 )
      {
        nc_assert( rtol > 0.0 );
        nc_assert( rtol < 0.5 );
        nc_assert(nc_is_grid(a));
        nc_assert(nc_is_grid(b));
        nc_assert(!a.empty());
        nc_assert(!b.empty());
        nc_assert(std::isfinite(rtol) && rtol > 0.0);
        const double fac = 1.0 + rtol;
        const double lo = ncmin(a.front(), b.front());
        const double hi = ncmax(a.back(), b.back());
        VectD g;
        g.reserve(a.size() + b.size());
        auto add = [&g](double x)
        {
          if (g.empty() || x != g.back())
            g.push_back(x);
        };
        auto farEnough = [fac](double x, double y)
        {
          if (x == y || x == 0.0 || y == 0.0)
            return x != y;
          if ((x < 0.0) != (y < 0.0))
            return true;
          const double ax = ncabs(x);
          const double ay = ncabs(y);
          return ncmax(ax, ay) / ncmin(ax, ay) > fac;
        };
        std::size_t ia = 0;
        for (std::size_t ib = 0; ib < b.size(); ++ib) {
          const double x = vectAt(b, ib);
          while (ia < a.size() && vectAt(a, ia) < x) {
            add(vectAt(a, ia));
            ++ia;
          }
          bool keep = x == lo || x == hi;
          if (!keep && ia < a.size())
            keep = farEnough(x, vectAt(a, ia));
          if (!keep && ia > 0)
            keep = farEnough(x, vectAt(a, ia - 1));
          if (keep)
            add(x);
        }
        while (ia < a.size()) {
          add(vectAt(a, ia));
          ++ia;
        }
        return g;
      }

      double combinedGnFctWeight( VDOSGn::Order n )
      {
        // We want to assign higher weight to lower n values, and highest of all to
        // n=1. w(n)=1/n is unstable since sum of all weights does not converge, so
        // we go for ~1/n^2 instead. We also choose a form where the sum to n=inf is
        // normalised, and the parameter we tune is the weight of n=1. It is a bit
        // of trial and error, and likely not super crucial how this is setup. But
        // the following seemed to give sensible results.
        constexpr double p1 = 0.3;//1st order phonon weight
        constexpr double a = 1.0/p1-2.0;
        constexpr double ap1 = a + 1.0;
        return ap1 / ((n.value()+a)*(n.value()+ap1));
      }

      // Calculates the trapezoidal integral of a non-negative piecewise-linear
      // function represented by x and y.
      // Removes as many points from the back as possible while keeping the
      // discarded integral at most frac times the original integral.
      // The vectors are modified in place; the original integral is returned.
      void trimTailByIntegral( VectD& x, VectD& y, double frac )
      {
#ifndef NDEBUG
        const auto npts = x.size();
        nc_assert(npts == y.size());
        nc_assert(npts >= 2);
        nc_assert(std::isfinite(frac) && frac >= 0.0);
        nc_assert( nc_is_grid(x) );
        for (auto& g : y) {
          nc_assert(std::isfinite(g));
          nc_assert(g >= 0.0);
        }
#endif
        double integral;
        {
          StableSumKahan sum;
          for (std::size_t i = 1; i < x.size(); ++i) {
            const double dx = vectAt(x, i) - vectAt(x, i - 1);
            sum.add( dx * (vectAt(y, i - 1) + vectAt(y, i)) );
          }
          integral = sum.sum()*0.5;
        }

        const double limit = frac * integral;
        double removed = 0.0;

        while (x.size() > 2) {
          const auto i = x.size() - 2;
          const double dx = vectAt(x, i + 1) - vectAt(x, i);
          const double area = 0.5 * dx * (vectAt(y, i) + vectAt(y, i + 1));
          if (removed + area > limit)
            break;
          removed += area;
          x.pop_back();
          y.pop_back();
        }
      }

      void trimEquidistantGridUpperEdge(EquidistantGrid& g, double xmax)
      {
#ifndef NDEBUG
        nc_assert(g.npts >= 2);
        nc_assert(std::isfinite(g.x0));
        nc_assert(std::isfinite(g.binWidth));
        nc_assert(g.binWidth > 0.0);
        nc_assert(std::isfinite(g.x1()));
        nc_assert(std::isfinite(xmax));
        EquidistantGrid expected = g;
        auto check = [&expected,&g] {
          return ( expected.x0 == g.x0
                   && expected.binWidth == g.binWidth
                   && expected.npts == g.npts );
        };
        nc_assert(check());
        while ( expected.npts > 2 && expected.x1()>xmax )
          --expected.npts;
#endif
        if ( g.npts <= 2 || g.x1() <= xmax ) {
          nc_assert(check());
          return;
        }

        if ( xmax <= g.x0 + g.binWidth ) {
          g.npts = 2;
          nc_assert(check());
          return;
        }

        const std::size_t orig_npts = g.npts;
        const double q = (xmax - g.x0) / g.binWidth;
        nc_assert( std::isfinite(q) && q >= 0.0 );
        if (q < 1.0) {
          g.npts = 2;
        } else if (q >= static_cast<double>(orig_npts - 1)) {
          nc_assert( g.npts == orig_npts );
        } else {
          const double k = std::floor(q);
          nc_assert(std::isfinite(k) && k >= 1.0);
          nc_assert(k < static_cast<double>(orig_npts - 1));
          const std::size_t npts = static_cast<std::size_t>(k) + 1;
          g.npts = orig_npts < npts ? orig_npts : npts;
        }

        // Correct rounding errors:
        while (g.npts > 2 && g.x1() > xmax)
          --g.npts;
        while (g.npts < orig_npts) {
          ++g.npts;
          if (g.x1() > xmax) {
            --g.npts;
            break;
          }
        }
        nc_assert(check());
      }

    }
  }
}

std::pair<NC::VectD,NC::VectD>
NC::VDOS::getCombinedGnFct( const GnExpansion& gnexpn )
{
  const auto& Gn = gnexpn.Gn;
  const auto nmax = Gn.maxOrder().value();
  const double invkT = 1.0 / Gn.kT();
  std::vector<PWLFct> fs;
  VectD ws;
  fs.reserve(nmax);
  ws.reserve(nmax);
  for ( auto nm1 : ncrange(nmax) ) {
    const auto n = nm1+1;
    fs.emplace_back();
    fs.back().x0 = Gn.eRange(n).first*invkT;
    fs.back().binWidth = Gn.binWidth(n)*invkT;
    fs.back().f = Gn.getRawSpectrum(n);
    ws.push_back(combinedGnFctWeight( n ));
  };
  std::vector<EquidistantGrid> individual_grids;
  individual_grids.reserve(nmax);
  {
    unsigned n=0;
    for ( auto& f : fs ) {
      ++n;
      double bw = f.binWidth;
      std::size_t npts = f.f.size();
      //Speed up (a lot!) by already at this stage thinning out point candidates
      //at higher n values.
      while ( npts > 1000 && n > 2 ) {
        bw *= 2;
        npts = (npts+1)/ 2;
      }
      while ( npts > 200 && n > 4 ) {
        bw *= 2;
        npts = (npts+1)/ 2;
      }
      while ( npts > 125 && n > 12 ) {
        bw *= 2;
        npts = (npts+1)/ 2;
      }
      while ( npts > 80 && n > 100 ) {
        bw *= 2;
        npts = (npts+1)/ 2;
      }
      if ( n >= 8
           && !individual_grids.empty()
           && individual_grids.back().binWidth == bw
           && intervalsOverlap( individual_grids.back().x0,
                                individual_grids.back().x1(),
                                f.x0, f.x1() ) ) {
        //Speed up by simply adding on to existing grid. Yes, this is not
        //strictly the same, since the thinning+trimming might have left the two
        //grids "out of phase", so it is not true that the two x0's are offset
        //by an integral number of bin widths. But at high n, we are anyway
        //looking at very smooth bell-curves whose width is much larger than the
        //binwidths.
        individual_grids.back()
          = coverEquidistantGrids( individual_grids.back(),
                                   EquidistantGrid{ f.x0, bw, npts } );

      } else {
        individual_grids.emplace_back();
        individual_grids.back().x0 = f.x0;
        individual_grids.back().binWidth = bw;
        individual_grids.back().npts = npts;
      }
      trimEquidistantGridUpperEdge(individual_grids.back(),
                                   gnexpn.sabRange.y1()+bw);

    }
  }

  VectD grid = makeCommonGrid( individual_grids );
  VectD vals = evalPWLSum(fs,grid,ws);
  std::pair<VectD,VectD> res;
  res.first = std::move(grid);
  res.second = std::move(vals);
  return res;
}

std::pair<NC::VectD,NC::VectD>
NC::VDOS::determineAlphaBetaGridFromGn( const GnExpansion& gnexpn,
                                        unsigned nalpha, unsigned nbeta )
{
  //fixme: parallelise parts with FactoryJobs?

  //fixme: Tune powspace powers for the alpha ranges (middle range might
  //        even be linspace)..

  const unsigned requested_nalpha = nalpha;
  const unsigned requested_nbeta = nbeta;

  //For high quality performance in the low E limit, we dedicate some points
  //just for that purpose:
  VectD e0grid;
  {
    unsigned npts_E0
      = std::max<unsigned>(12u,static_cast<unsigned>
                           (std::min<unsigned>(nalpha,nbeta)*0.2+0.5));
    VectD tmp;
    std::tie(e0grid,tmp) = setupE0ABGrid( gnexpn, npts_E0 );
    nc_assert_always(e0grid.size()>=5);
    nc_assert_always(nbeta > e0grid.size()+25);
    nc_assert_always(nalpha > e0grid.size()+20);
    nbeta -= static_cast<unsigned>(e0grid.size());
    nalpha -= static_cast<unsigned>(e0grid.size());
    nc_assert_always( nbeta > 30 && nbeta < 1000000 );
    nc_assert_always( nalpha > 20 && nalpha < 1000000 );
  }


  const double alpha2x = gnexpn.alpha2x;

  const PairDD betaRange = gnexpn.sabRange.yRange();
  const double alphaMax = gnexpn.sabRange.x1();
  nc_assert_always(alpha2x > 0.0);
  nc_assert_always(alphaMax > 0.0);

  VectD bvals, gnprojvals;
  std::tie(bvals, gnprojvals) = getCombinedGnFct(gnexpn);
  nc_assert_always( bvals.size() >= 5 );
  nc_assert_always( bvals.front() < 0.0 );
  nc_assert_always( bvals.back() > 0.0 );

  Span<const double> bvals_view(bvals);
  Span<const double> gnprojvals_view(gnprojvals);
  {
    //Discard pts outside betaRange, but occasionally keep one point extra, to
    //avoid an edge-effects due to an inadvertent extrapolation towards 0 in the
    //edge region.
    std::size_t i = 0;
    while ( bvals_view[i] < betaRange.first )
      ++i;
    if ( i > 0 && bvals_view[i] > betaRange.first )
      --i;//keep one point going over the edge.

    bvals_view = bvals_view.subspan(i);
    gnprojvals_view = gnprojvals_view.subspan(i);
    if ( bvals_view.back() > betaRange.second ) {
      //Do the same for the upper limit, although this is expected to happen
      //only extremely rarely in usual operations.
      auto newsize = bvals_view.size();
      nc_assert_always( bvals_view.size() >= 5 );
      while ( newsize > 2
              && ncmin(bvals_view[newsize-1],
                       bvals_view[newsize-2]) > betaRange.second ) {
        --newsize;
      }
      bvals_view = bvals_view.subspan(0,newsize);
      gnprojvals_view = gnprojvals_view.subspan(0,newsize);
    }

    nc_assert_always( bvals_view.size() >= 5 );
    nc_assert_always( bvals_view.front() < 0.0 );
    nc_assert_always( bvals_view.back() > 0.0 );
    nc_assert_always( bvals_view.size() == gnprojvals_view.size() );
  }

  //Reduce number of points:
  std::tie(bvals, gnprojvals)
    = reducePtsInDistribution( bvals_view, gnprojvals_view, nbeta );
  nc_assert_always( bvals.size() >= 2 );
  nc_assert_always( bvals.size() <= nbeta );
  nc_assert_always( bvals.size() == gnprojvals.size() );
  nc_assert_always( bvals.front() < 0.0 );
  nc_assert_always( bvals.back() > 0.0 );

  //Alpha is easier, since we know the behaviour of the formulas
  //f(x,n)=exp(-x)*x^n/n!. Investigations show that a powspace with p=2 gives
  //consistent relative errors over the entire range. This is the same as
  //linearly spacing in sqrt(x) values, and since the f(x,n) values approach
  //gaussians with sigma=sqrt(n) at high x values, this is perhaps not
  //surprising in hindsight, at least not at high x values.

  const double xmax = alpha2x*alphaMax;
  nc_assert_always(alpha2x>0.0&&std::isfinite(alpha2x));
  const double x2alpha = 1.0/alpha2x;
  nc_assert_always(x2alpha>0.0&&std::isfinite(x2alpha));

  constexpr double alphapow = 2.0;
  constexpr double alphapow_low = 4.0;
  unsigned nlow = static_cast<unsigned>(nalpha*0.08+4.5);
  unsigned nmid = static_cast<unsigned>(nalpha*0.12+10.5);
  const unsigned nzero = 1;//always begin with 0
  while ( !(nalpha > nlow+nmid+nzero+5) ) {
    nlow -= 1;
    nmid -= 1;
  }
  nc_assert_always(nlow >= 3);
  nc_assert_always(nmid >= 8);
  nc_assert_always(nalpha > nlow+nmid+nzero+5);

  const auto npow = nalpha-(nlow+nmid+nzero);
  const double xpow0 = ncmin(2.5,0.2*xmax);
  const double xmid0 = ncmin(0.15,0.2*xpow0);

  //We start xlow at 0.005 since the fall-back to linear interpolation over
  //[0.0,xlow0] is exactly the ~linear behaviour exhibited by the n=1 curve at
  //low x! However, note that the merging of E0 pts later might actually ruin
  //this. There is not really anything to be done about that.
  const double xlow0 = ncmin( 0.005, 0.1*xmid0 );

  auto xlow = powspace(xlow0,xmid0,nlow+1,alphapow_low);
  xlow.pop_back();
  auto xmid = powspace(xmid0,xpow0,nmid+1,1.5);
  xmid.pop_back();
  auto xpow = powspace( xpow0, xmax, npow, alphapow );

  VectD avals;
  avals.reserve( requested_nalpha );
  avals.push_back(0.0);
  for ( auto& x : xlow )
    avals.push_back( x*x2alpha );

  for ( auto& x : xmid )
    avals.push_back( x*x2alpha );
  for ( auto& x : xpow )
    avals.push_back( x*x2alpha );
  nc_assert_always(floateq(avals.back(),alphaMax));
  avals.back() = alphaMax;//force exact
  nc_assert_always(nc_is_grid(avals));
  nc_assert_always(avals.size()==nalpha);

  //Before returning we need to merge the e0grid points into both alpha and beta
  //grids, ensure we avoid pts too close, and to ensure the exact final
  //count. For alpha, it is possible for extreme values of alpha2x (e.g. in
  //low-temp uranium) to get a gap between the last e0 pt and the first xlow
  //point. If this happens, we try to fill it out.

  nc_assert( !e0grid.empty() );
  if ( e0grid.back() < xlow.front()*x2alpha ) {
    auto nn = static_cast<std::size_t>(e0grid.size() * 0.1+0.5);
    e0grid.resize( e0grid.size()-nn );
    auto tmp = linspace( e0grid.back(), xlow.front()*x2alpha, nn + 2 );
    for ( auto e : Span<const double>(tmp).subspan(1,nn) )
      e0grid.push_back( e );
  }
  {
    const double rtol = 0.01;//avoid pts too close to each other
    avals = mergeGridsWithTol( avals, e0grid, rtol );
    auto topOff = [rtol](VectD& v, std::size_t n)
    {
      double rtol_try = 10*rtol;
      while ( v.size() < n ) {
        topOffGrid(v, n, rtol_try );
        rtol_try *= 0.25;
      }
    };

    topOff(avals, requested_nalpha );

    const std::size_t bvals_size_withoute0 = bvals.size();
    bvals = mergeGridsWithTol( bvals, e0grid, rtol );
    nc_assert_always( bvals.size() >= bvals_size_withoute0 );
    nc_assert( nc_is_grid( bvals ) );
    topOff(bvals, requested_nbeta );
  }

  //Finally ensure that our grid values are sufficiently spaced out to not give
  //numerical issues when using the grid (e.g. when doing numerical integrations
  //in relative coords within a grid cell, or when exporting to ENDF where we
  //can't use full double precision for grid coordinates).
  constexpr double abGridSpaceOutTol = 0.0001;
  spaceOutGrid(avals,abGridSpaceOutTol);
  spaceOutGrid(bvals,abGridSpaceOutTol);

  //Final sanity checks on the guarantees:
  // -> both grids always contain 0.0:
  nc_assert( std::binary_search(bvals.begin(), bvals.end(), 0.0) );
  nc_assert_always( avals.front() == 0.0 );
  // -> both have the exact dimension requested:
  nc_assert_always( avals.size() == requested_nalpha );
  nc_assert_always( bvals.size() == requested_nbeta );
  // -> both are grids of course:
  nc_assert( nc_is_grid( avals ) );
  nc_assert( nc_is_grid( bvals ) );
  return {avals, bvals};
}

std::pair<unsigned,unsigned> NC::VDOS::gridDimFromLux( VDOSLux vdoslux )
{
  const unsigned override_nbins = ncgetenv_int("HACK_NBINS");
  if ( override_nbins )
    return { override_nbins, override_nbins };

  if ( vdoslux.isLegacy() ) {
    nc_assert_always( vdoslux.lvl() <= 5 );
    return { 50*(1<<vdoslux.lvl()), 100*(1<<vdoslux.lvl()) };
  }

  nc_assert_always( vdoslux.lvl() <= 6 );

  switch ( vdoslux.lvl() ) {
  case 0u: return {   60,  130 };//was 50x100 in legacy scheme
  case 1u: return {  110,  240 };//was 100x200 in legacy scheme
  case 2u: return {  170,  400 };//was 200x400 in legacy scheme
  default:
  case 3u: return {  280,  650 };//was 400x800 in legacy scheme
  case 4u: return {  470, 1080 };//was 800x1600 in legacy scheme
  case 5u: return {  770, 1770 };//was 1600x3200 in legacy scheme
  case 6u: return { 2000, 4000 };//not present in legacy scheme
  }
}

std::pair<NC::VectD,NC::VectD>
NC::VDOS::setupE0ABGrid( const GnExpansion& gnexpn, unsigned npts )
{
  nc_assert_always(npts >= 5);
  const VDOSGn& Gn = gnexpn.Gn;
  const double alpha2x = gnexpn.alpha2x;
  nc_assert_always(alpha2x>0.0);
  const double x2alpha = 1.0 / alpha2x;
  nc_assert_always(x2alpha>0.0);
  const auto nmax = Gn.maxOrder().value();
  const double kT = Gn.kT();
  const double invkT = 1.0/kT;
  std::vector<PWLFct> fcts;
  fcts.reserve(16);

  double minus_log_nfactorial = 0.0;// accumulate -ln(n!)
  for ( auto nm1 : ncrange(nmax) ) {
    const auto n = nm1+1;
    minus_log_nfactorial -= std::log(static_cast<double>(n));
    if ( n>1 ) {
      //check if we can break already
      const double betamax = Gn.eRange(n).second*invkT;
      if ( betamax <= 0.0 )
        break;
      constexpr double relcontriblvl = 1e-8;
      constexpr double accuracy = 1e-11;
      auto xrange = rangeXNexpMX(n, relcontriblvl, accuracy );
      if ( betamax < xrange.first * x2alpha )
        break;
    }

    //Collect contributions on the line alpha=beta, but modified with a factor
    //of sqrt(beta) to account for the relative width of the phase space at
    //different beta values.

    //shave off negative beta:
    fcts.emplace_back();
    auto& f = fcts.back();
    f.x0 = Gn.eRange(n).first * invkT;
    f.binWidth = Gn.binWidth(n)*invkT;
    nc_assert_always(f.binWidth>0.0);
    f.f = Gn.getRawSpectrum(n);
    //Discard non-positive values:
    f = pwlNarrowToPos(f);
    nc_assert_always( f.x0 > 0.0 && f.f.size() >= 2 );
    auto f_fmut = f.f_mutable();

    //Now add alpha and phase space factors:
    for (std::size_t i = 0; i < f.f.size(); ++i) {
      //note alpha=beta on the E->0 phasespace, so x=beta*alpha2x
      const double beta = f.x0 + i * f.binWidth;
      //relative phasespace width is proportional to sqrt(b) as E->0
      double factor = std::sqrt(beta);
      //Add also alpha factor: exp(-x)*x^n/n!:
      const double x = alpha2x*beta;//alpha = beta in the E->0 limit
      factor *= std::exp( -x + static_cast<double>(n)*std::log(x)
                          + minus_log_nfactorial);
      f_fmut[i] *= factor;
    }
  }

  //Find common grid for all of these functions:
  VectD grid;
  {
    //FIXME: should we instead use ncmax in the next call, and do additional
    //trimming when used for either beta or alpha?
    const double gridmax = ncmin( gnexpn.sabRange.x1(),
                                  gnexpn.sabRange.y1() );
    nc_assert_always( gridmax > 0.0 );
    std::vector<EquidistantGrid> allgrids;
    allgrids.reserve(fcts.size());
    for ( auto& f : fcts ) {
      allgrids.emplace_back();
      auto& g = allgrids.back();
      g.x0 = f.x0;
      g.binWidth = f.binWidth;
      nc_assert( g.binWidth > 0.0 );
      g.npts = f.f.size();
      trimEquidistantGridUpperEdge(g, gridmax);
      nc_assert( g.x1() <= gridmax );
    }

    grid = makeCommonGrid( allgrids );

    nc_assert_always(!grid.empty());
    nc_assert( grid.back() <= gridmax );
    nc_assert_always( grid.size() >= 10 );
    //ensure we have 0.0 in this:
    if ( grid.front() < 1e-3 * fcts.front().binWidth ) {
      grid.front() = 0.0;
    } else {
      //A bit expensive:
      grid.insert(grid.begin(), 0.0);
    }
  }
  nc_assert_always(grid.size() > 5 && grid.front() == 0.0 );
  nc_assert( nc_is_grid(grid) );

  //Add up all contributions on the common grid (here with no ad-hoc weights,
  //since we actually incorporated the real weight factors into the fct values
  //just above):
  auto contrib = evalPWLSum(fcts,grid);

  //Do not keep extreme tail points with negligible impact:
  trimTailByIntegral( grid, contrib, 1e-9 );

  nc_assert_always( grid.size() == contrib.size() );
  if ( npts < grid.size() )
    std::tie(grid, contrib) = reducePtsInDistribution( grid, contrib, npts );
  nc_assert_always( grid.size() <= npts );
  nc_assert_always( grid.size() == contrib.size() );
  nc_assert_always( grid.front() >= 0.0 );
  nc_assert_always( grid.front() >= gnexpn.sabRange.y0() );
  nc_assert_always( grid.back() <= gnexpn.sabRange.y1() );

  std::pair<VectD,VectD> res;
  res.first = std::move(grid);
  res.second = std::move(contrib);
  return res;
}

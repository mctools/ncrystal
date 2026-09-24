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
#include "NCrystal/internal/utils/NCMsg.hh"//SABXSDIAG-TEMPORARY
#include <iomanip>//SABXSDIAG-TEMPORARY

//TODO: reduce usage of nc_assert_always in this file once the model has been
//      used in production for a while.

namespace NC=NCrystal;

namespace NCRYSTAL_NAMESPACE {

  namespace VDOS {

    namespace {

      double combinedGnFctWeight( VDOSGn::Order n )
      {
        // We want to assign higher weight to lower n values, and highest of all
        // to n=1. w(n)=1/n is unstable since sum of all weights does not
        // converge, so we go for ~1/n^2 instead. We also choose a form where
        // the sum to n=inf is normalised, and the parameter we tune is the
        // weight of n=1. It is a bit of trial and error, and likely not super
        // crucial how this is setup. But the following seemed to give sensible
        // results.
        constexpr double p1 = 0.3;//1st order phonon weight
        constexpr double a = 1.0/p1-2.0;
        constexpr double ap1 = a + 1.0;
        return ap1 / ((n.value()+a)*(n.value()+ap1));
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

  //All Gn functions are on lattices which (mathematically) have a node at
  //beta=0, but numerically the node positions there are off by a few ulps
  //from 0, with random sign. Snap them to exactly 0, otherwise a random tiny
  //value (rather than exactly 0) can end up in the final beta grid, which
  //must contain exactly 0. Since points closer than 10% of the smallest bin
  //width are merged, there is at most a single point which can be affected.
  {
    double minbw = kInfinity;
    for ( auto& g : individual_grids )
      minbw = ncmin( minbw, g.binWidth );
    for ( auto& x : grid ) {
      if ( ncabs(x) < 1e-6 * minbw )
        x = 0.0;
    }
  }

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
    //
    //All comparisons with the edges have a tiny tolerance, since with all Gn on
    //lattices anchored at 0, points can be mathematically exactly at the edges
    //(which are given by nodes of the Gn), and whether they end up on one or
    //the other side of it should not depend on rounding errors.
    const double etol = 1e-9 * ( betaRange.second - betaRange.first );
    std::size_t i = 0;
    while ( bvals_view[i] < betaRange.first - etol )
      ++i;
    if ( i > 0 && bvals_view[i] > betaRange.first + etol )
      --i;//keep one point going over the edge.

    bvals_view = bvals_view.subspan(i);
    gnprojvals_view = gnprojvals_view.subspan(i);
    if ( bvals_view.back() > betaRange.second + etol ) {
      //Do the same for the upper limit, although this is expected to happen
      //only extremely rarely in usual operations.
      auto newsize = bvals_view.size();
      nc_assert_always( bvals_view.size() >= 5 );
      while ( newsize > 2
              && ncmin(bvals_view[newsize-1],
                       bvals_view[newsize-2]) > betaRange.second + etol ) {
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
    = reducePtsByEquidistribution( bvals_view, gnprojvals_view, nbeta );
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
        //Below this, no further points could be added anyway (and
        //1+rtol_try would eventually be indistinguishable from 1):
        if ( rtol_try < 1e-12 )
          NCRYSTAL_THROW(CalcError,"Unable to add enough points to"
                         " reach requested grid size");
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
  //SABXSDIAG-TEMPORARY: dedicated diagnostic commit, see tests/src/app_sabxsdiag.
  if ( std::getenv("NCRYSTAL_SABXS_DIAG") ) {
    NCRYSTAL_MSG( "SABXSDIAG alphaGrid size=" << avals.size()
                 << " betaGrid size=" << bvals.size() );
    for ( auto i : ncrange(avals.size()) )
      NCRYSTAL_MSG( "SABXSDIAG alphaGrid[" << i << "]="
                   << std::setprecision(17) << vectAt(avals,i) );
    for ( auto i : ncrange(bvals.size()) )
      NCRYSTAL_MSG( "SABXSDIAG betaGrid[" << i << "]="
                   << std::setprecision(17) << vectAt(bvals,i) );
  }
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
    if ( f.f.empty() ) {
      //Did not actually extend enough into positive values to have 2 grid
      //pts. Assume this will be true from now on.
      fcts.pop_back();
      break;
    }

    nc_assert_always( f.x0 > 0.0 && f.f.size() >= 2 );
    auto f_fmut = f.f_mutable();

    //Now add alpha and phase space factors:
    for (std::size_t i = 0; i < f.f.size(); ++i) {
      //note alpha=beta on the E->0 phasespace, so x=beta*alpha2x
      const double beta = f.xAt(i);
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
      nc_assert( g.x1() <= gridmax + 2e-9*g.binWidth );
    }

    grid = makeCommonGrid( allgrids );

    nc_assert_always(!grid.empty());
    //trimEquidistantGridUpperEdge has a tiny tolerance, so we might be
    //marginally above gridmax:
    if ( grid.back() > gridmax )
      grid.back() = gridmax;
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
    std::tie(grid, contrib) = reducePtsByEquidistribution( grid, contrib, npts );
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

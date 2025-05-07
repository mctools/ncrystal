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

#include "NCrystal/internal/vdos/NCVDOSToScatKnl.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/vdos/NCVDOSGn.hh"
#include "NCrystal/internal/phys_utils/NCKinUtils.hh"
#include "NCrystal/internal/utils/NCString.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/fact_utils/NCFactoryJobs.hh"
#include "NCrystal/internal/sab/NCSABUtils.hh"
namespace NC=NCrystal;

namespace NCRYSTAL_NAMESPACE {

  namespace V2SKDetail {

    namespace {

      static bool s_verbose = ncgetenv_bool("DEBUG_PHONON");
      inline double stirlingsSeriesSum9thOrder(double inv_n)
      {
        //Calculate and return the sum_{k=0}^9{ gk/x^k }, needed to estimate
        //Stirling's series (see eq. 5.11.3 and 5.11.4 in
        //https://dlmf.nist.gov/5.11): Thus, pass in inv_n = 1/n and get back
        //1+1/(12n)+1/(288n^2)+... (up to terms with 1/n^9) from Stirling's
        //series, The result must be multiplied by sqrt(2pi*n)*(n/e)^n to provide
        //an estimate of the faculty of n ("n!"). Factors ci are taken from:
        //https://oeis.org/A001163 and https://oeis.org/A001164
        //
        //For arguments n>=9 this can be used to evaluate n! to a relative
        //precision of O(1e-13) or better.
        constexpr double c1 = 1./12.;
        constexpr double c2 = 1./288.;
        constexpr double c3 = -139/51840.;
        constexpr double c4 = -571./2488320.;
        constexpr double c5 = 163879./209018880.;
        constexpr double c6 = 5246819./75246796800.;
        constexpr double c7 = -534703531./902961561600.;
        constexpr double c8 = -4483131259./86684309913600.;
        constexpr double c9 = 432261921612371./514904800886784000.;
        return 1.0 + inv_n*(c1+inv_n*(c2+inv_n*(c3+inv_n*(c4+inv_n*(c5+inv_n*(c6+inv_n*(c7+inv_n*(c8+inv_n*c9))))))));
      }

      constexpr double kTmsd_to_alpha2x(double kT,double msd)
      {
        //alpha2x is factor required to convert alpha value to x (aka "2W" in Sjolanders paper).
#if nc_cplusplus >= 201402L
        constexpr double fact = ( 2.0*const_neutron_mass_evc2/(constant_hbar*constant_hbar) );
        return fact*kT*msd;
#else
        //In c++11, constexpr functions can only consist of a single return statement.
        return (( 2.0*const_neutron_mass_evc2/(constant_hbar*constant_hbar) ))*kT*msd;
#endif
      }

      void nc_array_add_inplace( double * ncrestrict tgt,
                                 const double * ncrestrict src,
                                 std::size_t n )
      {
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] += src[i];
      }

      VectD fillSABFromVDOS( const VDOSGn& Gn_asym,
                             const double msd,
                             const VectD& alphaGrid,
                             const VectD& betaGrid,
                             ScaleGnContributionFct scaleGnContribFct,
                             unsigned min_order = 1,
                             unsigned max_order = std::numeric_limits<unsigned>::max() )
      {
        // Evaluate S(alpha,beta) from Sjolander's II.28, recasted to alpha/beta
        // and excluding sigma*kT/4E from the definition of S.

        const auto nalpha = alphaGrid.size();
        VectD sab( nalpha * betaGrid.size(), 0.0);
        const unsigned maxOrder = Gn_asym.maxOrder().value();
        const double kT = Gn_asym.kT();

        //Alpha-dependency is contained in f(x,n)=exp(-x)*x^n/n! where
        //x=alpha2x*alpha.  For n<stirling_threshold, f(x,n) are calculated
        //directly and recursively, using f(x,n)=f(x,n-1)*(x/n). For higher n,
        //this becomes numerically unstable and we instead use Stirling's series
        //to rewrite the factor of "n!" and evaluate f(x,n) directly. NB,
        //stirling_threshold should not be changed without validation: raise it
        //too much, and the direct method runs into issues at high n, decrease
        //it too much there will be precision (and speed) issues.
        //
        //The reason for picking 16 is it seems the stirling formula as
        //implemented here reaches "a few ulps" precision at n=16, and
        //presumably it is better to switch to it as early as possible for
        //numerical stability (although pure computation speed considerations
        //would imply switching as late as possible).
        constexpr unsigned stirling_threshold = 16;//See comment ^^^
        const double alpha2x = V2SKDetail::kTmsd_to_alpha2x(kT,msd);
        auto x_vals = vectorTrf(alphaGrid,[alpha2x](double alpha){ return alpha*alpha2x; } );
        auto expmhalfx_vals = vectorTrf(x_vals,[](double x){ return std::exp(-0.5*x); } );
        VectD logx_vals = ( maxOrder >= stirling_threshold ? vectorTrf(x_vals,[](double x){ return std::log(x); } ) : VectD() );
        auto fxn_cache = expmhalfx_vals;//We apply half of the exp(-x) factor before
        //the x^n/n! factor is added, and the other
        //half of exp(-x) we add later. This is done
        //to extend range of valid x-values.
        VectD alpha_factors;
        alpha_factors.resize(x_vals.size());

        //For reasons of numerical stability and efficiency, we always evaluate S
        //first for negative beta, and obtain the values at positive beta by the
        //identity ("detailed balance"): S(alpha,+beta)=S(alpha,-beta)*exp(-beta).
        //
        //Find range of betaGrid which should be flipped:
        auto itbeta_zero = findClosestValInSortedVector(betaGrid, 0.0);
        nc_assert( *itbeta_zero == 0.0 );
        auto itbeta_firstflip = findClosestValInSortedVector(betaGrid, -betaGrid.back());
        nc_assert( *itbeta_firstflip == -betaGrid.back() );
        const std::size_t idx_firstflip = std::distance(betaGrid.cbegin(),itbeta_firstflip);
        const std::size_t idx_zero = std::distance(betaGrid.cbegin(),itbeta_zero);
        nc_assert( idx_zero+1 <= betaGrid.size() );
        const VectD betaGridNonPositive( betaGrid.begin(), std::next(betaGrid.begin(),idx_zero+1) );//vector instead of span to simplify code below
        auto expbeta_vals_nonposbeta = vectorTrf( betaGridNonPositive, [](double beta){ return std::exp(beta); } );

        //Now we loop and fill the S-table. First over phonon order, n, next over
        //beta and finally alpha. We take care to keep as many calculations as
        //possible in the outer loops, while avoiding unnecessarily repeating
        //calculations or utilising enormous memory caches.

        for ( unsigned n = 1; n <= maxOrder; ++n) {
          if ( n < min_order )
            continue;
          if ( n > max_order )
            break;

          const double contribScaleFactor = scaleGnContribFct ? scaleGnContribFct(n) : 1.0;
          nc_assert(contribScaleFactor>=0.0);

          //Prepare for Gn(beta)-evaluations:
          auto betakT_Range = Gn_asym.eRange(n);
          const double invn = 1.0/n;

          //Preparation of f(x)=exp(-x)*x^n/n! is more tricky, due to reasons of
          //efficiency and numerical issues. As explained above, for orders below
          //stirling_threshold, we build up recursively, and above we use
          //Stirling's formula:
          if ( n < stirling_threshold ) {
            //for reasonably low orders we can build up slowly and cheaply (leaving
            //out of fxn_cache a final factor of expmhalfx for numerical stability):
            for (auto x : enumerate(x_vals) )
              vectAt(fxn_cache,x.idx) *= x.val*invn;
            for (auto fxn : enumerate(fxn_cache) )
              vectAt(alpha_factors,fxn.idx) = fxn.val * vectAt(expmhalfx_vals,fxn.idx) * kT;
          } else {
            //Very high order phonons, f(x) is non-zero at very high values of x, but
            //exp(-0.5*x) becomes 0, precluding the direct/cheap evaluation. Instead
            //we evaluate directly, with the help of Stirling's series for the
            //factorial, n!. Everything is suitably rearranged so cancellations happen
            //before the exponential is evaluated.
            const double gn = V2SKDetail::stirlingsSeriesSum9thOrder(invn);
            const double fact= kT * kInvSqrt2Pi/(std::sqrt(n)*gn);
            if (!fact)
              continue;//nothing can contribute at this order (should not really happen?)
            const double logn = std::log(n);
            for (auto x : enumerate(x_vals) ) {
              const double exparg = n * ( vectAt(logx_vals,x.idx) - logn + 1.0 ) - x.val;
              vectAt(alpha_factors,x.idx) = fact * std::exp(exparg);
            }
          }

          for ( auto beta : enumerate(betaGridNonPositive) ) {
            //Can evaluate more precisely at negative beta values and simply flip +
            //apply detailed balance factor for positive.
            nc_assert( beta.val<=0.0 );
            double energy = beta.val * kT;
            if (!valueInInterval(betakT_Range,energy))
              continue;//Gn(beta) zero here.

            double expMbeta(0.0);
            std::size_t posbeta_idx(0);
            if ( beta.idx >= idx_firstflip ) {
              posbeta_idx = idx_zero + ( idx_zero - beta.idx );
              expMbeta = vectAt(expbeta_vals_nonposbeta,beta.idx);
            }

            const double Gn_asym_eval = contribScaleFactor * Gn_asym.eval(n,energy);
            if ( !(Gn_asym_eval>0.0) )
              continue;
#if 0
            //readable version of innermost loop:
            const auto offset_alpharow = beta.idx*nalpha;
            const auto offset_alpharow_posbeta = posbeta_idx*nalpha;
            for ( auto alpha_fact : enumerate(alpha_factors) ) {
              if ( !(alpha_fact.val>0.0) )
                continue;
              double contrib_S_negbeta = alpha_fact.val * Gn_asym_eval;

              //Negative beta:
              vectAt( sab, offset_alpharow + alpha_fact.idx ) += contrib_S_negbeta;
              if ( expMbeta ) {
                //Expand to beta>0 by first copying S-values,
                //S(alpha,+beta)=S(alpha,-beta) and then applying the detailed
                //balance factor, exp(-beta/2) twice to the values at beta>0, thus
                //ensuring in the end that all entries end up with a correct factor:
                vectAt( sab, offset_alpharow_posbeta + alpha_fact.idx ) += contrib_S_negbeta * expMbeta;
              }
            }//alpha loop
#else
            //Attempt to super-streamline this innermost loop:
            nc_assert(nalpha==alpha_factors.size());
            double * itAlphaFactB = &alpha_factors[0];
            double * itAlphaFact = itAlphaFactB;
            double * itAlphaFactE = itAlphaFactB + nalpha;
            //Alpha-factors increase and then decrease. Thus, if we first skip
            //over any initial zeros in the alpha factors, we can break (rather
            //than just continue) in the final loop whenever we see a zero.

            //First skip over any initial zeros in alpha-factors:
            while( !(*itAlphaFact>0.0) && itAlphaFact!=itAlphaFactE )
              ++itAlphaFact;
            double * itSAB = &sab[0] + (beta.idx*nalpha + (itAlphaFact-itAlphaFactB));
            if ( expMbeta) {
              double * itSAB_posbeta = &sab[0] + (posbeta_idx*nalpha + (itAlphaFact-itAlphaFactB));
              for (;itAlphaFact!=itAlphaFactE;++itAlphaFact,++itSAB,++itSAB_posbeta) {
                if ( !(*itAlphaFact>0.0) )
                  break;
                double contrib_S_negbeta = *itAlphaFact * Gn_asym_eval;
                *itSAB += contrib_S_negbeta;
                *itSAB_posbeta += contrib_S_negbeta*expMbeta;
              }//alpha loop where +-|beta| is available
            } else {
              for (;itAlphaFact!=itAlphaFactE;++itAlphaFact,++itSAB) {
                if ( !(*itAlphaFact>0.0) )
                  break;
                *itSAB += *itAlphaFact * Gn_asym_eval;
              }//alpha loop where only -|beta| is available
            }
#endif
          }//beta loop
        }//phonon order loop

        return sab;
      }

      VectD fillSABFromVDOSConcurrent( const VDOSGn& Gn_asym,
                                       const double msd,
                                       const VectD& alphaGrid,
                                       const VectD& betaGrid,
                                       ScaleGnContributionFct scaleGnContribFct )
      {
        //Ideally, one would concurrently create the SAB from each order available
        //in Gn_asym, and then simply merge them afterwars. However, the memory
        //consumption required to do that might be very large. So instead, we can
        //put the processing of e.g. orders 1-10, 11-20, 21-30, etc. into separate
        //concurrent threads, and then add up the resulting SAB's
        //afterwards. However, if we do this we must ALWAYS use the same
        //subdivision of orders irrespective of how many threads are available,
        //since the summation of floating point numbers depend on the order in
        //which they are added: We do not want the number of threads available to
        //change results, and we want an increase in maxorder to only affect the
        //high-E region.
        //
        //Also, be aware that sab.size()=nalpha*nbeta might be huge, and each
        //concurrent job will need its own copy of it, so we should not make
        //njobs too huge!
        const unsigned norders = Gn_asym.maxOrder().value();
        const unsigned njobs = ( norders <= 16 ? 1 : norders / 16 );
        if ( njobs == 1 )
          return fillSABFromVDOS( Gn_asym, msd,
                                  alphaGrid, betaGrid, scaleGnContribFct );
        const unsigned norders_per_job = norders / njobs;
        nc_assert_always( norders_per_job >= 1 );

        SmallVector<VectD,16> results;
        results.resize(njobs);
        unsigned nextorder = 1;
        FactoryJobs jobs;
        for ( auto ijob : ncrange(njobs) ) {
          const unsigned min_order = nextorder;
          nextorder += norders_per_job;
          const unsigned max_order = std::min<unsigned>( nextorder-1, norders );
          nc_assert_always( max_order != norders || ijob+1 == njobs );
          nc_assert_always(min_order >= 1);
          nc_assert_always(max_order >= min_order);
          nc_assert_always(max_order <= norders);
          VectD * resptr = &results.at(ijob);
          jobs.queue([resptr,min_order,max_order,
                      &Gn_asym, msd, &alphaGrid, &betaGrid,&scaleGnContribFct]
                     ()
          {
            *resptr = fillSABFromVDOS(Gn_asym,msd,
                                      alphaGrid,betaGrid,scaleGnContribFct,
                                      min_order, max_order );
          });
        }
        jobs.waitAll();

        //Add up results (:
        VectD res = std::move(results.at(0));
        nc_assert_always( res.size() > 0 );
        for ( unsigned ijob = 1; ijob < njobs; ++ijob ) {
          VectD& src = results.at( ijob );
          nc_assert_always( res.size() == src.size() );
          nc_array_add_inplace( &*res.begin(), &*src.begin(), res.size() );
          src.clear();
        }
        return res;
      }

    }
  }
}

NC::PairDD NC::rangeXNexpMX(unsigned n, double eps, double accuracy ) {
  //Interval where f(x) = x^n*exp(-x) is above eps*fpeak.

  nc_assert(eps>0.0&&eps<1.0&&eps>1e-200&&n>0&&accuracy>0&&accuracy<=1e-2);

  //The function f(x) = x^n*exp(-x) peaks at x=n and falls off on both
  //sides. Returns the two solutions to f(x)= f(n)*eps, describing the central
  //range around x=n where the function is higher than eps times the peak
  //value.

  //Must solve for x:
  //  x^n*exp(-x) = eps*[n^n*exp(-n)]
  //Raise to 1/n power and get:
  //  x*exp(-x/n) = eps^(1/n)*n*exp(-1)
  //<=> (x/n)*exp(-x/n) =  (1/e)*eps^(1/n) = k
  //
  //Which can be solved numerically for x/n:

  const double fn = static_cast<double>(n);
  const double k = kInvE * std::pow(eps,1.0/fn);
  auto f = [k](double y) { return y*std::exp(-y)-k; };
  return { fn*findRoot2( f, 0.0,   1.0, accuracy ),
           fn*findRoot2( f, 1.0, 700.0, accuracy ) };
}

NC::PairDD NC::findExtremeSABPointWithinAlphaPlusCurve(double E_div_kT, PairDD alphaRange, PairDD betaRange)
{
  //Find the extreme (as in highest alpha, lowest beta) kinematically
  //accessible point in the provided rectangular region in (alpha,beta) space
  //for a neutron with energy/kT= Emax_div_kT. Returns {-1,-1} in case no point
  //is accessible. Note that we on purpose consider only the kinematic edge
  //given by the alpha+(beta) and beta=-E/kT curves, ignoring the alpha-(beta)
  //curve.
  nc_assert( alphaRange.second > alphaRange.first );
  nc_assert( alphaRange.first >= 0.0 );
  nc_assert( betaRange.second > betaRange.first );
  nc_assert( E_div_kT > 0.0 );

  NCRYSTAL_DEBUGONLY(const bool should_be_accessible = sabPointWithinAlphaPlusCurve(E_div_kT,alphaRange.first,betaRange.second));

  if ( betaRange.second <= -E_div_kT ) {
    nc_assert(!should_be_accessible);
    return { -1.0, -1.0 };//no accessible points in region
  }

  auto alphaPlus = [E_div_kT](double beta)
                   {
                     nc_assert( beta >= -E_div_kT );
                     return 2*E_div_kT + beta + 2 * std::sqrt( E_div_kT * ( E_div_kT + beta ) );
                   };
  const double apb1 = alphaPlus(betaRange.second);
  if ( apb1 <= alphaRange.first ) {
    nc_assert(!should_be_accessible);
    return { -1.0, -1.0 };//no accessible points in region
  }

  nc_assert(should_be_accessible);

  //Clip lower beta range at -E/kT:
  betaRange.first = ncmax( betaRange.first, -E_div_kT );

  const double apb0 = alphaPlus(betaRange.first);
  if ( apb0 >= alphaRange.second )
    return { alphaRange.second, betaRange.first };//entire rectangle is accessible

  //Cut away excess reach of rectangular region along alpha:
  alphaRange.second = ncmin(alphaRange.second,apb1);

  //Cut away excess reach of rectangular region along beta:
  if ( apb0 < alphaRange.first ) {
    //The next formula follows from inverting the formula
    //alphaPlus(betaRange.first) = alphaRange.first:
    betaRange.first = alphaRange.first - 2.0 * std::sqrt( E_div_kT * alphaRange.first );
    nc_assert( floateq( alphaPlus(betaRange.first), alphaRange.first ) );
  }

  //Rectangular region has no excess now, result is given by its extreme
  //corner:
  return { alphaRange.second, betaRange.first };
}

bool NC::sabPointWithinAlphaPlusCurve(double E_div_kT, double alpha, double beta )
{
  //Same as findExtremeSABPointWithinAlphaPlusCurve, but only testing whether or
  //not a given single point is accessible, in the sense that it has beta>=-E/kT
  //and alpha<alpha+(beta). Note that this deliberately ignores the alpha-(beta)
  //curve.
  nc_assert( alpha >= 0.0 );
  nc_assert( E_div_kT > 0.0 );
  const double c = E_div_kT;
  const double cpb = c + beta;
  if ( cpb < 0.0 )
    return false;
  //With c=E/kT, we check the alpha+ condition:
  //  alpha+(beta) = 2*c + beta + 2 * sqrt( c * ( c + beta ) ).
  //  So "within" means: alpha+(beta) >= alpha
  //  <=> sqrt( c * ( c + beta ) ) >= (alpha-beta)/2 - c
  const double t = 0.5 * ( alpha - beta ) - c;
  return t <= 0.0 || c*cpb >= t*t;
}

NC::VectD NC::setupAlphaGrid( double kT, double msd, double alphaMax, unsigned npts )
{
  nc_assert(npts>=20);

  //Strategy (region 0 ensures proper binning for upscattering of ultra-low
  //energy neutrons, and are simply merged onto the grid pts from the other 3
  //regions at the end, using the finalise_grid function further down):
  //
  // Region 0 (N=15% of pts): x=linspace from x=1e-3 to x=alpha_upscatmax*alpha2x
  // Region 1 (N=29% of pts): x=1e-50 + linspace[N-1] from x=1e-10 to x=1.0 (including 1.0).
  // Region 2 (N=23% of pts): linspace[N+2] from x=1.0 to x=15.0 (excluding both endpoints).
  // Region 3 (rest, N~=33% of pts): geomspace[N+1] from x=15.0 to final limit (including 15.0).
  //
  //We give special treatment to cases where alphaMax is not above alphaG15Maxx.
  //
  //alpha_upscatmax is roughly 10, but slightly higher or lower depending on
  //number of points.

  const double alpha2x = V2SKDetail::kTmsd_to_alpha2x(kT,msd);
  const double x2alpha = 1.0 / alpha2x;
  const double alphaMin = x2alpha*1e-50;
  const double alphaMin2 = x2alpha*1e-10;
  const double alphaG1Maxx = x2alpha*1.0;
  const double alphaG15Maxx = x2alpha*15.0;

  unsigned n0 = static_cast<unsigned>(npts*0.15+0.5);
  unsigned n1 = static_cast<unsigned>(npts*0.29+0.5);
  unsigned n2 = static_cast<unsigned>(npts*0.23+0.5);
  unsigned n3 = npts-(n1+n2+n0);
  unsigned npts_123 = n1+n2+n3;
  nc_assert(n1>=3&&n2>=3&&n3>=3&&n0>=3&&n0+n1+n2+n3==npts);

  const double alpha_upscatmax = (n0<10?6.0:(n0>50?14.0:10.0));
  const auto grid_region0 = linspace( ncmin(1e-3,alphaMax*0.01),
                                      ncmin(alpha_upscatmax,alphaMax*0.99),
                                      n0 );

  auto finalise_grid = [&grid_region0,npts](const VectD& grid) -> VectD
                       {
                         nc_assert_always(grid.size()+grid_region0.size()==npts);

                         //We want to merge the grid_region0 points into the
                         //grid points. However, the grid_region0 points will be
                         //moved a bit to not accidentally coincide with
                         //existing points. Thus, for each grid_region0 pt, find
                         //the two neighbouring values (looking in *both* grid
                         //and grid_region0). Place the point exactly between
                         //those two. I.e. just merge+sort all numbers, but make
                         //sure we know which are from region0. Those will be
                         //placed halfway between neighbours afterwards.
                         std::vector<std::pair<double,bool>> all_merged;
                         for (auto e: grid)
                           all_merged.emplace_back(e,false);
                         for (auto e: grid_region0)
                           all_merged.emplace_back(e,true);
                         std::stable_sort(all_merged.begin(),all_merged.end());
                         auto it = std::next(all_merged.begin());
                         auto itLast = std::prev(all_merged.end());
                         for (;it!=itLast;++it) {
                           if ( it->second == false )
                             continue;
                           it->first = 0.5*(std::prev(it)->first+std::next(it)->first);
                         }
                         VectD out;
                         out.reserve(npts);
                         for (const auto& e: all_merged)
                           out.push_back(e.first);
                         nc_assert_always(nc_is_grid(out));
                         nc_assert_always(out.size()==npts);
                         return out;
                       };

  if ( alphaMax <= alphaMin*100.0 ) {
    //special case #1, insanely small alphaMax
    return finalise_grid(linspace(alphaMax*0.001,alphaMax,npts_123));
  }

  VectD grid;
  grid.reserve(npts_123);
  grid.push_back(alphaMin);
  if ( alphaMax <= alphaG1Maxx*10.0 ) {
    //special case #2, small alphaMax.
    vectorAppend(grid, linspace(alphaMin2,alphaMax,npts_123-1) );
    nc_assert( grid.size() == npts_123 && nc_is_grid(grid) );
    return finalise_grid(grid);
  }

  //Region 1:
  vectorAppend(grid, linspace(alphaMin2,alphaG1Maxx,n1-1) );

  //Region 2:
  if ( alphaMax < 2.0 * alphaG15Maxx) {
    //Special case #3: somewhat small alphaMax => absorb region3 into region 2
    auto ls2 = linspace(alphaG1Maxx,alphaMax,n2+n3+2);
    grid.insert( grid.end(), std::next(ls2.begin()), std::prev(ls2.end()) );
    nc_assert( grid.size() == npts_123 && nc_is_grid(grid) );
    return finalise_grid(grid);
  }
  auto ls2 = linspace(alphaG1Maxx,alphaG15Maxx,n2+2);
  grid.insert( grid.end(), std::next(ls2.begin()), std::prev(ls2.end()) );
  //Region 3:
  vectorAppend(grid, geomspace(alphaG15Maxx,alphaMax,n3) );
  nc_assert( grid.size() == npts_123 && nc_is_grid(grid) );
  return finalise_grid(grid);
}

NC::VectD NC::setupBetaGrid( const NC::VDOSGn& Gn, double betaMax, unsigned vdoslux, unsigned override_nbins )
{
  nc_assert(Gn.maxOrder().value()>=1);
  nc_assert_always(vdoslux<=5);
  nc_assert_always(betaMax>0.0);

  //Reach of 1-phonon and 3-phonon spectrums (converted E->beta):
  const double invkT = 1.0/Gn.kT();
  const double G1 = ncabs(Gn.eRange(1).first*invkT);
  nc_assert(G1>0.0);
  double G3 = ncabs(Gn.eRange( std::min<unsigned>(3,Gn.maxOrder().value()) ).first*invkT);//nb: actually =G1 or =G2 in case of maxOrder<3.
  nc_assert( Gn.maxOrder().value()==1 || G3>G1 );
  if (G1==G3)
    G3 = G1*1.0001;


  //Construct betagrid which will run from [-betaMax,D] with
  //0<D<=betaMax. Assuming betaMax is large enough, D will always be large
  //enough to encompass the first phonon order, and energy transfers of 20*kT,
  //whichever is larger (at luxlevel 4 and 5 the 20*kT limit is increased to
  //30*kT and 40*kT respectively).
  //
  //We divide pts in regions:
  //Region0: The single point beta=0.
  //Region1: [-G1,0) (all points flipped and duplicated in (0,G1]).
  //Region2: [-D,-G1) (all points flipped and duplicated in (G1,D]).
  //Region3: [-betaMax,-D) (no duplication of points). Must be odd number of points.
  //
  //always encompass single-phonons, and clip G3/D to betaMax if needed, to
  //ensure no regions will have 0 range:
  betaMax = ncmax(betaMax,G1*1.01);
  G3 = ncmin(G3,betaMax*0.9999);
  const double D = ncmin(betaMax*0.9999,ncmax(G3,std::max<int>(2,vdoslux-1)*10.0));
  nc_assert( D>0.0 && G3 > G1 && D>=G3 && D<=betaMax );

  //How many points to use? The total is given as 100 * 2^(vdoslux), meaning
  //100 for vdoslux=0, 200 for vdoslux=1, and so on up to 3200 for vdoslux=5:
  const unsigned ntotal = override_nbins ? override_nbins : 100*(1<<vdoslux);//100 * 2^(vdoslux).

  VectD grid;
  grid.reserve(ntotal);

  //First do a bit of pre-analysis for the G1 spectrum in Region1:
  //Here we don't just space the points uniformly, but try to place them to
  //get the single-phonon curve best described. Additionally, we use 20% (but
  //at least 5) of the points to put some very small values around
  //beta=0. This is not needed to describe the actual data, but is used to
  //reduce artefacts in some integration/sampling algorithms (ideally those
  //algs should instead be updated to avoid such artefacts all by
  //themselves...).

  auto rawspec_G1 = Gn.getRawSpectrum(1);
  auto erange_G1 = Gn.eRange(1);
  //Convert E->beta:
  erange_G1.first *= invkT;
  erange_G1.second *= invkT;
  nc_assert(floateq(erange_G1.first,-G1));
  auto evals_G1 = linspace(erange_G1.first,erange_G1.second,rawspec_G1.size());
  nc_assert_always(evals_G1.front()<0.0);
  nc_assert_always(evals_G1.back()>0.0);
  //Ignore positive part and 0 (NB: value representing zero can actually be
  //slightly non-zero due to numerical issues, hence need for epsilon):
  const double epsilon = -0.1*Gn.binWidth(1);
  nc_assert_always(evals_G1.front()<epsilon);
  while ( evals_G1.back()>epsilon )
    evals_G1.pop_back();
  nc_assert(evals_G1.size()>=2);
  rawspec_G1.resize(evals_G1.size());

  //The absolute maximum number of points which can be useful in region 1 is
  //evals_G1.size() plus a bit for fine-grained grid near beta=0:
  const unsigned n1_max = static_cast<unsigned>(evals_G1.size()*1.25+6.5) & ~1;//& ~1 ensures even
  nc_assert(n1_max>=8);

  //Determine how many pts to assign to the different regions, depending on the
  //length of those regions:

  const unsigned n0 = 1;//the beta=0 point

  //Length of beta-scale covered by various regions:
  nc_assert_always(D>G1);
  nc_assert_always(betaMax>D);
  double L1 = 2*G1;
  double L2 = 2*(D-G1);
  double L3 = betaMax-D;
  //Ad-hoc boost of central parts where we expect more variation:
  L1 *= 4;
  L2 *= 2;
  double Lnorm = 1.0/(L1+L2+L3);
  double f1 = ncmax(0.2,L1*Lnorm);
  double f2 = L2*Lnorm;
  if (f1+f2>0.99) {
    double tmp = 0.99/(f1+f2);
    f1 *= tmp;
    f2 *= tmp;
  }

  unsigned n1 = std::min<unsigned>(n1_max,static_cast<unsigned>(f1*ntotal+0.5) / 2);//div 2 due to duplication
  unsigned n2 = static_cast<unsigned>(f2*ntotal+0.5) / 2;//div 2 due to duplication
  n1 = std::max<unsigned>(15,n1);
  n2 = std::max<unsigned>(1,n2);
  unsigned n012;
  do {
    n012 = 2*(n1+n2)+n0;
    if (n012 >= ntotal - 1) {
      if (n2>n1)
        --n2;
      else
        --n1;
    }
    else
      break;
  } while (true);
  nc_assert_always(n1>=10);

  unsigned n3 = ntotal-n012;

  nc_assert(n0==1);
  nc_assert(n1>=10);
  nc_assert(n2>=1);
  nc_assert(n3>=1);
  nc_assert( ntotal == n0 + n3 + 2 * ( n1 + n2) );
  nc_assert_always(n1 <= ntotal );
  nc_assert_always(n2 <= ntotal );
  nc_assert_always(n3 <= ntotal );
  nc_assert_always(D>G1);
  nc_assert_always(betaMax>D);

  //Time to add values to the grid! First Region3:
  {
    auto vals = linspace(-betaMax, -D, n3+1);
    grid.insert(grid.begin(),vals.begin(),std::prev(vals.end()));
    nc_assert(grid.size()==n3);
  }

  const auto idx_R2Start = grid.size();

  //Then Region 2 - negative values.
  {
    auto vals_r2 = linspace(-D,-G1,n2+1);
    vals_r2.resize(vals_r2.size()-1);
    vectorAppend(grid,vals_r2);
  }

  //Now Region 1 - negative values.
  {
    unsigned n1_near0 = std::max<unsigned>(5,static_cast<unsigned>(n1*0.2+0.5));
    unsigned n1_spectrum = n1-n1_near0;
    if (n1_spectrum >= evals_G1.size()) {
      //Great: can keep all points in G1 spectrum:
      n1_spectrum = evals_G1.size();
      n1_near0 = n1 - n1_spectrum;
    } else {
      //Remove least important points from G1 spectrum:
#if 0
      //Old way, just run reducePtsInDistribution. This might leave huge gaps if
      //for instance VDOS is zero over a wide internal region:
      std::tie(evals_G1,rawspec_G1) = reducePtsInDistribution( evals_G1,rawspec_G1, n1_spectrum );
#else
      //New way, avoid leaving huge gaps without pts. Doing so can leave large
      //gaps in the beta points, which when can lead to numerical artifacts when
      //a scattering kernel is later integrated and sampled. If at some point
      //(TODO!) we implement better integration/sampling algorithms, we can
      //hopefully stop doing this - and perhaps also stop needing the n1_near0
      //values.
      unsigned n_for_gaps = ( ( n1_spectrum > 30 && evals_G1.size()-n1_spectrum > 10 )
                              ? std::max<unsigned>(5,static_cast<unsigned>(n1_spectrum*0.1+0.5))
                              : 0 );
      n1_spectrum -= n_for_gaps;
      const double G1binwidth = evals_G1.at(1)-evals_G1.at(0);
      std::tie(evals_G1,rawspec_G1) = reducePtsInDistribution( evals_G1,rawspec_G1, n1_spectrum );
      if ( n_for_gaps > 0 ) {
        struct Gap {
          Gap(double bbb0, double bbb1)
            : b0(bbb0), b1(bbb1)
          {}
          double b0, b1;//edges of gap
          unsigned n = 0;//number of points allocated to fill the gap
          bool operator<( const Gap& o ) const {
            double a = ( b1 - b0 ) / ( n + 1 );
            double b = ( o.b1 - o.b0 ) / ( o.n + 1 );
            if ( floateq(a,b,1e-13,1e-13) )
              return b0 > o.b0;
            return a > b;//reverse sort
          }
        };
        //Find all gaps:
        std::vector<Gap> gaps;
        gaps.reserve(256);
        double dmin = 1.5 * G1binwidth;
        auto it = evals_G1.begin();
        auto itLast = std::prev(evals_G1.end());
        for ( ; it != itLast; ++it ) {
          double dd = *std::next(it)-*it;
          if ( dd > dmin )
            gaps.emplace_back( *it, *std::next(it) );
        }
        //Allocate pts one by one into largest remaining gap:
        while (n_for_gaps>0 && !gaps.empty() ) {
          std::stable_sort(gaps.begin(),gaps.end());
          gaps.front().n += 1;
          --n_for_gaps;
        }
        //Transfer points in gaps to evals_G1:
        for ( const auto& g : gaps ) {
          if ( g.n > 0 ) {
            double bw = (g.b1-g.b0)/(g.n+1.0);
            for ( auto i : ncrange(g.n) ) {
              evals_G1.push_back( g.b0 + (i+1)*bw );
            }
          }
        }
        std::sort(evals_G1.begin(),evals_G1.end());
      }
      if ( n_for_gaps > 0 ) {
        //unused gap filler points, put somewhere else:
        n1_near0 += n_for_gaps;
      }
#endif
    }

    for ( auto e : evals_G1 )
      grid.push_back(e);

    //To reduce artefacts in sampling/integration algs, add points near 0:
    nc_assert( n1_near0 > 0 );
    nc_assert( grid.back() < 0.0 );
    auto v = geomspace(ncmin(1e-50,-0.001*grid.back()),-grid.back()*0.1,n1_near0);
    std::reverse(v.begin(),v.end());
    for (auto e: v)
      grid.push_back( -e );
  }
  const auto idx_R1Back = grid.size()-1;

  //Then Region0, just the beta=0 point:
  grid.push_back( 0.0 );

  //Now replay region1 and region2 in flipped reverse order to generate
  //positive beta values:
  for ( auto i = idx_R1Back; i>=idx_R2Start; --i )
    grid.push_back( -vectAt(grid,i) );

  nc_assert( grid.front()==-betaMax );
  nc_assert( grid.back()==D );
  nc_assert( grid.size() == ntotal );
  nc_assert( nc_is_grid(grid) );
  return grid;
}

NC::ScatKnlData NC::createScatteringKernel( const VDOSData& vdosdata,
                                            unsigned vdoslux,
                                            double targetEmax_requested,
                                            VDOSGn::TruncAndThinningParams ttpars,
                                            ScaleGnContributionFct scaleGnContributionFct,
                                            Optional<unsigned> call_override_max_order )
{
  //Hidden unofficial env-vars used for special debugging purposes:
  const unsigned override_max_order = ( call_override_max_order.has_value()
                                        ? call_override_max_order.value()
                                        : static_cast<unsigned>(ncgetenv_int("HACK_MAXORDER")) );
  const double override_alphamax = ncgetenv_dbl("HACK_ALPHAMAX");
  const double override_betamax = ncgetenv_dbl("HACK_BETAMAX");
  const unsigned override_nbins = ncgetenv_int("HACK_NBINS");

  //Which Emax should we target (i.e. aim to cover the kinematic reachable area
  //for neutrons of that energy):
  nc_assert_always( vdoslux <= 5 );
  nc_assert_always(targetEmax_requested>=0.0);
  constexpr double lux2emax[6] = { 0.5, 1.0, 3.0, 5.0, 8.0, 12.0 };//Emax in eV for vdosluxs 0 to 5
  nc_assert_always( vdoslux < 6 );
  double targetEmax = targetEmax_requested>0.0 ? targetEmax_requested : lux2emax[vdoslux];

  if ( V2SKDetail::s_verbose )
    NCRYSTAL_MSG("VDOS2SK initialising with T="<<vdosdata.temperature()
                 <<", vdoslux="<<vdoslux
                 <<", aiming for Emax="<<targetEmax<<"eV"
                 <<(targetEmax_requested>0.0?" (as requested)":"")<<", ...");

  //Initialise evaluators:
  VDOSEval vdoseval(vdosdata);
  const double kT = vdoseval.kT();
  const double invkT = 1.0/kT;
  const double gamma0 = vdoseval.calcGamma0();
  const double msd = vdoseval.getMSD( gamma0 );
  double targetEmax_div_kT = targetEmax*invkT;
  unsigned max_phonon_order = std::max<unsigned>(override_max_order,4);
  VDOSGn Gn_asym(vdoseval,ttpars);
  Gn_asym.growMaxOrder(max_phonon_order);

  //What are the highest phonon order we allow? When user requested given target
  //Emax, allow higher order expansions to accommodate. Otherwise, keep low to
  //keep default initialisation times reasonable (with special extreme settings
  //for vdoslux 0 and 5).
  unsigned order_limit = 1000;
  if ( targetEmax_requested > 0.0 || vdoslux == 5 )
    order_limit *= 10;
  if ( vdoslux==0 )
    order_limit /= 10;

  const double emax_lowest_allowed = ( targetEmax_requested>0.0 ? targetEmax_requested : 1e-15 );

  //Now increase order dynamically until the last order only has contributions
  //to S(alpha,beta) outside the kinematic reach of Emax:
  const double relcontriblvl = std::pow(10.0,-(3.0+2.0*vdoslux));//e.g.: 1e-3 for vdoslux 0, 1e-9 for vdoslux 3, 1e-13 for vdoslux 5
  const double x2alpha = 1.0 / V2SKDetail::kTmsd_to_alpha2x(kT,msd);
  auto findAlphaBetaRangeOfOrder = [&Gn_asym,x2alpha,invkT,relcontriblvl](unsigned n) {
                                     auto eRange = Gn_asym.eRange(n, relcontriblvl);
                                     PairDD betaRange( eRange.first * invkT, eRange.second * invkT  );
                                     auto xRange = rangeXNexpMX( n, relcontriblvl );
                                     PairDD alphaRange( xRange.first * x2alpha, xRange.second * x2alpha  );
                                     return std::make_pair(alphaRange,betaRange);
                                   };
  while (true) {
    if (override_max_order>0)
      break;
    Gn_asym.growMaxOrder(max_phonon_order);
    PairDD alphaRange, betaRange;
    std::tie(alphaRange, betaRange) = findAlphaBetaRangeOfOrder(Gn_asym.maxOrder().value());
    if (sabPointWithinAlphaPlusCurve(targetEmax_div_kT,alphaRange.first,betaRange.second)) {
      ++max_phonon_order;//could consider larger stepsize, but need to carefully check usage in the following
    } else {
      break;
    }
    if (max_phonon_order>order_limit) {
      //Too slow - unfeasible to fill out S(alpha,beta) all the way out to the
      //kinematic curve for E=targetEmax. In this case it is better to reduce
      //targetEmax, to at least get a consistent table (and hope the free-gas
      //extrapolation mechanisms will be adequate already at this lower
      //threshold).
      double targetEmax_reduced  = targetEmax;
      do {
        targetEmax_reduced *= 0.99;
        if ( targetEmax_reduced < emax_lowest_allowed )
          NCRYSTAL_THROW2(CalcError,"VDOS expansion too slow - can not reach E="<<emax_lowest_allowed
                          <<"eV after "<<order_limit<<" phonon convolutions (likely causes: either the target energy"
                          " value is too high, vdoslux too low, the temperature too high, or the VDOS is very unusual).");
      } while (sabPointWithinAlphaPlusCurve(targetEmax_reduced*invkT,alphaRange.first,betaRange.second));
      if (V2SKDetail::s_verbose)
        NCRYSTAL_WARN("VDOS2SK Could only reach Emax="<<targetEmax_reduced<<"eV and not the requested Emax="<<targetEmax<<"K");
      targetEmax_div_kT = targetEmax_reduced * invkT;
      targetEmax = targetEmax_reduced;
      break;
    }
  }
  nc_assert_always( targetEmax_requested==0.0 || targetEmax_requested == targetEmax );
  Gn_asym.growMaxOrder(max_phonon_order);

  //Ok, we now know how many orders we need to reach targetEmax. Next step is to
  //look at the contribution of each order insided the kinematic reach of
  //targetEmax, and use it to determine alpha/beta limits:

  double betaMin = 0.0;
  double alphaMax = 0.0;
  for ( unsigned n = 1; n<=max_phonon_order; ++n ) {
    PairDD alphaRange, betaRange;
    std::tie(alphaRange, betaRange) = findAlphaBetaRangeOfOrder(n);
    auto ep = findExtremeSABPointWithinAlphaPlusCurve(targetEmax_div_kT, alphaRange, betaRange);
    alphaMax = ncmax(alphaMax,ep.first);
    betaMin = ncmin(betaMin,ep.second);
  }
  nc_assert_always(betaMin<0.0 && alphaMax > 0.0);
  double upper_beta = -betaMin*1.01;
  double upper_alpha = alphaMax*1.01;
  if (override_alphamax)
    upper_alpha = override_alphamax;
  if (override_betamax)
    upper_beta = override_betamax;
  nc_assert_always( upper_beta>0.0 && upper_alpha>0.0 );

  //Ok, time to setup the alpha/beta grids. The grid-spacing is not even, rather
  //it attempts to best accomodate features of the distributions:
  VectD betaGrid = setupBetaGrid( Gn_asym, upper_beta, vdoslux, override_nbins );
  const unsigned alpha_size = ( override_nbins ? override_nbins : betaGrid.size()/2 );
  VectD alphaGrid = setupAlphaGrid( kT, msd, upper_alpha, alpha_size );

  //All done, now all that remains is to go through the (alpha,beta) pts in the
  //grid and use Sjolander's II.28 equation to calculate S(alpha,beta) there as
  //the sum of individual phonon orders:
#if 0
  //old way, no concurrency:
  auto sab = V2SKDetail::fillSABFromVDOS( Gn_asym, msd, alphaGrid, betaGrid, scaleGnContributionFct );
#else
  auto sab = V2SKDetail::fillSABFromVDOSConcurrent( Gn_asym, msd, alphaGrid, betaGrid, scaleGnContributionFct );
#endif
  double suggestedEmax = ( override_max_order>0 ? 0.0 : targetEmax);
  if ( scaleGnContributionFct!=nullptr && scaleGnContributionFct(max_phonon_order) == 0.0 ) {
    //Caller might have essentially removed the last order(s), so it is unknown
    //how far the kernel can be used.
    suggestedEmax = 0.0;
  }

  if (V2SKDetail::s_verbose)
    NCRYSTAL_MSG("VDOS2SK created SK with vdos expansion order N="<<max_phonon_order
                 <<", Emax="<<targetEmax<<"eV, nalpha="<<alphaGrid.size()<< " nbeta="<<betaGrid.size());

  ScatKnlData out;
  out.alphaGrid = std::move(alphaGrid);
  out.betaGrid  = std::move(betaGrid);
  out.sab       = std::move(sab);
  out.temperature = vdoseval.temperature();
  out.boundXS = vdosdata.boundXS();
  out.elementMassAMU = vdosdata.elementMassAMU();
  out.knltype = ScatKnlData::KnlType::SAB;
  out.suggestedEmax = suggestedEmax;
  out.betaGridOptimised = true;//prevent beta-thickening code upon conversion to SABData
  return out;
}

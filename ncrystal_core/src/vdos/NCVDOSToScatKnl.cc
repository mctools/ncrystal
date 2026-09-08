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

#include "NCrystal/internal/vdos/NCVDOSToScatKnl.hh"
#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/vdos/NCVDOSGn.hh"
#include "NCrystal/internal/vdos/NCVDOSKnlGrid.hh"
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

      //sorts grid, removes duplicates (within floateq) and ensures that pts are
      //spaced at least a factor of 1.01 apart.
      void sortAndCleanupGrid(VectD& grid_candidates, double tol = 1.01, double tolmax = 0.0 )
      {
        //        NCRYSTAL_MSG("sortAndCleanupGrid tol="<<tol);
        //fixme: also add extra pts in regions where there is too far between points? For instance the pts [0,0.01,2,..] has a factor 200 between pt 2 and 3.
        //fixme: unit test
        VectD& v = grid_candidates;
        std::sort(v.begin(),v.end());
        const std::size_t n = v.size();
        nc_assert(n>=2);
        nc_assert(tol >= 1.0001);
        nc_assert(v.front() >= 0.0);
        nc_assert_always(v.back() > tol * v.front());
        auto check2=[tol](double a,double b)
        {
          nc_assert_always( b>=0.0 );
          if ( !( b>tol*a ) )
            NCRYSTAL_THROW(BadInput,"VDOS grid: Too tight spacing detected.");
        };
        if ( n==2 ) {
          check2(v.back(),v.front());
          return;
        }


        VectD tmp;
        tmp.reserve(n);
        tmp.push_back(v.front());


        const std::size_t nm1 = n-1;
        for (std::size_t k = 1; k < nm1; ++k) {
          const double curr = v[k];
          const double prev = tmp.back();
          const double required = tol * prev; // need next >= tol*prev
          if (curr >= required) {
            if (!floateq(curr, prev))
              tmp.push_back(curr);
            continue;
          }

          //violation, move up to required, if there is space.
          const double nextCand = v[k + 1];
          if ( nextCand >= tol * required ) {
            if (!floateq(required, prev))
              tmp.push_back(required);
          } else {
            // drop
          }
        }

        // Always append last point unchanged (possible discarding the current
        // last pt if duplicate):
        if ( floateq(v.back(), tmp.back()) ) {
          nc_assert(tmp.size()>=2);
          tmp.back() = v.back();
        } else {
          tmp.push_back(v.back());
        }

        std::swap(tmp,v);

        if ( v.size()<2 )
          NCRYSTAL_THROW(BadInput,"VDOS Grid trimmed too hard");
        if ( v.size()==2 ) {
          check2(v.back(),v.front());
          return;
        }

        //Verify tol is OK for last pt, otherwise wiggle or remove pt before:
        if ( !(v.back()>tol*vectAt(v,v.size()-2) ) ) {
          double t = v.back()/tol;
          if ( t > tol*vectAt(v,v.size()-3) ) {
            NCRYSTAL_MSG("Fix 2nd to last");
            vectAt(v,v.size()-2) = t;
          } else {
            NCRYSTAL_MSG("discard 2nd to last");
            vectAt(v,v.size()-2) = v.back();
            v.pop_back();
            if ( v.size()==2 ) {
              check2(v.back(),v.front());
              return;
            }
          }
        }
        v.shrink_to_fit();

        Optional<double> dbg_prev;
        for ( auto i : ncrange(v.size()) ) {
          auto e = vectAt(v,i);
          if (dbg_prev.has_value()&&dbg_prev.value()!=0.0) {
            NCRYSTAL_MSG("sort and cleaned ["<<i<<"/"<<v.size()<<"]: "<<e<<" rel prev: "<<e/dbg_prev.value());
          } else {
            NCRYSTAL_MSG("sort and cleaned ["<<i<<"/"<<v.size()<<"]: "<<e);
          }
          dbg_prev = e;
        }

#ifndef NDEBUG
        for (std::size_t i = 0; i + 1 < v.size(); ++i) {
          nc_assert( v[i + 1] >= tol * v[i]
                     || floateq(v[i + 1], tol * v[i]) );
        }
#endif

        if ( tolmax > 0.0 ) {
          nc_assert_always( tolmax > tol*tol );
          VectD u;
          u.reserve(v.size()+10);
          // auto vv = Span<const double>(v).subspan(1);
          // u.push_back(v.front());
          double prev = v.front();
          for ( auto e : v ) {
            while ( prev> 0.0 && prev*tolmax < e ) {
              prev = prev*tolmax;
              u.push_back(prev);
              NCRYSTAL_MSG("TKTEST adding tolmax pt: "<<prev);
            }
            prev = e;
            u.push_back(e);
          }
          std::swap(u,v);
        }

      }



      void merge_into(VectD& v, const VectD& new_points,
                      double tol ) {
        auto is_dup = [tol](double a, double b)
        {
          return std::abs(a - b) <= 0.05 * std::max(std::abs(a), std::abs(b));
        };

        auto itv = v.begin();
        for (double x : new_points) {
          while (itv != v.end() && *itv < x)
            ++itv;
          bool near_existing = false;
          if (itv != v.end() && is_dup(x, *itv))
            near_existing = true;
          if (!near_existing && itv != v.begin() && is_dup(x, *std::prev(itv)))
            near_existing = true;
          if (!near_existing)
            itv = v.insert(itv, x);
        }
      }

      VectD
      trim_x_by_equal_yintegral( const double* x_begin, const double* x_end,
                                 const double* y_begin, const double* y_end,
                                 std::size_t N )
      {
        const std::size_t M = static_cast<std::size_t>(x_end - x_begin);
        nc_assert_always( M == static_cast<std::size_t>(y_end - y_begin) );
        nc_assert_always( N >= 2 );
        if ( M == 0 )
          return {};
        if (N >= M)
          return VectD(x_begin, x_end);

        // CDF (trapezoidal):
        VectD cdf;
        cdf.reserve(M-1);
        StableSumKahan area;
        for (std::size_t i = 1; i < M; ++i) {
          area.add( (y_begin[i] + y_begin[i - 1])
                    * (x_begin[i] - x_begin[i - 1]) );
          cdf.push_back(area.sum());
        }
        //area.mult( 0.5 );//trapezoidal
        if (!(area.sum()>0.0))
          return {};

        const double norm( 1.0 / area.sum() );
        for ( auto& e : cdf )
          e *= norm;
        nc_assert_always( cdf.size()+1 == M );

        for ( auto& e : cdf )
          NCRYSTAL_MSG("TKTEST cdf[..]="<<fmt(e));
        // Pick grid x at evenly spaced CDF levels
        VectD out;
        out.reserve(N);

        std::size_t k = 0;//track first index with cdf[k] >= t
        const double invNm1 = 1.0 / (N-1.0);
        for (std::size_t j = 0; j < N; ++j) {
          //advance to next quantile:
          double t = j * invNm1;
#if 1
          //Enhance the tails (try powers of 1,2,3,4,... or even between):
          //t = 1.0 - ncsquare(ncsquare(1-t));//ok... seems to be best
          t = 1.0 - ncsquare(1.0-t);//okish
          //t = 1.0 - ncsquare(ncsquare(ncsquare(1-t)))//no ok;
#endif
#if 1
          constexpr double oneminuseps = 1.0-1e-14;
          if ( t >= oneminuseps ) {
            const double lastx = *std::prev(x_end);
            if ( out.empty() || lastx > out.back() )
              out.push_back(lastx);
            break;
          }
#endif
          while (k+1 < M && vectAt(cdf,k) < t) {
            ++k;
          }

          if (k+1 >= M)
            break;// cannot satisfy further quantiles

          // if ( !( out.empty() || x_begin[k]>out.back() ) ) {
          //   NCRYSTAL_MSG( "FAIL x_begin[k]="
          //                 << fmt(x_begin[k])
          //                 <<" out.back()="
          //                 << fmt(out.back())
          //                 << " x_begin[k]-out.back()="<<fmt(x_begin[k]-out.back())
          //                 );
          // }
          //          nc_assert( out.empty() || x_begin[k]>out.back() );
          if ( out.empty() || x_begin[k]>out.back() )
            out.push_back(x_begin[k]);
          //          NCRYSTAL_MSG("  quantile t="<<t<<" gave beta="<<out.back());

          // Move k one step forward, so won't repeat same x.
          if (k + 2 < M)
            ++k;
        }

        return out;
      }

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
        constexpr double c3 = -139./51840.;
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


      VectD setupE0ABGrid( const VDOSGn& Gn, double msd, unsigned nbins )
      {
        //In the limit E=0, the available phasespace approaches a narrowing
        //region around the line alpha=beta. Based on the Gn functions, find
        //grid points (same for both alpha and beta grids) which ensures a good
        //granularity in the region near E=0:

        //Fixme: Also consider higher order contributions, until the point where the
        //Gn functions have no positive components?

        auto rawspec_G1 = Gn.getRawSpectrum(1);
        auto erange_G1 = Gn.eRange(1);
        //Convert E->beta:
        const double kT = Gn.kT();
        const double invkT = 1.0/kT;
        erange_G1.first *= invkT;
        erange_G1.second *= invkT;
        NCRYSTAL_MSG("TKTEST erange G1 (beta): "<<erange_G1.first<<" "<<erange_G1.second);

        // nc_assert(floateq(erange_G1.first,-G1));
        auto evals_G1 = linspace(erange_G1.first,erange_G1.second,rawspec_G1.size());
        nc_assert_always(evals_G1.front()<0.0);
        nc_assert_always(evals_G1.back()>0.0);


        //Fixme: put this scoped block into a helper function

        //For low energy we have a phasespace tight around the line at
        //alpha~=beta.  For G1, the contribution to this alpha shape is
        //proportional to exp(-2W) which is proportional to exp(-alpha*constant)
        //and near alpha=beta this becomes a weight of exp(-beta*constant). But
        //in addition to this, the phasespace also has a width proportional to
        //sqrt(beta) (in the limit E->0). So in addition to G1 itself, the
        //contribution at a given beta value will be weighted by
        //sqrt(beta)*exp(-2W), with 2W given be beta (or alpha) times the
        //constant:
        const double alpha2x = V2SKDetail::kTmsd_to_alpha2x(kT,msd);
        //find first non-negative beta value in evals_G1:
        nc_assert(evals_G1.size()==rawspec_G1.size());
        auto itB_begin = std::lower_bound(evals_G1.begin(), evals_G1.end(), 0.0);
        auto itB_end = evals_G1.end();
        auto itG1AtB_begin = std::next(rawspec_G1.begin(),
                                       itB_begin-evals_G1.begin());
        //Now, find the integrated contributions at E=0 (up to some normalisation
        //factor).
        VectD tmp;
        tmp.reserve( static_cast<std::size_t>(itB_end-itB_begin) );
        auto itG1AtB = itG1AtB_begin;
        auto itB = itB_begin;
        for ( ; itB!=itB_end; ++itB, ++itG1AtB ) {
          const double twoW = alpha2x*(*itB);
          tmp.push_back( std::sqrt(*itB)*twoW*std::exp(-twoW)*(*itG1AtB) );
          // NCRYSTAL_MSG("TKTEST beta = "<<*itB<<" contrib = "<<tmp.back());
        }
        nc_assert( itG1AtB == rawspec_G1.end() );
        //const std::size_t nE0 = 100;//fixme
        VectD res
          = V2SKDetail::trim_x_by_equal_yintegral( &*itB_begin, &*itB_end,
                                                   tmp.data(),
                                                   tmp.data() + tmp.size(), nbins );//fixme: something more tail oriented, not just N equal quantiles.
        // fixme... der er noget galt med enhederne eller noget.
        // res.push_back(4);
        // res.push_back(5);
        // res.push_back(8);
        // res.push_back(10);
        // res.push_back(12);
        // res.push_back(16);
        // res.push_back(20);
        // res.push_back(24);
        // res = linspace(0.0,120.0,200.0);
        return res;
      }

      void nc_array_add_inplace( double * ncrestrict tgt,
                                 const double * ncrestrict src,
                                 std::size_t n )
      {
        for ( std::size_t i = 0; i < n; ++i )
          tgt[i] += src[i];
      }

#if 0
      //new ai:

#else
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

        const unsigned maxOrder = Gn_asym.maxOrder().value();
        const double kT = Gn_asym.kT();
        nc_assert( kT > 0.0 );
        const auto nbeta = betaGrid.size();
        const auto nalpha = alphaGrid.size();
        VectD sab( nalpha * nbeta, 0.0);
        VectD beta_kT( betaGrid );
        // for ( auto e : enumerate(beta_kT) )
        //   NCRYSTAL_MSG("TKTEST beta["<<e.idx<<"] = "<<fmt(e.val));
        //0.38171954446335293


        // nc_assert( nc_is_grid( beta_kT ) );
        for ( auto& e : beta_kT )
          e *= kT;
        // for ( auto e : enumerate(beta_kT) )
        //   NCRYSTAL_MSG("TKTEST beta_kT["<<e.idx<<"] = "<<fmt(e.val));
        // NCRYSTAL_MSG("TKTEST kT="<<kT);
        nc_assert( nc_is_grid( beta_kT ) );
        VectD workbuf_gneval;
        VectD gneval;
        gneval.reserve(nbeta);


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
        VectD logx_vals = ( maxOrder >= stirling_threshold
                            ? vectorTrf(x_vals,[](double x){ return ( x > 0.0 ? std::log(x) : -kInfinity ); } )
                            : VectD() );
        auto fxn_cache = expmhalfx_vals;//We apply half of the exp(-x) factor before
        //the x^n/n! factor is added, and the other
        //half of exp(-x) we add later. This is done
        //to extend range of valid x-values.
        VectD alpha_factors;
        alpha_factors.resize(x_vals.size());

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
          //          auto betakT_Range = Gn_asym.eRange(n);
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
              nc_assert_always(exparg <= 708.0);
              vectAt(alpha_factors,x.idx) = fact * std::exp(exparg);
              //fixme: break if the alpha_factor just calculated is below
              //alpha_factor_floor and x>n, since the rest will just be zero??
            }
          }

          //Fixme: replace the following with binary search?
          constexpr double alpha_factor_floor = 1e-100;//fixme tune? Or just 0?
          std::size_t ialphaB = 0;
          std::size_t ialphaE = nalpha;
          while ( ialphaB != ialphaE
                  && vectAt(alpha_factors,ialphaB) <= alpha_factor_floor  )
            ++ialphaB;
          while ( ialphaB != ialphaE
                  && vectAt(alpha_factors,ialphaE-1) <= alpha_factor_floor  )
            --ialphaE;

            // for ( std::size_t ialpha = 0; ialpha < alpha_factors.size(); ++ialpha ) {
            // // for ( auto alpha_fact : enumerate(alpha_factors) ) {
            //   double alpha_fact = vectAt(alpha_factors,ialpha);



          // double* sabIt = sab.data();

          //Evaluate Gn function:
          nc_assert( nc_is_grid( beta_kT ) );
          Gn_asym.evalMany( n, beta_kT, gneval, workbuf_gneval );

          for ( std::size_t ibeta = 0; ibeta < nbeta; ++ibeta ) {
            //double beta = vectAt(betaGrid,ibeta);
            //for ( auto beta : enumerate(betaGrid) ) {
            //double energy = beta * kT;
            // if (!valueInInterval(betakT_Range,energy)) {
            //   // sabIt += nalpha;
            //   continue;//Gn(beta) zero here.
            // }

            //const double Gn_asym_eval = contribScaleFactor * Gn_asym.eval(n,energy);
            const double Gn_asym_eval = contribScaleFactor * vectAt(gneval,ibeta);//fixme skip earlier if contribScaleFactor is 0
            if ( !(Gn_asym_eval>0.0) ) {
              // sabIt += nalpha;
              continue;
            }
            //const auto offset_alpharow = ibeta*nalpha;
            //            const auto offset_alpharow_posbeta = posbeta_idx*nalpha;
            //double* sabrow = &sab[offset_alpharow];
            //const std::size_t nalpha = alpha_factors.size();
            const double * itAFact = alpha_factors.data() + ialphaB;
            const double * itAFactE = alpha_factors.data() + ialphaE;
            double* itSab = &sab[ibeta*nalpha+ialphaB];
            for ( ; itAFact != itAFactE; ++itAFact, ++itSab )
              *itSab += ( *itAFact * Gn_asym_eval );

            // for ( std::size_t ialpha = ialphaB; ialpha < ialphaE; ++ialpha ) {
            // // for ( auto alpha_fact : enumerate(alpha_factors) ) {
            //   double alpha_fact = vectAt(alpha_factors,ialpha);
            //   ///FIXME REMOVED AS TEST if ( !(alpha_fact>0.0) )
            //   ///FIXME REMOVED AS TEST   continue;
            //   double contrib_S_beta = alpha_fact * Gn_asym_eval;
            //   *sabIt++ += contrib_S_beta;
            //   //*sabrow++ += contrib_S_beta;
            //   //FIXME REMOVE COMMENT //Negative beta:
            //   // vectAt( sab, offset_alpharow + ialpha ) += contrib_S_beta;
            // }//alpha loop
          }//beta loop
        }//phonon order loop

        return sab;
      }
#endif
      VectD fillSABFromVDOSConcurrent( const VDOSGn& Gn_asym,
                                       const double msd,
                                       const VectD& alphaGrid,
                                       const VectD& betaGrid,
                                       ScaleGnContributionFct scaleGnContribFct )
      {
        //FIXME: NEW PLAN FOR CONCURRENCY WHICH DOES NOT NEED HUGE SAB TABLE FOR EACH ORDER: Split up the betagrid or alpha grid instead, and operate on spans of the final sab table. Probably best to split up along alpha, so we don't needlessly calculate alpha factors. Also, if alpha range has zero overlap, we can avoid gn eval entirely for a given order.
        //
        //
        //





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
        const unsigned njobs = ( norders <= 16
                                 ? 1
                                 : std::min<unsigned>(1,norders / 16) );//fixme: upper limit of... 6?
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
          VectD tmp;
          src.swap( tmp );//immediate dealloc
        }
        return res;
      }

    }
  }
}

NC::PairDD NC::rangeXNexpMX(unsigned n, double eps, double accuracy ) {
  //FIXME: cache high-res results for lowest n=200 orders?


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

NC::Optional<NC::PairDD>
NC::findExtremeSABPointWithinAlphaPlusCurve(double E_div_kT,
                                            PairDD alphaRange,
                                            PairDD betaRange)
{
  //Find the extreme (as in highest alpha, lowest beta) kinematically
  //accessible point in the provided rectangular region in (alpha,beta) space
  //for a neutron with energy/kT= Emax_div_kT. Returns NullOpt in case no point
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
    return NullOpt;//no accessible points in region
  }

  auto alphaPlus = [E_div_kT](double beta)
                   {
                     nc_assert( beta >= -E_div_kT );
                     return 2*E_div_kT + beta + 2 * std::sqrt( E_div_kT * ( E_div_kT + beta ) );
                   };
  const double apb1 = alphaPlus(betaRange.second);
  if ( apb1 <= alphaRange.first ) {
    nc_assert(!should_be_accessible);
    return NullOpt;//no accessible points in region
  }

  nc_assert(should_be_accessible);

  //Clip lower beta range at -E/kT:
  betaRange.first = ncmax( betaRange.first, -E_div_kT );

  const double apb0 = alphaPlus(betaRange.first);
  if ( apb0 >= alphaRange.second )
    return PairDD( alphaRange.second, betaRange.first );//entire rectangle is accessible

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
  return PairDD( alphaRange.second, betaRange.first );
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

#define NTOTAL_HACK_EXTRA 1//fixme
#define ORDER_LIMIT_EXTRA 10//fixme
NC::VectD NC::setupAlphaGrid( double kT, double msd,
                              const VectD& ptsE0,
                              double alphaMax, unsigned npts )
{
  nc_assert(npts>=20);
  npts *= NTOTAL_HACK_EXTRA;
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
  NCRYSTAL_MSG("cpp alpha2x = "<<alpha2x);
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

  auto finalise_grid = [&grid_region0,npts,&ptsE0](const VectD& grid) -> VectD
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
                         nc_assert_always(nc_is_grid(out));//fixme: always
                         nc_assert_always(out.size()==npts);

                         (void)V2SKDetail::merge_into;
#if 0 // fixme
                         auto v = linspace( 1e-10, 25.0, 200 );
                         V2SKDetail::merge_into( out, v, 0.05 );
                         for ( auto e: v )
                           NCRYSTAL_MSG(" alphagrid: "<<e);
                         nc_assert_always( nc_is_grid(out) );
#endif

                         if (!ptsE0.empty() ) {
                           if ( ptsE0.front() > 0.0 ) {
                             V2SKDetail::merge_into(out,ptsE0, 0.05);
                           } else {
                             auto ite = ptsE0.begin();
                             auto iteE = ptsE0.end();
                             while ( ite!=iteE && !(*ite>0.0) )
                               ++ite;
                             V2SKDetail::merge_into(out,VectD(ite,iteE), 0.05);
                           }
                         }
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

  //Reach of 1-phonon and 3-phonon spectrums towards negative beta (converted
  //E->beta):
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
#if 0
  const double epsilon = 0.1*Gn.binWidth(1)*invkT;//FIXME: I ADDED invkT now!!!!!!!
  const double mepsilon = -epsilon;
  nc_assert_always(evals_G1.front()<mepsilon);
  while ( evals_G1.back()>mepsilon )
    evals_G1.pop_back();
#else
  const double epsilon = -0.1*Gn.binWidth(1);
  nc_assert_always(evals_G1.front()<epsilon);
  while ( evals_G1.back()>epsilon )
    evals_G1.pop_back();
#endif

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
      std::tie(evals_G1,rawspec_G1) = reducePtsInDistribution( evals_G1,
                                                               rawspec_G1,
                                                               n1_spectrum );
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

namespace NCRYSTAL_NAMESPACE {
  namespace V2SKDetail {
    namespace {
      struct Spectrum {
        VectD beta;
        VectD weight;
      };
      VectD getMostImportantSpectrumPts( Spectrum&& in,
                                         std::size_t n, bool fill_gap )
      {
        VectD out;
        if ( in.beta.size() <= n ) {
          nc_assert_always( in.beta.size() == n );//fixme, we should have reduced n already if this happens
          out = std::move( in.beta );
          return out;
        }
        VectD& x = in.beta;
        VectD& y = in.weight;
        //Remove least important points from spectrum:
        if ( !fill_gap ) {
          //Original (NCrystal 2.0) approach, just run
          //reducePtsInDistribution. This might leave huge gaps if for instance
          //VDOS is zero over a wide internal region (fixme: revisit if this is
          //again good enough after moving to new knl processing):
          VectD tmp;
          std::tie(out,tmp) = reducePtsInDistribution( x, y, n );
          return out;
        };

        //Improved (NCrystal 2.6) approach, avoid leaving huge gaps without
        //pts. Doing so can leave large gaps in the beta points, which when can
        //lead to numerical artifacts when a scattering kernel is later
        //integrated and sampled.
        std::size_t n_for_gaps
          = ( ( n > 30 && x.size()-n > 10 )
              ? std::max<std::size_t>(5,static_cast<std::size_t>(n*0.1+0.5))
              : 0 );
        n -= n_for_gaps;
        const double binwidth = x.at(1) - x.at(0);
        std::tie(x,y) = reducePtsInDistribution( x,y, n );
        if ( n_for_gaps > 0 ) {
          struct Gap {
            Gap(double bbb0, double bbb1)
              : b0(bbb0), b1(bbb1)
            {}
            double b0, b1;//edges of gap
            std::size_t n = 0;//number of points allocated to fill the gap
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
          const double dmin = 1.5 * binwidth;
          auto it = x.begin();
          auto itLast = std::prev(x.end());
          for ( ; it != itLast; ++it ) {
            const double dd = *std::next(it)-*it;
            if ( dd > dmin )
              gaps.emplace_back( *it, *std::next(it) );
          }
          //Allocate pts one by one into largest remaining gap:
          while (n_for_gaps>0 && !gaps.empty() ) {
            std::stable_sort(gaps.begin(),gaps.end());
            gaps.front().n += 1;
            --n_for_gaps;
          }
          //Transfer points in gaps to x:
          for ( const auto& g : gaps ) {
            if ( g.n > 0 ) {
              double bw = (g.b1-g.b0)/(g.n+1.0);
              for ( auto i : ncrange(g.n) ) {
                x.push_back( g.b0 + (i+1)*bw );
              }
            }
          }
          std::sort(x.begin(),x.end());
        }
// #if 0 //fixme: need final "extra pt inserter"
//         if ( n_for_gaps > 0 ) {
//           //unused gap filler points, put somewhere else:
//           n1_near0 += n_for_gaps;
//         }
// #endif
        out = std::move(x);
        return out;
      }

      std::pair<Spectrum,Spectrum> extractSpectrum( const VDOSGn& Gn,
                                                    VDOSGn::Order order=1 )
      {
        //extracts Gn spectrum + splits in positive and negative parts.
        const double invkT = 1.0 / Gn.kT();
        const VectD& y = Gn.getRawSpectrum(order);
        VectD x;
        {
          auto xrange = Gn.eRange(order);
          //Convert E->beta:
          xrange.first *= invkT;
          xrange.second *= invkT;
          //nc_assert(floateq(xrange.first,-G1));
          x = linspace(xrange.first,xrange.second,y.size());
        }
        nc_assert_always(x.front()<0.0);
        nc_assert_always(x.back()>0.0);
        //Split into negative ([-xmin,0]) and positive parts ((0,xmax]):
        std::pair<Spectrum,Spectrum> res;
        auto& neg = res.first;
        auto& pos = res.second;
        const std::size_t n = x.size();
        neg.beta.reserve(n);
        neg.weight.reserve(n);
        pos.beta.reserve(n);
        pos.weight.reserve(n);
        //Value at 0.0 might actually have ended up slightly positive or
        //negative. We "snap" this value to 0.0 again, by comparison with the
        //binwidth:
        const double eps = 0.1*Gn.binWidth(1)*invkT;
        for ( std::size_t i = 0; i < n; ++i ) {
          if ( x[i] <= eps ) {
            neg.beta.push_back(x[i]);
            neg.weight.push_back(y[i]);
          } else {
            pos.beta.push_back(x[i]);
            pos.weight.push_back(y[i]);
          }
        }
        nc_assert_always( !neg.beta.empty() && !pos.beta.empty() );
        if ( ncabs(neg.beta.back()) < eps )
          neg.beta.back() = 0.0;//snap ~0 to 0
        // for ( auto& e : neg.beta )
        //   NCRYSTAL_MSG(" neg.beta: "<<e);
        // for ( auto e : pos.beta )
        //   NCRYSTAL_MSG(" pos.beta: "<<e);
        nc_assert_always( neg.beta.size() == neg.weight.size() );
        nc_assert_always( pos.beta.size() == pos.weight.size() );
        return res;
      }

    }
  }
}


std::pair<NC::VectD,NC::VectD>
NC::setupBetaGridNEWNEW( const NC::VDOSGn& Gn, double betaMax,
                         unsigned vdoslux, std::size_t override_nbins )
{
#if 1
  //testing
  double binwidth_prev = -1.0;
  for ( unsigned order = 1; order <= Gn.maxOrder().value(); ++order ) {
    const double bw = Gn.binWidth(order);
    if ( bw == binwidth_prev )
      continue;
    binwidth_prev = bw;
    NCRYSTAL_MSG("TKTEST G"<<order<<" binwidth: "<<fmt(bw));
  }
  NCRYSTAL_MSG("TKTEST final G"<<Gn.maxOrder().value()<<" binwidth: "<<fmt(Gn.binWidth(Gn.maxOrder())));
#endif



  //Returns { betaGrid, alphaPtsForLowE }
  nc_assert(Gn.maxOrder().value()>=1);
  nc_assert_always(vdoslux<=5);
  nc_assert_always(betaMax>0.0);

  //Reach of 1-phonon and 3-phonon spectrums towards negative beta (converted
  //E->beta):
  const double invkT = 1.0/Gn.kT();
  const double G1 = ncabs(Gn.eRange(1).first*invkT);
  NCRYSTAL_MSG(" ncabs(Gn.eRange(1).first*invkT) = "<<ncabs(Gn.eRange(1).first*invkT));
  NCRYSTAL_MSG(" ncabs(Gn.eRange(2).first*invkT) = "<<ncabs(Gn.eRange(2).first*invkT));
  NCRYSTAL_MSG(" ncabs(Gn.eRange(3).first*invkT) = "<<ncabs(Gn.eRange(3).first*invkT));
  nc_assert(G1>0.0);
  unsigned order_g3 = std::min<unsigned>(3,Gn.maxOrder().value());
  NCRYSTAL_MSG(" order_g3="<<order_g3);
  double G3 = ncabs(Gn.eRange( order_g3 ).first*invkT);//nb: actually =G1 or =G2 in case of maxOrder<3.
  NCRYSTAL_MSG(" G1="<<fmt(G1));
  NCRYSTAL_MSG(" G3="<<fmt(G3));
  nc_assert( Gn.maxOrder().value()==1 || G3>G1 );
  if (G1==G3)
    G3 = G1*1.0001;
  NCRYSTAL_MSG(" G1="<<fmt(G1));
  NCRYSTAL_MSG(" G3="<<fmt(G3));


  //Construct betagrid which will run from [-betaMax,D] with
  //0<D<=betaMax. Assuming betaMax is large enough, D will always be large
  //enough to encompass the first phonon order (preferably the first 3).
  //
  //We divide pts in regions:
  //Region0: Pts to ensure good coverage of the E->0 phasespace. These pts must
  //          also be duplicated in the alpha grid.
  //Region1: [-G1,0]: Pts to ensure description of the structures in G1 + a few
  //         points near beta=0 to enhance lowE.
  //Region2: [-D,-G1) (D is normally the reach of G3)
  //Region3: [-betaMax,-D).
  //
  //always encompass single-phonons, and clip G3/D to betaMax if needed, to
  //ensure no regions will have 0 range:
  betaMax = ncmax(betaMax,G1*1.01);
  G3 = ncmin(G3,betaMax*0.9999);
  //const double D = ncmin(betaMax*0.9999,ncmax(G3,std::max<int>(2,vdoslux-1)*10.0));
  const double D = ncmin(betaMax*0.9999,G3);//fixme: allow empty g3 region instead of this!!!!
  nc_assert( D>0.0 && G3 > G1 && D>=G3 && D<=betaMax );

  //How many points to use? The total is given as 100 * 2^(vdoslux), meaning
  //100 for vdoslux=0, 200 for vdoslux=1, and so on up to 3200 for vdoslux=5:
  //FIXME: Add 40 pts for vdoslux=0 if it make sense (it most likely does)!
  const std::size_t ntotal
    = ( override_nbins ? override_nbins
        : ( vdoslux==0
            ? 140 //100 for vdoslux=0 is simply not enough (fixme revisit)
            : 100*(1<<vdoslux) ) );//100 * 2^(vdoslux).


  //First do a bit of pre-analysis for the G1 spectrum in Region1:
  //Here we don't just space the points uniformly, but try to place them to
  //get the single-phonon curve best described. Additionally, we use 10% (but
  //at least 5) of the points to put some very small values around
  //beta=0. This is not needed to describe the actual data, but is used to
  //reduce artefacts in some integration/sampling algorithms (ideally those
  //algs should instead be updated to avoid such artefacts all by
  //themselves...).

  auto spectrum = V2SKDetail::extractSpectrum( Gn, 1 );
  auto& specNeg = spectrum.first;
  auto& specPos = spectrum.second;

  //Always allocate a flat 10% for lowE upscat:
  const std::size_t n0 = std::min<std::size_t>(ntotal / 10, specPos.beta.size() );

  //The absolute maximum number of points which can be useful in region 1 is
  //evals_G1.size() plus a bit for fine-grained grid near beta=0:
  //  const double[] nearb0pts_eV = { -1e-3, 1e-6, 1e-9, 1e-8 }
  //fixme: tinyvector:
  SmallVector<double,6> nearb0pts_in_reg1 = { 1e-12, 1e-10, 1e-8, 1e-6, 1e-4 };
  for ( auto& e : nearb0pts_in_reg1 )
    e *= invkT;
  while ( !nearb0pts_in_reg1.empty() && nearb0pts_in_reg1.back() >= G1*0.9 )
    nearb0pts_in_reg1.pop_back();
  for ( auto& e : nearb0pts_in_reg1 )
    e *= -1.0;
  std::reverse(nearb0pts_in_reg1.begin(),nearb0pts_in_reg1.end());
  nc_assert( nc_is_grid(nearb0pts_in_reg1) );

  const std::size_t n1_max = static_cast<std::size_t>( 0.5 + specNeg.beta.size()
                                                       + nearb0pts_in_reg1.size() );
  nc_assert_always(n1_max>=8);

  //Determine how many pts to assign to regions 1..3, depending on the
  //length of those regions:

  //Length of beta-scale covered by various regions:
  nc_assert_always(D>G1);
  nc_assert_always(betaMax>D);
  NCRYSTAL_MSG(" Raw G1="<<fmtg(G1));
  NCRYSTAL_MSG(" Raw D="<<fmtg(D));
  NCRYSTAL_MSG(" Raw betaMax="<<fmtg(betaMax));


  double L1 = G1;
  double L2 = (D-G1);
  double L3 = betaMax-D;
  NCRYSTAL_MSG(" Raw L1="<<fmtg(L1));
  NCRYSTAL_MSG(" Raw L2="<<fmtg(L2));
  NCRYSTAL_MSG(" Raw L3="<<fmtg(L3));
  //Ad-hoc boost of central parts where we expect more variation:
  L1 *= 4;
  L2 *= 2;
  NCRYSTAL_MSG(" Boosted L1="<<fmtg(L1));
  NCRYSTAL_MSG(" Boosted L2="<<fmtg(L2));
  NCRYSTAL_MSG(" Boosted L3="<<fmtg(L3));
  double Lnorm = 1.0/(L1+L2+L3);
  double f1 = ncmax(0.2,L1*Lnorm);//at least 20% in region 1
  double f2 = L2*Lnorm;
  NCRYSTAL_MSG(" f1="<<fmtg(f1));
  NCRYSTAL_MSG(" f2="<<fmtg(f2));

  //FIXME: Why?!?
  // if (f1+f2>0.99) {
  //   double tmp = 0.99/(f1+f2);
  //   f1 *= tmp;
  //   f2 *= tmp;
  // }

  //Before we construct the grid, we find the values needed for region0 and
  //region1. This allows us to reallocated any missing pts to region 2+3 if
  //somehow we end up using a bit fewer than expected (should happen only
  //rarely).

  //Region0:
  VectD pts0;
  {
    nc_assert_always( n0 >= 2 );
    const bool fill_gap = false;//fixme check if important
    pts0 = getMostImportantSpectrumPts( std::move(specPos), n0, fill_gap );
    nc_assert( nc_is_grid(pts0) );
    nc_assert_always( pts0.size() == n0 );
  }

  std::size_t n1 = std::min<std::size_t>(n1_max,static_cast<std::size_t>(f1*(ntotal-n0)+0.5));


  //Region1:
  VectD pts1;
  nc_assert( n1 > 0 );
  {
    const std::size_t n1_near0 = nearb0pts_in_reg1.size();
    nc_assert( n1 > n1_near0 );
    const std::size_t n1_spec = n1-n1_near0;
    nc_assert_always( n1_spec >= 2 );
    const bool fill_gap = true;//fixme check if important
    NCRYSTAL_MSG("getMostImportantSpectrumPts( specNeg[size="<<specNeg.beta.size()<<"], n1_spec="<<n1_spec);
    pts1 = getMostImportantSpectrumPts( std::move(specNeg), n1_spec, fill_gap );
    NCRYSTAL_MSG("  getMostImportantSpectrumPts... return npts = "<<pts1.size());
    nc_assert_always( pts1.size() == n1_spec );
    nc_assert( nc_is_grid(pts1) );
    if (!nearb0pts_in_reg1.empty()) {
      nc_assert( nc_is_grid(nearb0pts_in_reg1) );
#if 1
      VectD tmp;
      const auto nn = nearb0pts_in_reg1.size() + pts1.size();
      tmp.reserve( nn );
      auto it1 = pts1.begin();
      auto it2 = nearb0pts_in_reg1.begin();
      auto it1E = pts1.end();
      auto it2E = nearb0pts_in_reg1.end();
      auto addPt = [&tmp]( double val ) {
        if ( tmp.empty() || !floateq( val, tmp.back(), 0.01, 1e-13 ) )
          tmp.push_back(val);
      };

      while ( true ) {
        if ( it1 == it1E ) {
          for ( ; it2 != it2E; ++it2 )
            addPt(*it2);
          break;
        }
        if ( it2 == it2E ) {
          for ( ; it1 != it1E; ++it1 )
            addPt(*it1);
          break;
        }
        //both active:
        if ( *it1 < *it2 )
          addPt( *it1++ );
        else
          addPt( *it2++ );
      }
      nc_assert_always( tmp.size() <= nn );
      nc_assert_always( tmp.size() >= 2 );
      std::swap( tmp, pts1 );
#else
      V2SKDetail::merge_into(pts1, VectD( nearb0pts_in_reg1.begin(),
                                          nearb0pts_in_reg1.end() ), 0.01 );
#endif
    }
    nc_assert_always( pts1.size() <= n1 );
    nc_assert( nc_is_grid(pts1) );
  }
  nc_assert_always( pts1.size() <= n1 );
  if ( pts1.size() < n1 ) {
    //nleft += ( n1 - pts1.size() );
    n1 = pts1.size();
  }

  nc_assert_always( n0 == pts0.size() );
  nc_assert_always( n1 == pts1.size() );
  //Find n2 and n3:
  std::size_t n2 = ( L2 > 0.0
                     ? static_cast<std::size_t>((ncmax(0.15,(L2/(L2+L3)))*(ntotal-n0-n1)+0.5))
                     : 0 );                                               ;
  nc_assert_always(n1>=8);
  std::size_t n3 = ntotal - (n0+n1+n2);//fixme: don't add a few points by accident if the region3 (or2) is empty!!

  nc_assert_always( n0 + n1 + n2 + n3 == ntotal );
  NCRYSTAL_MSG(" TKTEST n0 = "<<n0);
  NCRYSTAL_MSG(" TKTEST n1 = "<<n1);
  NCRYSTAL_MSG(" TKTEST n2 = "<<n2);
  NCRYSTAL_MSG(" TKTEST n3 = "<<n3);
  NCRYSTAL_MSG(" TKTEST ntotal = "<<ntotal);

  nc_assert(n0<32768);
  nc_assert(n1<32768);
  nc_assert(n2<32768);
  nc_assert(n3<32768);
  nc_assert(n1>=8);
  nc_assert( ntotal == n0 + n1 + n2 + n3 );
  nc_assert_always(n0 <= ntotal );
  nc_assert_always(n1 <= ntotal );
  nc_assert_always(n2 <= ntotal );
  nc_assert_always(n3 <= ntotal );
  nc_assert_always(D>G1);//fixme: allow equal
  nc_assert_always(betaMax>D);//fixme: allow equal

  //Time to add values to the grid!

  //Region3:
  VectD grid;
  grid.reserve(ntotal);
  if ( n3 > 0 ) {
    auto vals = linspace(-betaMax, -D, n3+1);
    grid.insert(grid.begin(),vals.begin(),std::prev(vals.end()));
    nc_assert(grid.size()==n3);
    nc_assert_always( grid.front()==-betaMax );
  }
  nc_assert( grid.size()<2 || nc_is_grid(grid) );

  //Region2:
  if ( n2 > 0 ) {
    auto vals = linspace(-D, -G1, n2+1);
    nc_assert( grid.empty() || vals.front() > grid.back() );
    grid.insert(grid.end(),vals.begin(),std::prev(vals.end()));
    nc_assert(grid.size()==n2+n3);
  }
  nc_assert( grid.size()<2 || nc_is_grid(grid) );
  //Region 1:
  nc_assert( n1 > 0 );
  {
      // for ( auto& e : grid )
      //   NCRYSTAL_MSG(" grid before reg1: "<<fmt(e));
      // for ( auto& e : pts1 )
      //   NCRYSTAL_MSG(" pts1: "<<fmt(e));
    if ( !( grid.empty() || pts1.front() > grid.back() ) ) {
      // for ( auto& e : grid )
      //   NCRYSTAL_MSG(" grid before reg1: "<<fmt(e));
      // for ( auto& e : pts1 )
      //   NCRYSTAL_MSG(" pts1: "<<fmt(e));
    }

    nc_assert( pts1.back() == 0.0 );
    nc_assert( grid.empty() || pts1.front() > grid.back() );
    grid.insert(grid.end(),pts1.begin(),pts1.end());
    nc_assert(grid.size()==n1+n2+n3);
  }
  nc_assert( nc_is_grid(grid) );
  //Region 0:
  nc_assert( n0 > 0 );
  {
    nc_assert( grid.empty() || pts0.front() > grid.back() );
    grid.insert(grid.end(),pts0.begin(),pts0.end());
  }
  nc_assert( nc_is_grid(grid) );
  nc_assert_always( grid.front()==-betaMax );
  nc_assert_always( pts0.empty() || grid.back() == pts0.back() );
  nc_assert_always( nc_is_grid(grid) );
  if ( !( grid.size() == ntotal ) ) {
    NCRYSTAL_MSG("TKTEST grid.size() = "<<grid.size()<<" ntotal = "<<ntotal);
  }
  nc_assert_always( grid.size() == ntotal );
  return { grid, pts0 };
}

NC::VectD NC::setupBetaGridNEW( const NC::VDOSGn& Gn,
                                double msd,
                                double betaMax,
                                unsigned vdoslux,
                                unsigned override_nbins )
{
  //Distribute pts into the following categories:
  //
  //  CatLowE: Needed to ensure good behaviour as E->0. These points are
  //           normally in the interval [0,G1up] and are based on contributions
  //           to the phasespace at E->0, as well as a few points near 0 (and 0
  //           itself).
  //    CatG1: Needed to describe the peaks and features in the G1 region. These
  //           points are normally in the interval [G1lower,G1up], where
  //           G1lower<0. We additionally try to avoid large empty regions
  //           (often seen e.g. in hydrogenic materials)
  //   CatG23: Needed to capture any lingering features in G2 and G3, but only
  //           outside CatG1. These points are normally in the interval
  //           [G3lower,G1low] where both are negative.
  // CatHighE: Needed to cover region [-betaMax,G3lower].
  //
  // After first finding points in CatLowE and CatG1, these are merged into one
  // symmetric set of points around beta=0. Then the leftover points are added
  // to CatG23 and CatHighE. For CatHighE we use a powspace instead of a
  // geomspace, since geomspace leaves too large gaps at the upper edges (fixme:
  // test this).
  //
  // During all this, we also make sure that points are always at least 10%
  // apart, to avoid numerical artifacts due to loss of effective floating point
  // precision in the relative positions between two points. We also try to
  // avoid them being more than a factor of 10 apart, for similar reasons
  // (fixme??).

  //fixme: unsigned -> std::size all over?

  const double kT = Gn.kT();
  nc_assert(kT>0.0);
  const double invkT = 1.0/kT;
  nc_assert(invkT>0.0);
  nc_assert(Gn.maxOrder().value()>=1);
  nc_assert_always(vdoslux<=5);
  nc_assert_always(betaMax>0.0);

  //VectD ptsE0 = V2SKDetail::setupE0ABGrid( Gn_asym, msd, 20+3*vdoslux*vdoslux );//fixme: hardwire N?!?
  const unsigned ntotal = NTOTAL_HACK_EXTRA*[override_nbins,vdoslux]()
  {
    if ( override_nbins )
      return override_nbins;
    //100 * 2^(vdoslux), but increased 40% for vdoslux=0 since it was hard to
    //get a ~convergent result with just a 100x50 grid.
    switch( vdoslux ) {
    case 0: return 140u;
    case 1: return 200u;
    case 2: return 400u;
    case 3: return 800u;
    case 4: return 1600u;
    case 5:
    default:
      nc_assert(vdoslux==5);
      return 3200u;
    };
  }();

  //Find relevant cross over points (as values > 0):
  const double B1 = ncmax( ncabs(Gn.eRange(1).first),
                           ncabs(Gn.eRange(1).second) ) * invkT;
  const double B2 = [B1,&Gn,invkT]()
  {
    if ( Gn.maxOrder().value() >= 3 )
      return ncabs(Gn.eRange(3).first) * invkT;
    if ( Gn.maxOrder().value() == 1 )
      return B1;
    nc_assert( Gn.maxOrder().value() == 2 );
    return ncabs(Gn.eRange(2).first) * invkT;
  }();
  nc_assert_always(B2>=B1);
  // if ( !(betaMax>B2) ) {
  //   B2 *= 0.999;
  // }
  nc_assert_always(betaMax>=B2);

  //Length of beta-scale covered by various regions:
  double L1 = 2*B1;
  double L2 = B2-B1;//might be 0
  double L3 = betaMax-B2;//might be 0

  //Ad-hoc boost of more important central parts (fixme revisit):

  L1 *= 8;
  L2 *= 2;

  //Use this to figure out fractions of pts we want to pour into different
  //regions.

  double Lnorm = 1.0/(L1+L2+L3);
  double f1 = ncmax(0.2,L1*Lnorm);
  //double f2 = L2*Lnorm;
  // if (f1+f2>0.99) {
  //   nc_assert_false(0&&"fixme i forgot what happens here");
  //   double tmp = 0.99/(f1+f2);
  //   f1 *= tmp;
  //   f2 *= tmp;
  // }

  unsigned n1 = static_cast<unsigned>(f1*ntotal+0.5) / 2;//div 2 due to later
                                                         //duplication of pts
  //ensure n1 even and at least 14:
  n1 = std::max<unsigned>(14,(n1 % 2) ? (n1 - 1) : n1);

  NCRYSTAL_MSG("TKTEST n1 = "<<n1);

  VectD CatLowEG1_pospts_data;
  {
    //Figure out actual pts in n1 region. We aim for 40% (but at least 8) for the
    //low-E part:
    VectD ptsE0 = V2SKDetail::setupE0ABGrid( Gn, msd,
                                             std::max<unsigned>(8,n1*0.4+0.5) );
    NCRYSTAL_MSG("TKTEST ptsE0: {"<<ptsE0.front()<<" .. "<<ptsE0.back()<<", n="<<ptsE0.size());
    //fixme: if we also return ptsE0, it could be used be setupAlphaGridNEW ?
    //Then we can also input the rawspec etc. to this function?

    auto rawspec_G1 = Gn.getRawSpectrum(1);
    auto erange_G1 = Gn.eRange(1);
    erange_G1.first *= invkT;//E->beta
    erange_G1.second *= invkT;//E->beta
    auto evals_G1 = linspace( erange_G1.first,
                              erange_G1.second,
                              rawspec_G1.size() );
    NCRYSTAL_MSG("TKTEST evals_G1: {"<<evals_G1.front()<<" .. "<<evals_G1.back()<<", n="<<evals_G1.size());
    nc_assert_always(evals_G1.front()<0.0);
    nc_assert_always(evals_G1.back()>0.0);

    //Ignore positive part and 0 (NB: value representing zero can actually be
    //slightly non-zero due to numerical issues, hence need for epsilon):
    const double epsilon = 0.1*Gn.binWidth(1)*invkT;//FIXME: I ADDED invkT now!!
    const double mepsilon = -epsilon;
    nc_assert_always(evals_G1.front()<mepsilon);
    while ( evals_G1.back()>mepsilon )
      evals_G1.pop_back();
    nc_assert_always(evals_G1.size()>=2);
    rawspec_G1.resize(evals_G1.size());
    NCRYSTAL_MSG("TKTEST evals_G1[neg]: {"<<evals_G1.front()<<" .. "<<evals_G1.back()<<", n="<<evals_G1.size());

    //Remove least important points from G1 spectrum:
    unsigned n1_spectrum = n1 - ptsE0.size();
    nc_assert_always( n1_spectrum >= 6 );
#if 0
    //Old way, just run reducePtsInDistribution. This might leave huge gaps if
    //for instance VDOS is zero over a wide internal region:
    std::tie(evals_G1,rawspec_G1) = reducePtsInDistribution( evals_G1,rawspec_G1,
                                                             n1_spectrum );
#else
    //New way, avoid leaving huge gaps without pts. Doing so can leave large
    //gaps in the beta points, which when can lead to numerical artifacts when
    //a scattering kernel is later integrated and sampled. FIXME revisit.
    unsigned n_for_gaps = ( ( n1_spectrum > 30 && evals_G1.size()-n1_spectrum > 10 )
                            ? std::max<unsigned>(5,static_cast<unsigned>(n1_spectrum*0.1+0.5))
                            : 0 );
    n1_spectrum -= n_for_gaps;
    const double G1binwidth = evals_G1.at(1)-evals_G1.at(0);
    std::tie(evals_G1,rawspec_G1) = reducePtsInDistribution( evals_G1,rawspec_G1, n1_spectrum );
    NCRYSTAL_MSG("TKTEST evals_G1[reduced]: {"<<evals_G1.front()<<" .. "<<evals_G1.back()<<", n="<<evals_G1.size());

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
      gaps.reserve(evals_G1.size()-1);
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
    }
#endif
    //Now add in the ptsE0 pts, after first flipping:
    for ( auto& e : evals_G1 ) {
      nc_assert(e<0.0);
      e *= -1.0;
    }
    evals_G1.reserve(evals_G1.size()+ptsE0.size());
    for ( auto& e : ptsE0 )
      evals_G1.push_back(e);//fixme vectorAppend
    CatLowEG1_pospts_data = std::move(evals_G1);

    // double firstNonZero = [&CatLowEG1_pospts_data]()
    // {
    //   for ( auto e : CatLowEG1_pospts_data )
    //     if ( e > 0.0 )
    //       return e;
    //   nc_assert_always(false);
    //   return 1.0;
    // }();
    // vectorAppend(CatLowEG1_pospts_data,
    //              powspace(firstNonZero*1e-4,firstNonZero,7,2.0));
    CatLowEG1_pospts_data.push_back(1e-20);
    CatLowEG1_pospts_data.push_back(1e-10);
    const double tolmax = 2.5 - 0.26*vdoslux;
    //const double tolmax = 2.0 - 0.16*vdoslux;
    V2SKDetail::sortAndCleanupGrid(CatLowEG1_pospts_data, 1.02,tolmax);
    NCRYSTAL_MSG("TKTEST CatLowEG1_pospts_data[gapped,cleaned,sorted]: {"<<CatLowEG1_pospts_data.front()<<" .. "<<CatLowEG1_pospts_data.back()<<", n="<<CatLowEG1_pospts_data.size());
  }
  nc_assert( nc_is_grid(CatLowEG1_pospts_data) );
  Span<double> CatLowEG1_pospts(CatLowEG1_pospts_data);
  if ( CatLowEG1_pospts.front() == 0.0 )
    CatLowEG1_pospts = CatLowEG1_pospts.subspan(1);
  // for ( auto& e : CatLowEG1_pospts )
  //   NCRYSTAL_MSG("TKTEST CatLowEG1_pospts final "<<e);

  nc_assert_always( CatLowEG1_pospts.front() > 0.0 );

  //Now we can divide pts between L2 and L3:
  nc_assert( ntotal > ( CatLowEG1_pospts.size()*2+1) );
  auto nleft = ntotal-( CatLowEG1_pospts.size()*2+1);

  unsigned n2, n3;
  if ( !(L2>0.0) ) {
    if ( !(L3>0.0) )
      NCRYSTAL_THROW(BadInput,"Beta grid for cases with only G1 is not"
                     "yet implemented");
    n2 = 0;
    n3 = nleft;
  } else if ( !(L3>0.0) ) {
    n2 = nleft;
    n3 = 0;
  } else {
    n2 = static_cast<unsigned>(nleft*(L2/(L2+L3))+0.5);
    n3 = nleft-n2;
  }
  NCRYSTAL_MSG("TKTEST n2="<<n2<<" n3="<<n3);
  nc_assert_always(ntotal<10000);
  nc_assert_always(n2<10000);
  nc_assert_always(n3<10000);

  //Putting it all together, first n3 pts from [-betamax,-B2), using a powscale
  //(geomscale was too harsh.. FIXME ORIGINAL ACTUALLY USED LINSPACE):
  VectD res;
  res.reserve(ntotal);
  if ( n3 > 0 ) {
    //auto p = powspace(B2,betaMax,n3+1, 2.0 );//fixme: tune 2.0
    //fixme: linspace, so no need for reverse iter
    auto p = linspace(B2,betaMax,n3+1 );
    for (auto it = p.crbegin(); it != std::prev(p.crend()); ++it)
      res.push_back( - *it );
  }
  NCRYSTAL_MSG("TKTEST after n3: res.size() : "<<res.size()<<" (expected n3="<<n3<<")");
  NCRYSTAL_MSG("TKTEST after n3: res.back() : "<<res.back());
  nc_assert_always(  res.size() == n3 );
  nc_assert_always( nc_is_grid(res) );//fixme _always

  //Next n2 pts from [-B2,-B1):
  if ( n2 > 0 )  {
    auto v = linspace(-B2,-B1,n2+1);
    v.pop_back();
    vectorAppend(res,v);
  }
  nc_assert_always( nc_is_grid(res) );//fixme _always
  nc_assert_always( res.size() == n2 + n3 );//fixme _always

  //Now we replay CatLowEG1_pospts symmetrically around 0:
  for (auto it = CatLowEG1_pospts.rbegin();
       it != CatLowEG1_pospts.rend(); ++it)
    res.push_back( - *it );
  nc_assert(res.back()<0.0);
  res.push_back(0.0);
  for ( auto e : CatLowEG1_pospts )
    res.push_back(e);
  nc_assert_always(res.size()==ntotal);
  nc_assert_always( nc_is_grid(res) );//fixme _always
  nc_assert( res.front()==-betaMax );
  nc_assert( floateq(res.back(),B1) );

  // for ( auto e : res )
  //   NCRYSTAL_MSG("TKTEST final beta : "<<e);

  return res;
}

//fixme: place in VDOS namespace:

NC::GnExpansion NC::expandVDOSToGnFcts( const VDOSData& vdosdata,
                                        unsigned vdoslux,
                                        double targetEmax_requested,
                                        VDOSGn::TruncAndThinningParams ttpars,
                                        //ScaleGnContributionFct scaleGnContributionFct,
                                        Optional<unsigned> call_override_max_order )
{
  //Hidden unofficial env-vars used for special debugging purposes:
  const unsigned override_max_order = ( call_override_max_order.has_value()
                                        ? call_override_max_order.value()
                                        : static_cast<unsigned>(ncgetenv_int("HACK_MAXORDER")) );
  const double override_alphamax = ncgetenv_dbl("HACK_ALPHAMAX");
  const double override_betamax = ncgetenv_dbl("HACK_BETAMAX");
  // const unsigned override_nbins = ncgetenv_int("HACK_NBINS");

  //Which Emax should we target (i.e. aim to cover the kinematic reachable area
  //for neutrons of that energy):
  nc_assert_always( vdoslux <= 5 );
  nc_assert_always(targetEmax_requested>=0.0);
  //FIXME: playing around
  //  constexpr double lux2emax[6] = { 0.5, 1.0, 3.0, 5.0, 8.0, 12.0 };//Emax in eV for vdosluxs 0 to 5
  //  constexpr double lux2emax[6] = { 0.5, 1.0, 3.0, 5.0, 8.0, 12.0 };//Emax in eV for vdosluxs 0 to 5
  //constexpr double lux2emax[6] = { 2.0, 3.0, 4.0, 5.0, 5.0, 5.0 };//Emax in eV for vdosluxs 0 to 5
  constexpr double lux2emax[6] = { 5.0, 5.0, 5.0, 5.0, 5.0, 5.0 };//Emax in eV for vdosluxs 0 to 5
  //constexpr double lux2emax[6] = { 5.0, 5.0, 5.0, 5.0, 5.0, 5.0 };//Emax in eV for vdosluxs 0 to 5
  nc_assert_always( vdoslux < 6 );
  double targetEmax = targetEmax_requested>0.0 ? targetEmax_requested : lux2emax[vdoslux];

  if ( V2SKDetail::s_verbose )
    NCRYSTAL_MSG("VDOS2SK initialising with T="<<vdosdata.temperature()//fixme: VDOS2SK -> ?? VDOS2Gn ??
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
#if 0
  //fixme: old was always at least 4
  unsigned max_phonon_order = std::max<unsigned>(override_max_order,4);
#else
  //fixme: new, but unverified
  unsigned max_phonon_order = std::max<unsigned>(override_max_order,1);
#endif

  GnExpansion res{ VDOSGn(vdoseval,ttpars), 0.0, 0.0,
                   V2SKDetail::kTmsd_to_alpha2x(kT,msd), msd };
  auto& Gn_asym = res.Gn;

  Gn_asym.growMaxOrder(max_phonon_order);
  // NCRYSTAL_MSG("TKTEST MAX ORDER AFTER INITIAL GROW: "<<Gn_asym.maxOrder().value());

  //What are the highest phonon order we allow? When user requested given target
  //Emax, allow higher order expansions to accommodate. Otherwise, keep low to
  //keep default initialisation times reasonable (with special extreme settings
  //for vdoslux 0 and 5).
  unsigned order_limit = 1000;
  if ( targetEmax_requested > 0.0 || vdoslux == 5 )
    order_limit *= 10;
#if 0//fixme document if we change
  //old:
  if ( vdoslux==0 )
    order_limit /= 10;
#else
  if ( vdoslux==0 )
    order_limit /= 2;
#endif

  order_limit *= ORDER_LIMIT_EXTRA;

  const double emax_lowest_allowed = ( targetEmax_requested>0.0 ? targetEmax_requested : 1e-15 );

  //Now increase order dynamically until the last order only has contributions
  //to S(alpha,beta) outside the kinematic reach of Emax:
  //const double relcontriblvl = std::pow(10.0,-(3.0+2.0*vdoslux));//e.g.: 1e-3 for vdoslux 0, 1e-9 for vdoslux 3, 1e-13 for vdoslux 5
  const double relcontriblvl = 1e-13;//fixme
  const double x2alpha = 1.0 / res.alpha2x;
  auto findAlphaBetaRangeOfOrder
    = [&Gn_asym,x2alpha,invkT,relcontriblvl](unsigned n)
    {
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
      // NCRYSTAL_MSG("TKTEST Order was outside targetEmax: n="<<Gn_asym.maxOrder().value());
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

  // NCRYSTAL_MSG("TKTEST MAX ORDER AFTER EXPANSION LOOP: "<<Gn_asym.maxOrder().value());

  //Ok, we now know how many orders we need to reach targetEmax. Next step is to
  //look at the contribution of each order insided the kinematic reach of
  //targetEmax, and use it to determine alpha/beta limits:

  res.suggestedEmax = ( override_max_order>0 ? 0.0 : targetEmax);

  double betaMin = 0.0;
  double alphaMax = 0.0;
  for ( unsigned n = 1; n<=max_phonon_order; ++n ) {
    PairDD alphaRange, betaRange;
    std::tie(alphaRange, betaRange) = findAlphaBetaRangeOfOrder(n);
    auto ep = findExtremeSABPointWithinAlphaPlusCurve(targetEmax_div_kT, alphaRange, betaRange);
    if ( ep.has_value() ) {
      alphaMax = ncmax(alphaMax,ep.value().first);
      betaMin = ncmin(betaMin,ep.value().second);
      // NCRYSTAL_MSG("TKTEST after n="<<n<<" betaMin is "<<betaMin);
    }
  }
  nc_assert_always(betaMin<0.0 && alphaMax > 0.0);

  // if ( res.suggestedEmax != 0.0 ) {
  //   //fixme: this is new... and perhaps not needed?!?!??
  //   NCRYSTAL_MSG("TKTEST MOVING UP BETAMIN from "<<betaMin);
  //   betaMin = ncmax( -res.suggestedEmax*1.01*invkT, betaMin );
  //   NCRYSTAL_MSG("TKTEST MOVING UP BETAMIN to "<<betaMin);
  // }

#if 0//fixme we used to have:
  res.upper_beta = -betaMin*1.01;
  res.upper_alpha = alphaMax*1.01;
#else
  res.upper_beta = -betaMin;
  res.upper_alpha = alphaMax;
#endif
  if (override_alphamax)
    res.upper_alpha = override_alphamax;
  if (override_betamax)
    res.upper_beta = override_betamax;
  nc_assert_always( res.upper_beta>0.0 && res.upper_alpha>0.0 );

  // NCRYSTAL_MSG("TKTEST MAX ORDER ON RETURN: "<<res.Gn.maxOrder().value());
  return res;//{ std::move(Gn_asym), upper_alpha, upper_beta };


  // res.upper_alpha = upper_alpha;
  // res.upper_beta = upper_beta;
  // return res;
}


NC::ScatKnlData NC::createScatteringKernel( const VDOSData& vdosdata,
                                            unsigned vdoslux,
                                            double targetEmax_requested,
                                            VDOSGn::TruncAndThinningParams ttpars,
                                            ScaleGnContributionFct scaleGnContributionFct,
                                            Optional<unsigned> call_override_max_order )
{
  const auto gnexpn = expandVDOSToGnFcts( vdosdata, vdoslux, targetEmax_requested,
                                          ttpars, call_override_max_order );
  const auto& Gn_asym = gnexpn.Gn;

///  //Hidden unofficial env-vars used for special debugging purposes:
///  const unsigned override_max_order = ( call_override_max_order.has_value()
///                                        ? call_override_max_order.value()
///                                        : static_cast<unsigned>(ncgetenv_int("HACK_MAXORDER")) );
///  const double override_alphamax = ncgetenv_dbl("HACK_ALPHAMAX");
///  const double override_betamax = ncgetenv_dbl("HACK_BETAMAX");
///  const unsigned override_nbins = ncgetenv_int("HACK_NBINS");
///
///  //Which Emax should we target (i.e. aim to cover the kinematic reachable area
///  //for neutrons of that energy):
///  nc_assert_always( vdoslux <= 5 );
///  nc_assert_always(targetEmax_requested>=0.0);
///  //FIXME: playing around
///  //constexpr double lux2emax[6] = { 0.5, 1.0, 3.0, 5.0, 8.0, 12.0 };//Emax in eV for vdosluxs 0 to 5
///  //  constexpr double lux2emax[6] = { 5.0, 5.0, 5.0, 5.0, 5.0, 5.0 };//Emax in eV for vdosluxs 0 to 5
///  constexpr double lux2emax[6] = { 5.0, 5.0, 5.0, 5.0, 5.0, 5.0 };//Emax in eV for vdosluxs 0 to 5
///  nc_assert_always( vdoslux < 6 );
///  double targetEmax = targetEmax_requested>0.0 ? targetEmax_requested : lux2emax[vdoslux];
///
///  if ( V2SKDetail::s_verbose )
///    NCRYSTAL_MSG("VDOS2SK initialising with T="<<vdosdata.temperature()
///                 <<", vdoslux="<<vdoslux
///                 <<", aiming for Emax="<<targetEmax<<"eV"
///                 <<(targetEmax_requested>0.0?" (as requested)":"")<<", ...");
///
///  //Initialise evaluators:
///  VDOSEval vdoseval(vdosdata);
///  const double kT = vdoseval.kT();
///  const double invkT = 1.0/kT;
///  const double gamma0 = vdoseval.calcGamma0();
///  const double msd = vdoseval.getMSD( gamma0 );
///  double targetEmax_div_kT = targetEmax*invkT;
///  unsigned max_phonon_order = std::max<unsigned>(override_max_order,4);
///  VDOSGn Gn_asym(vdoseval,ttpars);
///  double upper_beta, upper_alpha;
///
///  {
///  Gn_asym.growMaxOrder(max_phonon_order);
///
///  //What are the highest phonon order we allow? When user requested given target
///  //Emax, allow higher order expansions to accommodate. Otherwise, keep low to
///  //keep default initialisation times reasonable (with special extreme settings
///  //for vdoslux 0 and 5).
///  unsigned order_limit = 1000;
///  if ( targetEmax_requested > 0.0 || vdoslux == 5 )
///    order_limit *= 10;
///#if 0//fixme document if we change
///  //old:
///  if ( vdoslux==0 )
///    order_limit /= 10;
///#else
///  if ( vdoslux==0 )
///    order_limit /= 2;
///#endif
///
///  order_limit *= ORDER_LIMIT_EXTRA;
///
///  const double emax_lowest_allowed = ( targetEmax_requested>0.0 ? targetEmax_requested : 1e-15 );
///
///  //Now increase order dynamically until the last order only has contributions
///  //to S(alpha,beta) outside the kinematic reach of Emax:
///  const double relcontriblvl = std::pow(10.0,-(3.0+2.0*vdoslux));//e.g.: 1e-3 for vdoslux 0, 1e-9 for vdoslux 3, 1e-13 for vdoslux 5
///  //const double relcontriblvl = 1e-13;//fixme
///  const double x2alpha = 1.0 / V2SKDetail::kTmsd_to_alpha2x(kT,msd);
///  auto findAlphaBetaRangeOfOrder = [&Gn_asym,x2alpha,invkT,relcontriblvl](unsigned n) {
///                                     auto eRange = Gn_asym.eRange(n, relcontriblvl);
///                                     PairDD betaRange( eRange.first * invkT, eRange.second * invkT  );
///                                     auto xRange = rangeXNexpMX( n, relcontriblvl );
///                                     PairDD alphaRange( xRange.first * x2alpha, xRange.second * x2alpha  );
///                                     return std::make_pair(alphaRange,betaRange);
///                                   };
///  while (true) {
///    if (override_max_order>0)
///      break;
///    Gn_asym.growMaxOrder(max_phonon_order);
///    PairDD alphaRange, betaRange;
///    std::tie(alphaRange, betaRange) = findAlphaBetaRangeOfOrder(Gn_asym.maxOrder().value());
///    if (sabPointWithinAlphaPlusCurve(targetEmax_div_kT,alphaRange.first,betaRange.second)) {
///      ++max_phonon_order;//could consider larger stepsize, but need to carefully check usage in the following
///    } else {
///      break;
///    }
///    if (max_phonon_order>order_limit) {
///      //Too slow - unfeasible to fill out S(alpha,beta) all the way out to the
///      //kinematic curve for E=targetEmax. In this case it is better to reduce
///      //targetEmax, to at least get a consistent table (and hope the free-gas
///      //extrapolation mechanisms will be adequate already at this lower
///      //threshold).
///      double targetEmax_reduced  = targetEmax;
///      do {
///        targetEmax_reduced *= 0.99;
///        if ( targetEmax_reduced < emax_lowest_allowed )
///          NCRYSTAL_THROW2(CalcError,"VDOS expansion too slow - can not reach E="<<emax_lowest_allowed
///                          <<"eV after "<<order_limit<<" phonon convolutions (likely causes: either the target energy"
///                          " value is too high, vdoslux too low, the temperature too high, or the VDOS is very unusual).");
///      } while (sabPointWithinAlphaPlusCurve(targetEmax_reduced*invkT,alphaRange.first,betaRange.second));
///      if (V2SKDetail::s_verbose)
///        NCRYSTAL_WARN("VDOS2SK Could only reach Emax="<<targetEmax_reduced<<"eV and not the requested Emax="<<targetEmax<<"K");
///      targetEmax_div_kT = targetEmax_reduced * invkT;
///      targetEmax = targetEmax_reduced;
///      break;
///    }
///  }
///  nc_assert_always( targetEmax_requested==0.0 || targetEmax_requested == targetEmax );
///  Gn_asym.growMaxOrder(max_phonon_order);
///
///  //Ok, we now know how many orders we need to reach targetEmax. Next step is to
///  //look at the contribution of each order insided the kinematic reach of
///  //targetEmax, and use it to determine alpha/beta limits:
///
///  double betaMin = 0.0;
///  double alphaMax = 0.0;
///  for ( unsigned n = 1; n<=max_phonon_order; ++n ) {
///    PairDD alphaRange, betaRange;
///    std::tie(alphaRange, betaRange) = findAlphaBetaRangeOfOrder(n);
///    auto ep = findExtremeSABPointWithinAlphaPlusCurve(targetEmax_div_kT, alphaRange, betaRange);
///    alphaMax = ncmax(alphaMax,ep.first);
///    betaMin = ncmin(betaMin,ep.second);
///  }
///  nc_assert_always(betaMin<0.0 && alphaMax > 0.0);
///  upper_beta = -betaMin*1.01;
///  upper_alpha = alphaMax*1.01;
///  if (override_alphamax)
///    upper_alpha = override_alphamax;
///  if (override_betamax)
///    upper_beta = override_betamax;
///  nc_assert_always( upper_beta>0.0 && upper_alpha>0.0 );
///  }
///

  //Ok, time to setup the alpha/beta grids. The grid-spacing is not even, rather
  //it attempts to best accomodate features of the distributions:

  // const unsigned override_nbins = ncgetenv_int("HACK_NBINS");

  unsigned nalpha, nbeta;
  std::tie(nalpha, nbeta) = VDOS::gridDimFromLux( vdoslux );
  VectD alphaGrid, betaGrid;
  std::tie( alphaGrid, betaGrid )
    = VDOS::determineAlphaBetaGridFromGn( gnexpn, nalpha, nbeta );
  nc_assert_always(nc_is_grid(alphaGrid));
  nc_assert_always(nc_is_grid(betaGrid));
  // for ( auto e : alphaGrid )
  //   NCRYSTAL_MSG("TKTEST alphaGrid "<<fmt(e));//fixme: some of the leading alphagrid bins are TOO CLOSE!
  nc_assert_always(alphaGrid.front()==0.0);


  // {
  //   std::pair<unsigned,unsigned> gridDimFromLux( unsigned vdoslux );

  //   const unsigned nbeta = override_nbins ? override_nbins : 100*(1<<vdoslux);//100 * 2^(vdoslux).
  //   const unsigned nalpha = ( override_nbins ? override_nbins : nbeta/2 );
  //   auto fixme = VDOS::determineAlphaBetaGridFromGn( gnexpn,
  //                                                    // Gn_asym, msd,
  //                                                    // upper_beta, upper_alpha,
  //                                                    nbeta, nalpha );
  // }

#if 0

#  if 0
  VectD ptsE0;//empty
  VectD betaGrid = setupBetaGrid( Gn_asym, upper_beta, vdoslux, override_nbins );
#  elif 0
  VectD betaGrid = setupBetaGridNEW( Gn_asym, msd, upper_beta, vdoslux, override_nbins);
#  else
  VectD betaGrid, ptsE0;
  std::tie(betaGrid, ptsE0)
    = setupBetaGridNEWNEW( Gn_asym, gnexpn.upper_beta, vdoslux, override_nbins);
#  endif
  const unsigned alpha_size = ( override_nbins ? override_nbins : betaGrid.size()/2 );
  VectD alphaGrid = setupAlphaGrid( gnexpn.Gn.kT(), gnexpn.msd, ptsE0, gnexpn.upper_alpha, alpha_size );

#endif


  //All done, now all that remains is to go through the (alpha,beta) pts in the
  //grid and use Sjolander's II.28 equation to calculate S(alpha,beta) there as
  //the sum of individual phonon orders:
#if 0
  //old way, no concurrency:
  auto sab = V2SKDetail::fillSABFromVDOS( Gn_asym, msd, alphaGrid, betaGrid, scaleGnContributionFct );
#else
  auto sab = V2SKDetail::fillSABFromVDOSConcurrent( Gn_asym, gnexpn.msd, alphaGrid, betaGrid, scaleGnContributionFct );
#endif
  double suggestedEmax = gnexpn.suggestedEmax;
  const auto max_phonon_order = Gn_asym.maxOrder().value();//FIXME: Verify that our thread strategy does not increase this beyond the actual algorithm
  if ( scaleGnContributionFct!=nullptr && scaleGnContributionFct(max_phonon_order) == 0.0 ) {
    //Caller might have essentially removed the last order(s), so it is unknown
    //how far the kernel can be used.
    suggestedEmax = 0.0;
  }

  if (V2SKDetail::s_verbose)
    NCRYSTAL_MSG("VDOS2SK created SK with vdos expansion order N="<<max_phonon_order
                 <<", Emax="<<suggestedEmax<<"eV, nalpha="<<alphaGrid.size()<< " nbeta="<<betaGrid.size());

  ScatKnlData out;
  out.alphaGrid = std::move(alphaGrid);
  out.betaGrid  = std::move(betaGrid);
  out.sab       = std::move(sab);
  out.temperature = vdosdata.temperature();
  out.boundXS = vdosdata.boundXS();
  out.elementMassAMU = vdosdata.elementMassAMU();
  out.knltype = ScatKnlData::KnlType::SAB;
  out.suggestedEmax = suggestedEmax;
  out.betaGridOptimised = true;//prevent beta-thickening code upon conversion to SABData
  return out;
}

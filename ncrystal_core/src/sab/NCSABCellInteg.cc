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

#include "NCrystal/internal/sab/NCSABCellInteg.hh"

namespace NC = NCrystal;
namespace NCS = NCrystal::SABUtils;

//Fixme: A lot of cleanup and consolidation still to happen in this file, once
//we are in a position to benchmark the effect in realistic usage conditions.

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {
    namespace {

      enum class AlphaInterpMethod { LOG, LIN };

      class BetaEdgeData final : private NoCopyMove {
        //Helper class which is used to keep track of S and log(S) at the
        //alpha-edges during loop over regions. This helps reduce number of
        //std::log calls.
        double m_a, m_S, m_logS;
        const CellData* m_cell;
        int m_soffset;
        void calcUpdateSlogS() {
          auto& c = *m_cell;
          auto s_and_logs = interpolate_loglin_fallbacklinlin_fast2
            ( c.a1, c.S[m_soffset], c.a2, c.S[m_soffset+1], m_a,
              c.logS[m_soffset], c.logS[m_soffset+1]);
          m_S = s_and_logs.first;
          m_logS = s_and_logs.second;
        }
      public:
        BetaEdgeData( const CellData* c, int soffset, double alpha )
          : m_a(alpha), m_cell(c), m_soffset(soffset)
        {
          nc_assert( alpha <= c->a2 );
          nc_assert( alpha > c->a1 );
          if ( m_a == c->a2 ) {
            m_S = c->S[m_soffset+1];
            m_logS = c->logS[m_soffset+1];
          } else {
            calcUpdateSlogS();
          }
        }
        void updateCurrentAlpha(double alpha)
        {
          auto& c = *m_cell;
          nc_assert( alpha < m_a );
          nc_assert( alpha < c.a2 );
          nc_assert( alpha >= c.a1 );
          m_a = alpha;
          if ( alpha == c.a1 ) {
            m_S = c.S[m_soffset];
            m_logS = c.logS[m_soffset];
          } else {
            calcUpdateSlogS();
          }
        }
        double getAlpha() const { return m_a; }
        double getS() const { return m_S; }
        double getLogS() const { return m_logS; }
      };

      struct SOfAlphaGrid final : private NoCopyMove {
        //Class which sets up an alpha grid like linspace(a1,a2,n) with
        //associated interpolated values of S and logS. In case of loglin
        //interpolation, a naive implementation would use n-2 std::exp calls to
        //achieve this, but the implementation here instead gets by with just
        //log2(n-1) calls to std::sqrt (i.e. 4 at n=17).
        using Method = AlphaInterpMethod;
        SOfAlphaGrid( Method meth,
                      double a1, double s1, double a2, double s2, unsigned nn )
          : n(nn)
        {
          nc_assert(a1>=0);
          nc_assert(a2>a1);
          nc_assert(s1>=0.0);
          nc_assert(s2>=0.0);
          nc_assert( ncmin(s1,s2)>0.0 || meth == Method::LIN );
          nc_assert( isOneOf(n,2u,3u,5u,9u,17u,33u) );
          a[0] = a1;
          S[0] = s1;
          const unsigned nm1 = n-1;
          a[nm1] = a2;
          S[nm1] = s2;
          if ( n == 2 )
            return;

          const double inv_nm1 = 1.0 / nm1;
          const double da = (a2-a1)*inv_nm1;
          for ( unsigned i = 1; i < nm1; ++i )
            a[i] = a1 + da * i;

          if ( meth == Method::LIN ) {
            //linear
            final_k = s2-s1;
            double ds = final_k*inv_nm1;
            for ( unsigned i = 1; i < nm1; ++i )
              S[i] = s1 + ds * i;
            return;
          }

          //Example of the following alg for n=9, how power of k can be achieved
          //by multiplying together sqrt^n(k):
          //
          //k^(1/8) = sqrt(sqrt(sqrt(k)))
          //0 1 2 3 4 5 6 7 8 //n=9, nbins=8
          //        @ @ @ @   //k^(1/2)   4 1 times, stepsize = 4, areas=1
          //    @ @     @ @   //k^(1/4)   2 2 times, stepsize = 2, areas=2
          //  @   @   @   @   //k^(1/8)   1 4 times, stepsize = 1, areas=4
          double k = s2/s1;
          for ( unsigned i = 1; i < nm1; ++i )
            S[i] = s1;
          unsigned stepsize = nm1;
          unsigned areas = 1;
          while ( stepsize /= 2 ) {
            k = std::sqrt(k);
            unsigned i = stepsize;
            for ( unsigned ia = 0; ia<areas; ++ia ) {
              for ( unsigned j=0 ; j < stepsize; ++j, ++i )
                S[i] *= k;
              i += stepsize;
            }
            areas *= 2;
          }
          final_k = k;
        }
        static constexpr unsigned nmax = 33;
        double a[nmax];
        double S[nmax];
        std::size_t n;
        double final_k;//needed for adaptive alg
      };

      class IntegrandOfA final : private NoCopyMove {
      public:

        struct Input {
          AlphaInterpMethod interpAtB1, interpAtB2;
          double a1, a2, b1, b2, s11, s12, s21, s22;
          double E_div_kT;
          bool is_bounded_by_betaminus;
          bool is_bounded_by_betaplus;
        };
        IntegrandOfA( const Input& i )
          : m_interpAtB1(i.interpAtB1),
            m_interpAtB2(i.interpAtB2),
            m_4e(4.0 * i.E_div_kT),
            m_a2minusa1( i.a2 - i.a1 ),
            m_b1(i.b1), m_b2(i.b2), m_invdb(1.0/(i.b2-i.b1)),
            m_is_bounded_by_betaminus(i.is_bounded_by_betaminus),
            m_is_bounded_by_betaplus(i.is_bounded_by_betaplus),
            m_is_bounded_by_both( i.is_bounded_by_betaminus
                                  && i.is_bounded_by_betaplus )
        {
          //only for crossed cells:
          nc_assert( m_is_bounded_by_betaminus ||  m_is_bounded_by_betaplus );
          nc_assert( (i.b2-i.b1) > 0.0 );
          //We always initialise to 4 levels (16 bins, 17 pts) with the points
          //in reverse order (the reverse order makes it easy to ignore the
          //point at a2 by just skipping the first element - even after we grow
          //and append additional points).

          //Fixme: we are simply using SOfAlphaGrid for the initialisation here,
          //not sure if it is too wasteful.
          SOfAlphaGrid s_b1( m_interpAtB1, i.a1, i.s11, i.a2, i.s12, 17 );
          SOfAlphaGrid s_b2( m_interpAtB2, i.a1, i.s21, i.a2, i.s22, 17 );
          m_stepcache_at_b1 = s_b1.final_k;
          m_stepcache_at_b2 = s_b2.final_k;
          for ( unsigned k = 0; k < 17; ++k ) {
            unsigned j = 16-k;
            m_data.emplace_back(s_b1.a[j],s_b1.S[j],s_b2.S[j]);
          }
          nc_assert(m_data.size()==17);
          nc_assert(m_data.at(0).alpha>m_data.at(1).alpha);
        }

        void evalInitialF17(double* fvals) const {
          nc_assert(m_data.size()>=17);
          for ( unsigned i = 0; i < 17; ++i )
            fvals[i] = contrib( m_data[16-i] );
        }

        double growNextLevelAndCollectContribSum()
        {
          auto nold = m_data.size();
          grow();
          auto nnew = m_data.size();
          StableSum sum;
          for ( auto i = nold; i < nnew; ++i )
            sum.add(contrib(vectAt(m_data,i)));
          return sum.sum();
        }

      private:
        struct AlphaSlice {
          AlphaSlice(double a,double sb1,double sb2)
            : alpha(a), s_at_b1(sb1), s_at_b2(sb2) {}
          double alpha;
          double s_at_b1;
          double s_at_b2;
        };

        double contrib(const AlphaSlice& slice) const
        {
          const double a = slice.alpha;
          const double dbpm = std::sqrt( m_4e * a );//NB: Most expensive
                                                    //per-point calc in this
                                                    //line? (fixme: revisit with profiling)
          const double bl( m_is_bounded_by_betaminus ? a - dbpm : m_b1 );
          const double bu( m_is_bounded_by_betaplus ? a + dbpm : m_b2 );
          //To find the contribution we integrate S(a,b) over [bl,bu]. This is
          //easy, since we always interpolate linearly in b:
          const double bmiddle( m_is_bounded_by_both ? a : (bu+bl)*0.5 );
          const double rb = (bmiddle-m_b1)*m_invdb;
          const double smiddle = slice.s_at_b1*(1.0-rb)+slice.s_at_b2*rb;
          const double bumbl( m_is_bounded_by_both ? 2.0*dbpm : bu-bl );
          return bumbl * smiddle;
        }

        double growK( AlphaInterpMethod aim,
                      double inv_newnbins,
                      double& stepcache ) const
        {
          if ( aim == AlphaInterpMethod::LOG )
            return stepcache = std::sqrt(stepcache);
          else
            return stepcache*inv_newnbins;
        }

        void grow() {
          const std::size_t Nold = m_data.size();
          nc_assert( Nold < 100000 );

          //Figure out "deltaS" for this step (multiplicative or additive)
          double inv_newnbins = 0.5/(Nold-1);
          double k1 = growK(m_interpAtB1,inv_newnbins,m_stepcache_at_b1);
          double k2 = growK(m_interpAtB2,inv_newnbins,m_stepcache_at_b2);
          nc_assert(m_a2minusa1>0.0);
          const double da = m_a2minusa1*inv_newnbins;
          nc_assert(da>0.0);
          m_data.reserve_hint( 1025 );//reduce slow reallocs
          for ( std::size_t i = 1; i < Nold; ++i ) {
            const auto& ref = vectAt(m_data,i);
            m_data.emplace_back( ref.alpha + da,
                                 ( m_interpAtB1 == AlphaInterpMethod::LOG
                                   ? ref.s_at_b1*k1
                                   : ref.s_at_b1+k1 ),
                                 ( m_interpAtB2 == AlphaInterpMethod::LOG
                                   ? ref.s_at_b2*k2
                                   : ref.s_at_b2+k2 ) );
          }
        }

        SmallVector<AlphaSlice,33> m_data;
        AlphaInterpMethod m_interpAtB1;
        AlphaInterpMethod m_interpAtB2;
        double m_4e;
        double m_a2minusa1;
        double m_b1, m_b2, m_invdb;
        double m_stepcache_at_b1, m_stepcache_at_b2;
        bool m_is_bounded_by_betaminus;
        bool m_is_bounded_by_betaplus;
        bool m_is_bounded_by_both;
      };

      class R17Adaptive final : public Romberg {
        mutable IntegrandOfA m_iofa;
        unsigned m_minlvl;
        unsigned m_maxlvl;
        double m_prec;
      public:
        R17Adaptive( const IntegrandOfA::Input& inp,
                     unsigned minlvl, unsigned maxlvl, double prec)
          : m_iofa(inp), m_minlvl(minlvl), m_maxlvl(maxlvl), m_prec(prec) {}

        double evalFunc(double) const override
        {
          nc_assert_always(false);
          return 0.0;
        }

        void evalFuncMany(double* fvals, unsigned n,
                          double, double) const override
        {
          (void)n;
          nc_assert(n==17);
          m_iofa.evalInitialF17(fvals);
        }

        double evalFuncManySum(unsigned, double, double) const override
        {
          return m_iofa.growNextLevelAndCollectContribSum();
        }

        bool accept( unsigned level, double prev_est, double est,
                     double, double ) const override
        {
          return ( level==m_maxlvl
                   || ( level>=m_minlvl
                        && ncabs(prev_est-est) <= est*m_prec ) );
        }

      };

      struct DecodedIntegScheme
      {
        unsigned npts;
        bool is_flex;
        bool is_romberg;
        bool is_simpson;
        DecodedIntegScheme( StdLogLinCellIntegrator::IntegrationScheme s )
          : npts( static_cast<std::uint_fast32_t>(s) & 0x00FFFF ),
            is_flex( static_cast<std::uint_fast32_t>(s) & 0x100000 ),
            is_romberg( static_cast<std::uint_fast32_t>(s) & 0x010000 ),
            is_simpson( static_cast<std::uint_fast32_t>(s) & 0x040000 )
        {
          constexpr static std::uint_fast32_t test17 = 0x000011;
          static_assert( test17 == 17, "");
          //check we decoded exactly one type (trapezoidal when all others are
          //unset).
#ifndef NDEBUG
          const bool is_trapez = static_cast<std::uint_fast32_t>(s) & 0x020000;
#endif
          nc_assert( (int(is_romberg)+int(is_simpson)
                      +int(is_flex)+int(is_trapez))==1 );
          //only trapez has n=2 and only flex has n=65:
          nc_assert(isOneOf(npts,2u,3u,5u,9u,17u,33u,65u));
          nc_assert( npts!=65u || is_flex );
          nc_assert( npts!=2u || is_trapez );
          nc_assert( npts!=3u || (is_trapez||is_simpson) );
        }
        StdLogLinCellIntegrator::IntegrationScheme encode() const//fixme: used?
        {
          std::uint_fast32_t v = static_cast<std::uint_fast32_t>(npts);
          if (is_romberg) {
            v |= 0x010000;
            nc_assert(isOneOf(npts,5u,9u,17u,33u));
          } else if (is_flex) {
            v |= 0x100000;
            nc_assert(isOneOf(npts,5u,9u,17u,33u,65u));
          } else if (is_simpson) {
            v |= 0x040000;
            nc_assert(isOneOf(npts,3u,5u,9u,17u,33u));
          } else {
            v |= 0x020000;//trapez
            nc_assert(isOneOf(npts,2u,3u,5u,9u,17u,33u));
          }
          return static_cast<StdLogLinCellIntegrator::IntegrationScheme>(v);
        }
      };

      static void impl_numIntRegion( const CellData& entire_cell,
                                     const CellData& subcell,
                                     double E_div_kT,
                                     StdLogLinCellIntegrator::IntegrationScheme
                                     scheme_encoded,
                                     bool is_bounded_by_betaminus,
                                     bool is_bounded_by_betaplus,
                                     StableSum& tgt )
      {
        nc_assert(is_bounded_by_betaminus||is_bounded_by_betaplus);
        const DecodedIntegScheme scheme{scheme_encoded};

        const CellData& cs = subcell;
        const CellData& ce = entire_cell;

        const SOfAlphaGrid::Method method_b1( ncmin(ce.S[0],ce.S[1])
                                              ? SOfAlphaGrid::Method::LOG
                                              : SOfAlphaGrid::Method::LIN );
        const SOfAlphaGrid::Method method_b2( ncmin(ce.S[2],ce.S[3])
                                              ? SOfAlphaGrid::Method::LOG
                                              : SOfAlphaGrid::Method::LIN );

        //Guard against enourmous log-scale differences along alpha.
        bool use_romberg_adaptive(false);
        bool use_romberg_fixed(scheme.is_romberg);
        if ( scheme.is_flex ) {
          if ( scheme.npts >= 33 ) {
            use_romberg_adaptive = true;
          } else {
            //fixme: in some extreme cases (s1=1,s2=1e200), it might be more
            //appropriate to instead split the cell, and do the more expensive
            //integration only very near s2=1e200 where all the contributions
            //are.
            constexpr double thr = 1e4;//fixme: tune!
            constexpr double recthr = 1/thr;
            use_romberg_adaptive
              = ( ( method_b1==SOfAlphaGrid::Method::LOG &&
                    !valueInInterval(cs.S[0]*recthr,cs.S[0]*thr,cs.S[1]) )
                  || ( method_b2==SOfAlphaGrid::Method::LOG &&
                       !valueInInterval(cs.S[2]*recthr,cs.S[2]*thr,cs.S[3]) ) );
            use_romberg_fixed = !use_romberg_adaptive;
          }
        }

        nc_assert( cs.b1 == ce.b1 );
        nc_assert( cs.b2 == ce.b2 );
        nc_assert( (cs.b2-cs.b1)>0.0 );

        //No matter the integration scheme, we must find the contribution at
        //each point of the alpha grid (skip this in case of adaptive Romberg
        //integration).
        double contrib_at_a[SOfAlphaGrid::nmax];
        if (!use_romberg_adaptive) {
          nc_assert(scheme.npts<=SOfAlphaGrid::nmax);
          const double invdb = 1.0/(cs.b2-cs.b1);
          SOfAlphaGrid sofa_at_b1( method_b1, cs.a1, cs.S[0], cs.a2, cs.S[1],
                                   scheme.npts );
          SOfAlphaGrid sofa_at_b2( method_b2, cs.a1, cs.S[2], cs.a2, cs.S[3],
                                   scheme.npts );
          double * itC = contrib_at_a;
          double * itCE = itC + scheme.npts;
          const double * itSb1 = sofa_at_b1.S;
          const double * itSb2 = sofa_at_b2.S;
          const double * itA = sofa_at_b2.a;
          double bl(cs.b1), bu(cs.b2);
          const bool is_bounded_on_both_sides ( is_bounded_by_betaminus
                                                && is_bounded_by_betaplus );
          for ( ; itC!=itCE; ++itC ) {
            double Sb1 = *(itSb1++);
            double Sb2 = *(itSb2++);
            double a = *(itA++);
            double dbpm = 2.0 * std::sqrt( E_div_kT * a );//fixme: can we avoid this??
            if ( is_bounded_by_betaminus )
              bl = a - dbpm;
            if ( is_bounded_by_betaplus )
              bu = ncmax(bl,a + dbpm);//ncmax as a safeguard against FP issues
            //To find the contribution we integrate S(a,b) over [bl,bu]. This is
            //easy, since we always interpolate linearly in b:
            const double bmiddle( is_bounded_on_both_sides ? a : (bu+bl)*0.5 );
#if 0
            double smiddle = Sb1 + (Sb2-Sb1)*(bmiddle-cs.b1)*invdb;
#else
            double rb = (bmiddle-cs.b1)*invdb;
            double smiddle = Sb1*(1.0-rb)+Sb2*(rb);
#endif
            nc_assert_always(bu-bl > -1e-6);//could be slightly negative due to
                                            //numerical instabilities
            const double bumbl( is_bounded_on_both_sides
                                ? 2.0*dbpm
                                : ncmax(0.0,bu-bl) );

            *itC = bumbl*smiddle;
            nc_assert_always(*itC >= 0.0);
            nc_assert_always(std::isfinite(*itC));
          }
        }

        if ( use_romberg_adaptive ) {
          //fixme: divert to other function (before the contrib_at_a is
          //initialised above)
          double contrib;
          IntegrandOfA::Input i;
          i.interpAtB1 = method_b1;
          i.interpAtB2 = method_b2;
          i.a1 = cs.a1;
          i.a2 = cs.a2;
          i.b1 = cs.b1;
          i.b2 = cs.b2;
          i.s11 = cs.S[0];
          i.s12 = cs.S[1];
          i.s21 = cs.S[2];
          i.s22 = cs.S[3];
          i.E_div_kT = E_div_kT;
          i.is_bounded_by_betaminus = is_bounded_by_betaminus;
          i.is_bounded_by_betaplus = is_bounded_by_betaplus;
          unsigned minlvl, maxlvl;
          double prec;
          //fixme: revisit values below
          switch( scheme.npts ) {
          case 5: prec = 5e-3; minlvl=2; maxlvl=8; break;
          case 9: prec = 5e-4; minlvl=3; maxlvl=9; break;
          case 17: prec = 1e-6; minlvl=4; maxlvl=10; break;
          case 33: prec = 1e-9; minlvl=5; maxlvl=12; break;
          default:
            nc_assert(false);
          case 65: prec = 1e-12; minlvl=6; maxlvl=14; break;
          };
          R17Adaptive r17adapt(i,minlvl,maxlvl,prec);
          contrib = r17adapt.integrate(0.0,1.0);
          tgt.add( contrib*cs.a2 );
          tgt.add( -contrib*cs.a1 );
        } else if ( use_romberg_fixed ) {
          double contrib;
          nc_assert_always(contrib_at_a[0] >= 0.0);
          if ( scheme.npts<17 ) {
            if ( scheme.npts==5 ) {
              contrib = Romberg::fixedOrderIntegration5pts(contrib_at_a);
            } else {
              nc_assert(scheme.npts==9);
              contrib = Romberg::fixedOrderIntegration9pts(contrib_at_a);
            }
          } else {
            if ( scheme.npts==17 ) {
              contrib = Romberg::fixedOrderIntegration17pts(contrib_at_a);
            } else {
              nc_assert(scheme.npts==33);
              contrib = Romberg::fixedOrderIntegration33pts(contrib_at_a);
            }
          }
          tgt.add( contrib*cs.a2 );
          tgt.add( -contrib*cs.a1 );
        } else if ( scheme.is_simpson ) {
          nc_assert( !use_romberg_adaptive );
          nc_assert_always(contrib_at_a[0] >= 0.0);
          nc_assert( scheme.npts%2==1 && scheme.npts>=3 );
          StableSum ss;
          const unsigned nbins = scheme.npts-1;
          const double k = 1.0/(3.0*nbins);
          const double k2 = k + k;
          const double k4 = k2 + k2;
          const double * itCB = contrib_at_a;
          const double * itCL = itCB + nbins;
          ss.add( k * (*itCB) );
          for ( const double * itC = itCB+1;
                itC < itCL;
                itC += 2 )
            ss.add( k4 * (*itC) );
          for ( const double * itC = itCB+2;
                itC < itCL;
                itC += 2 )
            ss.add( k2 * (*itC) );
          ss.add( k * (*itCL) );
          double contrib = ss.sum();
          tgt.add( contrib*cs.a2 );
          tgt.add( -contrib*cs.a1 );
        } else {
          //Trapezoidal
          nc_assert( !use_romberg_adaptive );
          nc_assert_always(contrib_at_a[0] >= 0.0);
          unsigned nbins = scheme.npts-1;
          double da = (cs.a2-cs.a1)/nbins;
          nc_assert_always(da>0.0);
          const double * itC = contrib_at_a;
          const double * itCL = itC + nbins;
          nc_assert_always(itCL < contrib_at_a + scheme.npts);
          nc_assert_always(*itC >= 0.0);
          tgt.add( 0.5 * da * (*itC++) );
          for (; itC!=itCL; ++itC )
            tgt.add( da * (*itC) );
          tgt.add( 0.5 * da * (*itCL) );
        }
      }
    }
  }
}

void NCS::StdLogLinCellIntegrator::integrateWithinKB( const CellData& c,
                                                      double E_div_kT,
                                                      StdLogLinCellIntegrator::
                                                      IntegrationScheme scheme,
                                                      StableSum& tgt )
{
  SABCellSurvey surv( c.a1, c.a2, c.b1, c.b2, E_div_kT );//fixme: in
                                                         //SABProcessor.cc we
                                                         //call this function
                                                         //from a context where
                                                         //we already have a
                                                         //SABCellSurvey!
  auto& regions = surv.regions();
  if ( regions.empty() )
    return;
  //We must find S and log(S) at beta=b1 or beta=b2 at the region's alpha
  //bounds. We take advantage of the fact that the regions are sorted from
  //high to low alpha, and that neighbouring regions in the list share an
  //alpha point.

  auto r0 = regions.front();
  BetaEdgeData atb1( &c, 0, r0.alpha_up );
  BetaEdgeData atb2( &c, 2, r0.alpha_up );

  CellData subcell;//subcell belonging to a given region.
  //subcells ALWAYS have the full beta-range (we can't split cells at fixed beta
  //in the given interpolation scheme).
  subcell.b1 = c.b1;
  subcell.b2 = c.b2;

  for ( auto iregion : ncrange(surv.regions().size() ) ) {
    auto& r = regions.at(iregion);
    nc_assert( r.alpha_up > r.alpha_low );
    //Update subcell data @ upper alpha edge:
    nc_assert( r.alpha_up == atb1.getAlpha() );
    nc_assert( r.alpha_up == atb2.getAlpha() );
    subcell.a2 = r.alpha_up;
    subcell.S[1] = atb1.getS();
    subcell.S[3] = atb2.getS();
    subcell.logS[1] = atb1.getLogS();
    subcell.logS[3] = atb2.getLogS();
    atb1.updateCurrentAlpha( r.alpha_low );
    atb2.updateCurrentAlpha( r.alpha_low );
    //Update subcell data @ lower alpha edge:
    subcell.a1 = r.alpha_low;
    subcell.S[0] = atb1.getS();
    subcell.S[2] = atb2.getS();
    subcell.logS[0] = atb1.getLogS();
    subcell.logS[2] = atb2.getLogS();

    //Now evaluate region:
    if ( !(r.is_bounded_by_betaminus || r.is_bounded_by_betaplus ) ) {
      //region is fully rectangular, extending over entire [b1,b2] range!
      integrateFullCell( subcell, tgt );
      continue;
    }

    impl_numIntRegion( c, subcell, E_div_kT, scheme,
                       r.is_bounded_by_betaminus,
                       r.is_bounded_by_betaplus,
                       tgt );
  }
}

NCS::StdLogLinCellIntegrator::IntegrationScheme
NCS::StdLogLinCellIntegrator::str2IntegScheme( StrView v )
{
  using IS = IntegrationScheme;
  if ( v.startswith("Flex") ) {
    if ( v=="Flex9" )
      return IS::Flex9;
    if ( v=="Flex5" )
      return IS::Flex5;
    if ( v=="Flex17" )
      return IS::Flex17;
    if ( v=="Flex33" )
      return IS::Flex33;
    if ( v=="Flex65" )
      return IS::Flex65;
  } else if ( v.startswith("Romberg") ) {
    if ( v=="Romberg9" )
      return IS::Romberg9;
    if ( v=="Romberg5" )
      return IS::Romberg5;
    if ( v=="Romberg17" )
      return IS::Romberg17;
    if ( v=="Romberg33" )
      return IS::Romberg33;
  } else if ( v.startswith("Trapez") ) {
    if ( v=="Trapez2" )
      return IS::Trapez2;
    if ( v=="Trapez3" )
      return IS::Trapez3;
    if ( v=="Trapez5" )
      return IS::Trapez5;
    if ( v=="Trapez9" )
      return IS::Trapez9;
    if ( v=="Trapez17" )
      return IS::Trapez17;
    if ( v=="Trapez33" )
      return IS::Trapez33;
  } else if ( v.startswith("Simpson") ) {
    if ( v=="Simpson3" )
      return IS::Simpson3;
    if ( v=="Simpson5" )
      return IS::Simpson5;
    if ( v=="Simpson9" )
      return IS::Simpson9;
    if ( v=="Simpson17" )
      return IS::Simpson17;
    if ( v=="Simpson33" )
      return IS::Simpson33;
  }
  NCRYSTAL_THROW2(BadInput,"Invalid integration scheme: \""<<v
                  <<"\" (should be one of \""
                  << allIntegSchemesAsStr() <<"\")");
  return IS::Default;
}

const char * NCS::StdLogLinCellIntegrator::allIntegSchemesAsStr()
{
  return "Flex5;Flex9;Flex17;Flex33;Flex65;"
    "Romberg5;Romberg9;Romberg17;Romberg33;"
    "Trapez2;Trapez3;Trapez5;Trapez9;Trapez17;Trapez33;"
    "Simpson3;Simpson5;Simpson9;Simpson17;Simpson33";
}
const char * NCS::StdLogLinCellIntegrator::integSchemeToStr( IntegrationScheme v )
{
  using IS = IntegrationScheme;
  switch ( v ) {
    //fixme: something shorter, so might be used in cfg strings? r33?
  case IS::Flex5:  return "Flex5";
  case IS::Flex9:  return "Flex9";
  case IS::Flex17: return "Flex17";
  case IS::Flex33: return "Flex33";
  case IS::Flex65: return "Flex65";
  case IS::Romberg5:  return "Romberg5";
  case IS::Romberg9:  return "Romberg9";
  case IS::Romberg17: return "Romberg17";
  case IS::Romberg33: return "Romberg33";
  case IS::Trapez2:   return "Trapez2";
  case IS::Trapez3:   return "Trapez3";
  case IS::Trapez5:   return "Trapez5";
  case IS::Trapez9:   return "Trapez9";
  case IS::Trapez17:  return "Trapez17";
  case IS::Trapez33:  return "Trapez33";
  case IS::Simpson3:  return "Simpson3";
  case IS::Simpson5:  return "Simpson5";
  case IS::Simpson9:  return "Simpson9";
  case IS::Simpson17: return "Simpson17";
  case IS::Simpson33: return "Simpson33";
  default:
    nc_assert_always(false&&"invalid IntegrationScheme");
    return "";
  };
}

std::ostream& NCS::operator<<( std::ostream& os, const CellData& c )
{
  os << "CellData(a1=" << fmt(c.a1)
     << ", a2=" << fmt(c.a2)
     << ", b1=" << fmt(c.b1)
     << ", b2=" << fmt(c.b2)
     << ", S11=" << fmt(c.S[0])
     << ", S12=" << fmt(c.S[1])
     << ", S21=" << fmt(c.S[2])
     << ", S22=" << fmt(c.S[3])
     << ", logS11=" << fmt(c.logS[0])
     << ", logS12=" << fmt(c.logS[1])
     << ", logS21=" << fmt(c.logS[2])
     << ", logS22=" << fmt(c.logS[3]) << ')';
  return os;
}

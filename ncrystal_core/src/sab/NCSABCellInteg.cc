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

#include "NCrystal/internal/sab/NCSABCellInteg.hh"//fixme: reconsider filename?

namespace NC = NCrystal;
namespace NCS = NCrystal::SABUtils;

namespace NCRYSTAL_NAMESPACE {
  namespace SABUtils {
    namespace {

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

      struct SOfAlphaGrid final : NoCopyMove {
        //Class which sets up an alpha grid like linspace(a1,a2,n) with
        //associated interpolated values of S and logS. In case of loglin
        //interpolation, a naive implementation would use n-2 std::exp calls to
        //achieve this, but the implementation here instead gets by with just
        //log2(n-1) calls to std::sqrt (i.e. 4 at n=17).
        SOfAlphaGrid( double a1, double s1, double a2, double s2, unsigned nn )
          : n(nn)
        {
          nc_assert(a1>=0);
          nc_assert(a2>a1);
          nc_assert(s1>=0.0);
          nc_assert(s2>=0.0);
          nc_assert(isOneOf(n,2,3,5,9,17,33));
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

          if ( ncmin(s1,s2) == 0.0 ) {
            //linear
            double ds = (s2-s1)*inv_nm1;
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
        }
        static constexpr unsigned nmax = 33;
        double a[nmax];
        double S[nmax];
        std::size_t n;
      };

      struct DecodedIntegScheme
      {
        unsigned npts;
        bool is_romberg;
        bool is_simpson;
        DecodedIntegScheme( StdLogLinCellIntegrator::IntegrationScheme s )
          : npts( static_cast<std::uint_fast32_t>(s) & 0x00FFFF ),
            is_romberg( static_cast<std::uint_fast32_t>(s) & 0x010000 ),
            is_simpson( static_cast<std::uint_fast32_t>(s) & 0x040000 )
        {
          constexpr static std::uint_fast32_t test17 = 0x000011;
          static_assert( test17 == 17, "");
          nc_assert( !(is_romberg&&is_simpson) );
          nc_assert(isOneOf(npts,2,3,5,9,17,33));
        }
        StdLogLinCellIntegrator::IntegrationScheme encode() const//fixme: used?
        {
          std::uint_fast32_t v = static_cast<std::uint_fast32_t>(npts);
          if (is_romberg) {
            v |= 0x010000;
            nc_assert(isOneOf(npts,5,9,17,33));
          } else if (is_simpson) {
            v |= 0x040000;
            nc_assert(isOneOf(npts,3,5,9,17,33));
          } else {
            v |= 0x020000;//trapez
            nc_assert(isOneOf(npts,2,3,5,9,17,33));
          }
          return static_cast<StdLogLinCellIntegrator::IntegrationScheme>(v);
        }
      };

      static void impl_numIntRegion( const CellData& c,
                                     double E_div_kT,
                                     StdLogLinCellIntegrator::IntegrationScheme
                                     scheme_encoded,
                                     bool is_bounded_by_betaminus,
                                     bool is_bounded_by_betaplus,
                                     StableSum& tgt )
      {
        nc_assert(is_bounded_by_betaminus||is_bounded_by_betaplus);
        const DecodedIntegScheme scheme{scheme_encoded};

        SOfAlphaGrid sofa_at_b1( c.a1, c.S[0], c.a2, c.S[1], scheme.npts );
        SOfAlphaGrid sofa_at_b2( c.a1, c.S[2], c.a2, c.S[3], scheme.npts );

        const double invdb = 1.0/(c.b2-c.b1);

        //No matter the integration scheme, we must find the contribution at
        //each point of the alpha grid.
        double contrib_at_a[SOfAlphaGrid::nmax];
        nc_assert_always(scheme.npts<=SOfAlphaGrid::nmax);//fixme _always
        {
          double * itC = contrib_at_a;
          double * itCE = itC + scheme.npts;
          const double * itSb1 = sofa_at_b1.S;
          const double * itSb2 = sofa_at_b2.S;
          const double * itA = sofa_at_b2.a;
          double bl(c.b1), bu(c.b2);
          for ( ; itC!=itCE; ++itC ) {
            double Sb1 = *(itSb1++);
            double Sb2 = *(itSb2++);
            double a = *(itA++);
            double dbpm = 2.0 * std::sqrt( E_div_kT * a );
            if ( is_bounded_by_betaminus )
              bl = a - dbpm;
            if ( is_bounded_by_betaplus )
              bu = a + dbpm;
            //To find the contribution we integrate S(a,b) over [bl,bu]. This is
            //easy, since we always interpolate linearly in b:
            double bmiddle = (bu+bl)*0.5;
            double smiddle = Sb1 + (Sb2-Sb1)*(bmiddle-c.b1)*invdb;
            *itC = (bu-bl)*smiddle;
          }

          //fixme
          // for ( auto i : ncrange(scheme.npts) )
          //   NCRYSTAL_MSG("TKTEST cpp contrib(i="<<i<<",a="<<fmtg(sofa_at_b2.a[i])<<") = "<<fmtg(contrib_at_a[i]))

        }
        if ( scheme.is_romberg ) {
          double contrib;
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
          tgt.add( contrib*(c.a2-c.a1) );
        } else if ( scheme.is_simpson ) {
          nc_assert( scheme.npts%2==1 && scheme.npts>=3 );
          const unsigned nbins = scheme.npts-1;
          const double k = (c.a2-c.a1)/(3.0*nbins);
          const double k2 = k + k;
          const double k4 = k2 + k2;
          const double * itCB = contrib_at_a;
          const double * itCL = itCB + nbins;
          tgt.add( k * (*itCB) );
          for ( const double * itC = itCB+1;
                itC < itCL;
                itC += 2 )
            tgt.add( k4 * (*itC) );
          for ( const double * itC = itCB+2;
                itC < itCL;
                itC += 2 )
            tgt.add( k2 * (*itC) );
          tgt.add( k * (*itCL) );
        } else {
          //Trapezoidal
          unsigned nbins = scheme.npts-1;
          double da = (c.a2-c.a1)/nbins;
          const double * itC = contrib_at_a;
          const double * itCL = itC + nbins;
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
  SABCellSurvey surv( c.a1, c.a2, c.b1, c.b2, E_div_kT );
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

    //Needs careful numerical integration!
    //fixme double pre = tgt.sum();
    impl_numIntRegion( subcell, E_div_kT, scheme,
                       r.is_bounded_by_betaminus,
                       r.is_bounded_by_betaplus,
                       tgt );
    //fixme NCRYSTAL_MSG("TKTEST region numint gives: "<<fmtg(tgt.sum()-pre));
  }
}

NCS::StdLogLinCellIntegrator::IntegrationScheme
NCS::StdLogLinCellIntegrator::str2IntegScheme( StrView v )
{
  using IS = IntegrationScheme;
  if ( v.startswith("Romberg") ) {
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
  return "Romberg5;Romberg9;Romberg17;Romberg33;"
    "Trapez2;Trapez3;Trapez5;Trapez9;Trapez17;Trapez33;"
    "Simpson3;Simpson5;Simpson9;Simpson17;Simpson33";
}
const char * NCS::StdLogLinCellIntegrator::integSchemeToStr( IntegrationScheme v )
{
  using IS = IntegrationScheme;
  switch ( v ) {
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

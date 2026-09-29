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

#include "NCrystal/internal/phys_utils/NCSCTUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"

namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {

    Temperature sct_checkedTeff( Temperature t, Temperature teff )
    {
      //Physically Teff>=T always holds for harmonic systems, but allow
      //(and clamp) a tiny numerical dip below T:
      t.validate();
      teff.validate();
      if ( teff.dbl() >= t.dbl() )
        return teff;
      if ( teff.dbl() >= 0.999 * t.dbl() )
        return t;
      NCRYSTAL_THROW2(BadInput,"SCT model got effective temperature "
                      <<teff<<" significantly below the actual material "
                      "temperature "<<t<<" which is unphysical.");
    }

    //Geometric integration segments for the upscatter correction
    //integral (integrand vanishes at beta=0 and decays at least as fast
    //as exp(-beta), so [0,60] is beyond double precision exhaustive and
    //geometric segments keep low-order quadrature accurate):
    //(the finely divided first segments also resolve the narrow
    //1/(r-1) rise scale of the suppression factor at large r):
    constexpr double sct_corr_edges[11]
      = { 0.0, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 60.0 };

    inline double sct_corr_integrand( double c, double beta,
                                      double invA, double rm1 )
    {
      return - SCTXSProvider::evalFGBetaKernel( c, beta, invA )
        * std::expm1( -rm1*beta );
    }

    inline double sct_corr_betascale( double c, double invA )
    {
      //Characteristic beta-scale of the integrand: for heavy targets
      //the upscatter kernel decays like exp(-beta*(1+A/4)) (from
      //optimising the Gaussian over the allowed alpha range at large
      //beta), plus a Doppler width term ~2*sqrt(c/A). Scaling the
      //integration segments by s=min(1,4/(4+A)+2*sqrt(c/A)) both
      //resolves the narrow heavy-target integrand and preserves the
      //coverage guarantee, since s*decayrate >= 1 by construction (so
      //the covered range always corresponds to at least the exp(-60)
      //level):
      const double A = 1.0 / invA;
      return ncmin( 1.0, 4.0/(4.0+A) + 2.0*std::sqrt( c*invA ) );
    }

  }
}

double NC::SCTXSProvider::evalFGBetaKernel( double c, double beta,
                                            double invA )
{
  //See the header for the derivation and conventions (u=sqrt(alpha_E),
  //value independent of alpha convention, only needed for beta>=0):
  nc_assert( c > 0.0 && beta >= 0.0 && invA > 0.0 );
  auto alim = getAlphaLimits( c, beta );
  if ( !( alim.second > alim.first ) )
    return 0.0;
  const double um = std::sqrt( alim.first * invA );
  const double up = std::sqrt( alim.second * invA );
  auto yfct = [beta]( double u )
  {
    return 0.5 * ( u + beta / u );
  };
  auto zfct = [beta]( double u )
  {
    return 0.5 * ( u - beta / u );
  };
  //NB: erfcdiff(a,b) = erf(b)-erf(a). y(u) is not monotone (minimum at
  //u=sqrt(beta)), but erfcdiff needs no argument ordering. Both
  //erfcdiff and ncerf are the portable libm-free implementations:
  const double term_y = ( um > 0.0
                          ? erfcdiff( yfct(um), yfct(up) )
                          : ncerf( yfct(up) ) );
  const double term_z = ( um > 0.0
                          ? erfcdiff( zfct(um), zfct(up) )
                          : ncerf( zfct(up) ) );
  return 0.5 * ( term_y + std::exp( -beta ) * term_z );
}

NC::SCTXSProvider::SCTXSProvider( Temperature temp_k,
                                  Temperature teff,
                                  AtomMass mass,
                                  SigmaBound sigma,
                                  unsigned table_npts )
  : m_xsprovider( sct_checkedTeff(temp_k,teff), mass, sigma ),
    m_teff( sct_checkedTeff(temp_k,teff) ),
    m_r( m_teff.dbl() == temp_k.dbl()
         ? 1.0 : m_teff.dbl() / temp_k.dbl() ),
    m_invkTeff( 1.0 / m_teff.kT() ),
    m_invA( 1.0 / mass.relativeToNeutronMass() ),
    m_sb( sigma.dbl() ),
    m_csat( 35.0 * ncmax( m_invA, 1.0/m_invA ) ),
    m_clow( 1e-3 )
{
  nc_assert( m_r >= 1.0 );
  if ( m_r > 1.0 ) {
    //Precalculate the upscatter correction integral versus ln(c) over
    //the non-saturated region (below m_clow direct quadrature is used,
    //which in practice never happens). Node values are evaluated with
    //cheap fixed-order Romberg (accuracy far beyond the SCT model's
    //own percent-level physics accuracy) and interpolated with a
    //natural cubic spline. NB: a PCHIP variant was tried and rejected
    //on structural grounds: its shape-limiting is nonlinear in the
    //node data, so a limiter-branch flip during an otherwise smooth
    //temperature/parameter scan could introduce kink artifacts at the
    //interpolation-error scale (fine-grained scan probes showed both
    //variants artifact-free in the cases tried, but only the spline,
    //being linear in the data, excludes the failure mode by
    //construction -- and it needs ~3x fewer nodes for a given
    //accuracy anyway). The saturated endpoint is pinned to the exact
    //analytic value and zero slope, making the c=csat branch seam in
    //evalUpscatterCorr exactly continuous:
    m_lutlna = std::log( m_clow );
    m_lutlnb = std::log( m_csat );
    nc_assert_always( table_npts >= 8 && table_npts <= 100000 );
    const unsigned npts = table_npts;
    VectD yvals;
    yvals.reserve( npts );
    const double dx = ( m_lutlnb - m_lutlna ) / ( npts - 1 );
    for ( auto k : ncrange(npts) ) {
      const double x = m_lutlna + k * dx;
      yvals.push_back( corrIntegralFixedOrder( std::exp( x ) ) );
    }
    yvals.back() = ( m_r - 1.0 ) / m_r;
    const double h = 1e-3 * dx;
    const double fprime_a
      = ( corrIntegralFixedOrder( std::exp( m_lutlna + h ) )
          - yvals.front() ) / h;
    m_corrlut.set( yvals, m_lutlna, m_lutlnb, fprime_a, 0.0,
                   "sctupscattercorr",
                   "SCT upscatter correction integral vs ln(c)" );
  }
}

double NC::SCTXSProvider::corrIntegralFixedOrder( double c ) const
{
  //Same integral as evalUpscatterCorrIntegral, but with a fixed
  //17-point Romberg rule per segment (deterministic, adaptivity-free
  //cost for the lookup-table construction):
  nc_assert( c > 0.0 && m_r > 1.0 );
  const double rm1 = m_r - 1.0;
  const double bscale = sct_corr_betascale( c, m_invA );
  double fvals[17];
  StableSum sum;
  for ( auto i : ncrange(10) ) {
    const double blo = bscale * sct_corr_edges[i];
    //NB: no tail early-exit here (unlike evalUpscatterCorrIntegral):
    //a data-dependent segment count would make the node values
    //discontinuous functions of (r,A,c), spoiling the smoothness of
    //parameter scans for a negligible saving:
    const double w = bscale * sct_corr_edges[i+1] - blo;
    for ( auto k : ncrange(17) )
      fvals[k] = sct_corr_integrand( c, blo + w * ( k / 16.0 ),
                                     m_invA, rm1 );
    sum.add( w * Romberg::fixedOrderIntegration17pts( fvals ) );
  }
  return sum.sum();
}

NC::SCTXSProvider::~SCTXSProvider() = default;

double NC::SCTXSProvider::evalUpscatterCorrIntegral( double c,
                                                     double prec,
                                                     unsigned maxlvl ) const
{
  nc_assert( c > 0.0 && m_r > 1.0 && prec > 0.0 );
  const double rm1 = m_r - 1.0;
  const double invA = m_invA;
  auto integrand = [c,rm1,invA]( double beta )
  {
    return sct_corr_integrand( c, beta, invA, rm1 );
  };
  const double bscale = sct_corr_betascale( c, invA );
  StableSum sum;
  for ( auto i : ncrange(10) ) {
    //The remaining contribution is bounded by 2*exp(-beta_lo) (the
    //kernel obeys I(beta)<=1.5*exp(-beta) for beta>=1), so the tail
    //segments can be skipped once they cannot matter at the requested
    //precision:
    const double blo = bscale * sct_corr_edges[i];
    if ( i > 1 && 2.0 * std::exp( -blo ) < prec * sum.sum() )
      break;
    sum.add( integrateRombergFlex( integrand, blo,
                                   bscale * sct_corr_edges[i+1], prec,
                                   3, maxlvl ) );
  }
  return sum.sum();
}

double NC::SCTXSProvider::evalUpscatterCorr( double c ) const
{
  nc_assert( c > 0.0 );
  if ( m_r == 1.0 )
    return 0.0;
  if ( c >= m_csat ) {
    //Kernel deviations from its exact c->infinity limit exp(-beta) are
    //O(erfc(sqrt(c*A))+erfc(sqrt(c/A))), i.e. negligible here, and the
    //correction integral has the closed form (r-1)/r (cf. header):
    return ( m_r - 1.0 ) / m_r;
  }
  if ( c >= m_clow )
    return m_corrlut.eval( ncclamp( std::log(c),
                                    m_lutlna, m_lutlnb ) );
  return evalUpscatterCorrIntegral( c );
}

NC::CrossSect NC::SCTXSProvider::crossSection( NeutronEnergy ekin ) const
{
  //sigma_SCT(E) = sigma_FG(Teff)(E) - (sigma_b*A/(4c)) * corr(c), cf.
  //the header (c=E/kTeff, and the S-integral is carried out in the
  //mass-scaled alpha convention, hence the factor of A):
  const auto xsfg = m_xsprovider.crossSection(ekin);
  if ( m_r == 1.0 )
    return xsfg;
  const double c = ekin.dbl() * m_invkTeff;
  nc_assert( c > 0.0 );
  const double corr = evalUpscatterCorr( c );
  const double xs = xsfg.dbl() - m_sb * corr / ( 4.0 * c * m_invA );
  nc_assert( xs > 0.0 && corr >= 0.0 );
  return CrossSect{ xs };
}

NC::SCTSampler::SCTSampler( NeutronEnergy ekin,
                            Temperature temp_k,
                            Temperature teff,
                            AtomMass mass )
  : m_fgs( ekin, sct_checkedTeff(temp_k,teff), mass ),
    m_r( sct_checkedTeff(temp_k,teff).dbl() == temp_k.dbl()
         ? 1.0 : teff.dbl() / temp_k.dbl() )
{
  nc_assert( m_r >= 1.0 );
}

NC::SCTSampler::~SCTSampler() = default;

NC::PairDD NC::SCTSampler::sampleAlphaBeta( RNG& rng ) const
{
  //Free-gas sampling at Teff (in its own kTeff-based units), scaled by
  //r=Teff/T to actual-kT units, with upscatter candidates additionally
  //accepted with the detailed-balance-at-T probability
  //exp(-beta*(r-1)) -- cf. the header:
  if ( m_r == 1.0 )
    return m_fgs.sampleAlphaBeta( rng );
  const double rm1 = m_r - 1.0;
#ifndef NDEBUG
  int nloop = 0;
#endif
  while ( true ) {
    nc_assert( nloop++ < 10000 );
    auto ab = m_fgs.sampleAlphaBeta( rng );
    if ( ab.second <= 0.0
         || rng.generate() < std::exp( -ab.second * rm1 ) )
      return { ab.first * m_r, ab.second * m_r };
  }
}

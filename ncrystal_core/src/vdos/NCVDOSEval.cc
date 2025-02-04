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

#include "NCrystal/internal/vdos/NCVDOSEval.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCString.hh"

namespace NC=NCrystal;

namespace NCRYSTAL_NAMESPACE {
  namespace {
    static std::atomic<bool> s_verbose_vdoseval( ncgetenv_bool("DEBUG_PHONON") );

    constexpr double detail_xcothx_taylor_threshold = 0.1;

    double safe_xcothx( double x )
    {
      // The function x*coth(x) = x/tanh(x) must be evaluated with a Taylor expansion near x=0.
      // For the record we simply got the Taylor coefficients with sagemath:
      // > sage: f=x*coth(x)
      // First investigating number of orders needed for u=0.1 with command (with a few orders added for safety):
      // > sage: u=0.1;(f-f.taylor(x,0,12))(x=u).n()
      // Then generate the coefficients with:
      // > sage:  print( '\n'.join(('constexpr double c%i = %s;'%(c[1],str(c[0]))).replace('/','./') for c in (f.taylor(x,0,14)).coefficients()))
      if ( x < detail_xcothx_taylor_threshold ) {
        constexpr double c0 = 1;
        constexpr double c2 = 1./3;
        constexpr double c4 = -1./45;
        constexpr double c6 = 2./945;
        constexpr double c8 = -1./4725;
        constexpr double c10 = 2./93555;
        constexpr double c12 = -1382./638512875.;
        constexpr double c14 = 4./18243225.;
        const double y = x*x;
        return c0+y*(c2+y*(c4+y*(c6+y*(c8+y*(c10+y*(c12+y*c14))))));
      } else {
        return x / std::tanh(x);
      }
    }

    constexpr double detail_x3cothx_taylor_threshold = 0.1;
    double safe_x3cothx( double x ) {
      // The function x^3*coth(x) = x^3/tanh(x) must be evaluated with a Taylor
      // expansion near x=0. See safe_xcothx for remarks on how to get the
      // Taylor coefficients easily.

      if ( x < detail_x3cothx_taylor_threshold ) {
        constexpr double c2 = 1;
        constexpr double c4 = 1./3;
        constexpr double c6 = -1./45;
        constexpr double c8 = 2./945;
        constexpr double c10 = -1./4725;
        constexpr double c12 = 2./93555;
        constexpr double c14 = -1382./638512875;
        constexpr double c16 = 4./18243225;
        constexpr double c18 = -3617./162820783125;
        constexpr double c20 = 87734./38979295480125;
        const double y = x*x;
        return y*(c2+y*(c4+y*(c6+y*(c8+y*(c10+y*(c12+y*(c14+y*(c16+y*(c18+y*(c20))))))))));
      } else {
        return x*x*x / std::tanh(x);
      }
    }

    constexpr double detail_xdivsinhx_taylor_threshold = 0.07;
    double safe_xdivsinhx( double x ) {
      // The function x/sinh(x) = x / csch(x) must be evaluated with a Taylor
      // expansion near x=0.
      if ( x < detail_xdivsinhx_taylor_threshold ) {
        constexpr double c2 = -1./6.;
        constexpr double c4 = 7./360.;
        constexpr double c6 = -31./15120.;
        constexpr double c8 = 127./604800.;
        const double y = x*x;
        return 1.0+y*(c2+y*(c4+y*(c6+y*c8)));
      } else {
        return x / std::sinh(x);
      }
    }
#ifndef NDEBUG
    static bool dummy_val_xcothx = []() {
      //Sanity check overlap at the threshold:
      nc_assert( floateq( safe_xcothx(detail_xcothx_taylor_threshold*(1-1e-15)),safe_xcothx(detail_xcothx_taylor_threshold*(1+1e-15)),1e-14,1e-99) );
      nc_assert( floateq( safe_x3cothx(detail_x3cothx_taylor_threshold*(1-1e-15)),safe_x3cothx(detail_x3cothx_taylor_threshold*(1+1e-15)),1e-14,1e-99) );
      nc_assert( floateq( safe_xdivsinhx(detail_xdivsinhx_taylor_threshold*(1-1e-15)),safe_xdivsinhx(detail_xdivsinhx_taylor_threshold*(1+1e-15)),1e-14,1e-99) );
      return true;
    }();
#endif
  }

}

void NC::VDOSEval::enableVerboseOutput(bool status)
{
#ifndef NDEBUG
  (void)dummy_val_xcothx;
#endif
  s_verbose_vdoseval = status;
}

bool NC::VDOSEval::verboseOutputEnabled()
{
  return s_verbose_vdoseval;
}

template <class Fct, class TStableSum>
void NC::VDOSEval::integrateBinsWithFunction( Fct f, TStableSum& sum ) const
{
  //////////////////////////////////////////////////////////////////////////////
  // With VDOS=rho(E), integrate rho(E)*f(E) over [emin,emax] (a.k.a. the "binned
  // part" of the density function). This differs from an integral over
  // [0,infinity] only in that it is missing the parabolic part over
  // [0,emin], which should be handled separately by the calling code.
  //
  // Since all callers will have to do their own custom integration over the
  // parabolic part [0,emin] where the density is given by a parabola: rho(E) =
  // m_k * E^2 with a custom function f(E), we derive the general formula for
  // doing that in a dimensionless variable u=E/2kT here. So with ulim =
  // m_emin/2kT:
  //
  // parabolic_part = integral_0^(m_emin){ m_k*E^2 * f(E) }dE
  //                = integral_0^(ulim){ m_k* (2kT)^2u^2*f(E=u*2kt)2kT }du
  //                = m_k* (2kT)^3 * integral_0^(ulim){ u^2*f(E=u*2kT)}du
  //////////////////////////////////////////////////////////////////////////////

  static_assert( std::is_same<TStableSum,StableSum>::value, "" );

  nc_assert( m_emin > 0.0 );


  //Contribution in [m_emin,m_emax], considering individual bins (since the
  //integrand will then be smooth each time we invoke the Romberg algorithm):
  unsigned nbins = m_density.size()-1;
  for ( unsigned ibin = 0; ibin < nbins; ++ibin ) {
    const double d0=m_density.at(ibin);
    const double d1=m_density.at(ibin+1);
    const double e0=m_emin + m_binwidth*ibin;
    const double e1= (ibin+1==nbins ? m_emax : m_emin + m_binwidth*(ibin+1));
    //In this bin, rho(E) = (E-e0)*(d1-d0)/(e1-e0) + d0 = A*E+B,
    //with A=(d1-d0)/(e1-e0), B=d0-e0 *A
    const double A = (d1-d0)*m_invbinwidth;
    const double B = d0-e0*A;
    //NB: Reduced from integrateRomberg33 to integrateRomberg17 since it showed
    //up in profiling:
    double bincontrib = integrateRomberg17([&f,A,B](double e){return f(e)*(A*e+B);}, e0, e1);
    sum.add(bincontrib);
  }
}

NC::VDOSEval::~VDOSEval() = default;

NC::VDOSEval::VDOSEval(const VDOSData& vd)
  : m_density(vd.vdos_density()),
    m_emin(vd.vdos_egrid().first),
    m_emax(vd.vdos_egrid().second),
    m_kT(constant_boltzmann*vd.temperature().get()),
    m_temperature( DoValidate, vd.temperature()),
    m_elementMassAMU(vd.elementMassAMU())
{
  if ( s_verbose_vdoseval )
    NCRYSTAL_MSG("VDOSEval constructed ("<<m_density.size()
                 <<" density pts on egrid spanning ["<<fmt(m_emin,"%.14g")
                 <<", "<<fmt(m_emax,"%.14g")<<"])");

  nc_assert( m_elementMassAMU.dbl()>0.5 && m_elementMassAMU.dbl()<2000.0 );
  nc_assert_always(m_density.size()<static_cast<std::size_t>(std::numeric_limits<int>::max()-2));
  nc_assert(m_emin>=0.0&&m_emax>m_emin);

  if (m_emin<1e-5) {
    //some integrands will diverge for E->0 if m_emin is too small (below m_emin
    //we approximate the function as k*E^2, which can handle the
    //divergencies). We could of course silently discard the first bin and give
    //a WARNING, instead of throwing an exception. Note that the specific value
    //of the cutoff, 1e-5, is also used elsewhere in the class so it should not
    //just be changed here.
    NCRYSTAL_THROW(BadInput,"VDOS energy range should not be specified for values less than 1e-5eV = 0.01meV");
  }

  double emax_corrected = checkIsRegularVDOSGrid( PairDD(m_emin,m_emax), m_density );
  if ( !emax_corrected )
    NCRYSTAL_THROW(BadInput,"Received non-regularised VDOS. The VDOSEval class expects regularised"
                   " equidistant grid which can be extended downwards and exactly coincide with 0.");

  if ( s_verbose_vdoseval && emax_corrected != m_emax ) {
    NCRYSTAL_MSG("VDOSEval Correcting emax slightly for completely regular"
                 " grid: " << m_emax << " -> "<< emax_corrected << " (relative"
                 " change: " << ((emax_corrected-m_emax)/m_emax) << ")");
    m_emax = emax_corrected;
  }
  const double binwidth = ( m_emax - m_emin ) / ( m_density.size() - 1 );
  const unsigned long nbins_below_emin = static_cast<unsigned long>( m_emin/binwidth  + 0.5 );
  nc_assert( floateq(nbins_below_emin*binwidth,m_emin) );//regularised and corrected emax, so should be ~exact.

  m_nptsExtended = m_density.size() + nbins_below_emin;// [0, ..., emin, ..., emax]

  //below emin, f(e) = k*e^2, with continuity condition at e=emin:
  const double eminsq = m_emin*m_emin;
  nc_assert(eminsq>0.0);
  m_k = m_density.front()/eminsq;

  m_binwidth = binwidth;
  nc_assert_always(m_binwidth>0.0);
  m_invbinwidth = 1.0/m_binwidth;

  //Normalise by first calculating the integral. The integration range
  //(0,m_emin) where the density fct is parabolic is handled analytically:
  StableSum sum_integral;
  constexpr double onethird = 1.0 / 3.0;
  sum_integral.add( onethird * m_density.front() * m_emin );
  integrateBinsWithFunction( [](double){return 1.0;}, sum_integral );
  m_originalIntegral = sum_integral.sum();

  nc_assert_always(m_originalIntegral>0.0);
  const double scalefact = 1.0 / m_originalIntegral;
  for (auto& e : m_density)
    e *= scalefact;
  m_k *= scalefact;
}

double NC::VDOSEval::eval(double energy) const
{
  nc_assert(energy>=0.0);
  if (energy<=m_emin)
    return m_k * energy * energy;
  double relpos = ncclamp((energy-m_emin)*m_invbinwidth, -0.5, m_density.size()+0.5);
  int ibin = static_cast<int>(relpos);
  nc_assert(ibin>=0);
  if (ibin>=static_cast<int>(m_density.size()-1))
    return 0.0;
  relpos = ncclamp(relpos - ibin, 0.0, 1.0);
  return (1.0-relpos)*m_density.at(ibin)+relpos*m_density.at(ibin+1);
}

NC::VDOSEval::GridInfo NC::VDOSEval::getGridInfo() const
{
  nc_assert(std::numeric_limits<unsigned>::max()>m_density.size());
  return { m_emin, m_emax, static_cast<unsigned>(m_density.size()), m_nptsExtended };
}

double NC::VDOSEval::calcGamma0() const
{
  //Evaluate Sjolander1958 eq. II.3 with t=0 (NB: Sjolander is missing a factor
  //of emax, since he uses unit-less energies).
  //
  // We integrate the density with the function:
  //
  // f(E) = coth(E/2kT)/E = 1/E*tanh(E/2kT) = (1/2kt) / ( u*tanh(u) )
  //
  // The parabolic part of the integral, over [0,m_emin] (cf. comments in
  // integrateBinsWithFunction) thus becomes:
  //
  //   m_k* (2kt)^2 * integral_0^(m_emin/2kT){ u/(tanh(u))du }
  //
  // We use the safe_xcothx function to evaluate the integrand (needs a Taylor
  // expansion near u~=0).

  //First the parabolic part [0,emin]:
  const double twokT = 2.0*m_kT;
  const double inv2kT = 1.0 / twokT;
  StableSum sum;
  sum.add( m_k * ncsquare( twokT ) * integrateRomberg33( safe_xcothx, 0.0, m_emin * inv2kT ) );

  //And now the binned part [emin,emax]:
  auto f_binnedpart =  [ inv2kT ]( double e )
  {
    //NB: No Taylor expansion is possible near E/2kT~=0 here, but the precision
    //of 1/tanh(x) goes down for very small x, so for very large T and small
    //m_emin, we might introduce some imprecision here. For instance, if T=1e6K
    //(the most extreme we claim to support) and emin = 1e-5eV (the smallest we
    //claim to support), we get a precision of around 5e-11, which is likely
    //good enough. For reference, if we were to allow emin=1e-7eV we would get a
    //precision of just 3e-8 and with emin=1e-10eV it would become
    //4e-5. However, this is only at the extreme T=1e6K, at T=1000K and
    //emin=1e-10eV we would get a precision of 3e-8. So it would most likely be
    //OK for this integration to loosen the emin limit from 1e-5 to 1e-10 eV.
    return 1.0 / ( e * std::tanh( e * inv2kT ) );
  };
  integrateBinsWithFunction( f_binnedpart, sum );
  return m_emax * sum.sum();
}

double NC::VDOSEval::calcEffectiveTemperature() const
{
  //Defining the effective temperature Teff via
  //k*Teff=mean-energy-per-phonon-state, we can integrate over the energy per
  //phonon-state (Sjolander I.6) which is (1/2)*E*coth(E/2kT) and get:
  //
  // k*Teff = integral( rho(E) * (1/2)*E*coth(E/2kT) )
  //
  //  <=>
  //
  // Teff = (1/2k) * integral( rho(E) * E / tanh(E/2kT) )
  //
  //Two notes:
  //
  //   1) The factor of 3 in Sjolander's eq. I.7 is due to the three dimensions.
  //   2) The average kinetic energy of a particle in a free-gas Maxwell
  //      distribution is 3kT/2, or kT/2 per degree of freedom. Thus, which in
  //      the absence of internal degrees of freedom, similar logic as above
  //      applied to a free gas would yield an effective temperature of k*Teff =
  //      kT/2, or half of the actual temperature.
  //
  //   => All in all, one might wonder if we have a factor of 2 or 3 wrong in
  //   our result. However, plotting the results for a few materials shows that
  //   we get Teff -> T at large temperatures, making it unlikely (todo: it
  //   would of course be nice to clarify this further, but doesn't seem
  //   urgent).

  // We implement the integration in a similar manner as to what is done in
  // ::calcGamma0, with the difference of course being the function f(E) which
  // is instead f(E) = (1/2k) * E / tanh(E/2kT) = T * u / tanh(u). Correspondingly
  // the integrand of the parabolic part of the integral becomes:
  //
  //   parabolic part :  m_k* (2kT)^3 * integral_0^(ulim){ u^3/(tanh(u) }du
  //
  // We use the safe_x3cothx function to evaluate the integrand (needs a Taylor
  // expansion near u~=0).

  //First the parabolic part [0,emin] (leaving out the overall factor of T
  //until the return statement):
  const double twokT = 2.0 * m_kT;
  const double inv2kT = 1.0 / twokT;
  StableSum sum;
  sum.add( m_k * nccube( twokT ) * integrateRomberg33( safe_x3cothx, 0.0, m_emin * inv2kT ) );

  //Binned part [emin,emax] (again leaving out a factor of T until the return
  //statement):
  integrateBinsWithFunction( [inv2kT](double E){ return safe_xcothx( E * inv2kT ); }, sum );

  return m_temperature.dbl() * sum.sum();
}

double NC::VDOSEval::getMSD( double gamma0 ) const
{
  //On p326 of Sjolander1958, he gives the formula Q*gamma0=2W, with W the
  //Debye-Waller function 2W=msd*q^2 and (eq II.1) Q=gamma0*hbar^2*q^2/2M (NB:
  //Sjolander left out a factor of hbar in II.1, which we have inserted here, and
  //Sjolander uses the notation chi^2 instead of q^2).
  //
  //Thus, all in all:
  //
  //  msd = gamma0 * hbar^2 / (2M)
  //
  //Finally, we also need to apply a factor of 1/emax to the resulting msd,
  //since Sjolander works with unit-less energy and length, normalised to
  //respectively emax or 1/emax.

  constexpr double convfact = 0.5 * ncsquare(constant_hbar) * ncsquare(constant_c) / constant_dalton2eVc2;
  return convfact * gamma0 / ( m_elementMassAMU.dbl() * m_emax );
}

NC::PairDD NC::VDOSEval::evalG1AsymmetricAtEPair( double energy, double gamma0 ) const
{
  nc_assert( ! ( energy < 0.0 ) );
  if ( energy < numericallySafeG1SymmetricELimitInUnitsOfKT * m_kT ) {
    double G1sym = evalG1Symmetric( energy, gamma0 );
    if (!G1sym)
      return { 0.0, 0.0 };
    const double dbfact = std::exp( energy / (2*m_kT) );
    nc_assert( dbfact > 0.0 );
    nc_assert( G1sym > 0.0 );
    return { G1sym*dbfact, G1sym/dbfact };
  } else {
    nc_assert( energy * gamma0 > 0.0 );
    const double kkk = eval( energy ) * m_emax / ( energy * gamma0 );
    if ( !kkk )
      return { 0.0, 0.0 };
    return { -kkk / std::expm1( -energy / m_kT ),//-energy (absorb sign on kkk as well)
             kkk / std::expm1( +energy / m_kT ) };
  }
}

double NC::VDOSEval::evalG1Asymmetric( double energy, double gamma0 ) const
{
  const double absE = ncabs( energy );
  if ( absE > numericallySafeG1SymmetricELimitInUnitsOfKT * m_kT ) {
    //evaluating at point where the detailed balance factor is more than
    //exp(+-100), leading to numerically unstable cancellations between the
    //detailed balance factor and G1Sym. Thus, we use separate analytical
    //formulas for +energy and -energy. In principle we could always use the
    //formula below, except for the fact that the G1
    nc_assert( absE * gamma0 > 0.0 );
    const double kkk = eval( absE ) * m_emax / ( energy * gamma0 );
    if ( !kkk )
      return 0.0;
    return kkk / std::expm1( energy / m_kT );
  }
  double G1sym = evalG1Symmetric(absE,gamma0);
  return G1sym ? G1sym * std::exp( -energy / (2*kT()) ) : 0.0;
}

double NC::VDOSEval::evalG1Symmetric( double energy, double gamma0 ) const
{
  if (energy<0)
    return evalG1Symmetric(-energy,gamma0);

  //Sjolander1958 II.29 (but must add another factor of emax, since Sjolander
  //use unit-less energies):
  //
  //G1 = f(E)/(E*2*gamma0*sinh(E/2kT))

  const double twokT=2*m_kT;
  const double u = energy / twokT;
  if (energy<=m_emin) {
    //Here f(E) = m_k*E^2, so with u=E/2kT: G1 = (m_k*2kT/2*gamma0) * u / (sinh(u))
    return ( m_k * m_kT * m_emax / gamma0 ) * safe_xdivsinhx( u );
  } else {
    return eval(energy) * m_emax / ( energy * 2.0 * gamma0 * std::sinh(u) );
  }
}

namespace NCRYSTAL_NAMESPACE {
  namespace {
    bool isLinearlySpacedGrid(const VectD& grid, double tolerance ) {
      nc_assert(nc_is_grid(grid));
      nc_assert(grid.size()>=2);
      if (grid.size()==2)
        return true;
      double binwidth = ( grid.back()-grid.front() ) / ( grid.size() - 1.0 );
      double eps = tolerance*binwidth;
      for ( auto e : enumerate(grid) ) {
        if ( ncabs( (grid.front() + e.idx * binwidth) - e.val ) > eps )
          return false;
      }
      return true;
    }
  }
}

double NC::checkIsRegularVDOSGrid( const PairDD& egrid, const VectD& density, double tolerance )
{
  nc_assert_always(egrid.first>=1e-5);//checking 1e-5 threshold in VDOSEval constructor

  const unsigned long npts = density.size();
  nc_assert_always(npts>=2);
  double emax = egrid.second;
  double emin = egrid.first;
  nc_assert_always(emax>emin);
  const double binwidth_approx = (emax-emin)/(npts-1);
  const double nbins_below_emin_flt = emin/binwidth_approx;
  if ( nbins_below_emin_flt < 0.99 || ncabs( nbins_below_emin_flt - std::round(nbins_below_emin_flt) ) > tolerance )
    return 0.0;

  const unsigned long nbins_below_emin = static_cast<unsigned long>(nbins_below_emin_flt+0.5);
  const double  binwidth = emin / nbins_below_emin;
  const double emax_corrected = emin + (npts-1)*binwidth;
  nc_assert( emax_corrected>0.0 );
  return emax_corrected;
}

double NC::checkIsRegularVDOSGrid( const VectD& egrid, const VectD& density, double tolerance )
{
  nc_assert( nc_is_grid(egrid) );
  if ( egrid.size()!=2 && egrid.size() != density.size() )
    NCRYSTAL_THROW(BadInput,"VDOS energy grid vector must be 2 or have same size as density vector");
  if ( !isLinearlySpacedGrid(egrid,tolerance) )
    return 0.0;
  return checkIsRegularVDOSGrid( PairDD(egrid.front(),egrid.back()), density, tolerance );
}

std::pair<NC::VectD,NC::VectD> NC::regulariseVDOSGrid( const VectD& orig_egrid, const VectD& orig_density )
{
  static bool s_verbose_vdosregul = ncgetenv_bool("DEBUG_VDOSREGULARISATION");
  const bool extra_verbose = s_verbose_vdosregul;
  nc_assert_always( orig_density.size() > 2) ;
  nc_assert_always( orig_density.size() < 4000000000) ;
  nc_assert_always( orig_egrid.size()==2
                    || orig_egrid.size() == orig_density.size() );
  nc_assert_always( nc_is_grid(orig_egrid) );
  nc_assert_always( orig_egrid.front() >= 0.0 );

  if ( orig_egrid.front() < 1e-5 )
    NCRYSTAL_THROW(BadInput,"VDOS energy range can not be specified"
                   " for values less than 1e-5eV = 0.01meV");

  if ( extra_verbose )
    NCRYSTAL_MSG("Called regulariseVDOSGrid(["<<orig_egrid.front()<<",..,"
                 <<orig_egrid.back()<<"], "<<orig_density.size()<<" density pts");

  const double tolerance = 1e-6;//NB: Should be same as checkIsRegularVDOSGrid
                                //    default value!
  double emax_corrected = NC::checkIsRegularVDOSGrid( orig_egrid, orig_density, tolerance );
  if ( emax_corrected ) {
    //Is already OK within tolerance! Return regularised egrid, with potential
    //slight correction to emax to correct for numerical imprecision within
    //the allowed tolerance:
    if ( s_verbose_vdoseval ) {
      std::ostringstream msg;
      msg<<"regulariseVDOSGrid Grid was already regular within tolerance of "
         <<tolerance;
      if ( orig_egrid.back() != emax_corrected )
        msg<<" (corrected emax slightly "<<orig_egrid.back()<<" -> "
           <<emax_corrected<<", a relative change of "
           <<(emax_corrected/orig_egrid.back()-1.0)<<")";
      NCRYSTAL_MSG(msg.str());
    }
    return { VectD({orig_egrid.front(),emax_corrected}), orig_density };
  }

  //Ok, we need to regularise! We do this enforcing a uniform binwidth and
  //requiring that [0,emin] can be divided into a whole number of bins,
  //k. This will result in m whole bins in the interval [emin,emax], plus a
  //remainder, epsilon. Starting from the k-value which gives 2000 bins in
  //[emin,emax], we try increasing k-values and look for one which minimises
  //epsilon (to minimise re-parameterisation artifacts). The reason for
  //requiring at least 2000 bins, is to try to minimise the chance of
  //re-parameterisation artifacts)

  const double emin = orig_egrid.front();
  nc_assert(emin>0.0 );
  const double oldEmax = orig_egrid.back();
  nc_assert_always(oldEmax > emin );
  const double oldEmaxMinusEmin = oldEmax-emin;
  nc_assert_always(oldEmaxMinusEmin > 0.0 );
  const double oldEmaxMinusEminDivEmin = oldEmaxMinusEmin/emin;
  nc_assert_always(oldEmaxMinusEminDivEmin > 0.0 );

  const double kBegin = 2000.0/oldEmaxMinusEminDivEmin;//at least 2000 bins
                                                       //inside [emin,emax]
  double k = ncmax(1.0,std::round(kBegin));//k is double, to avoid conversions
                                           //below.
  PairDD best = { kInfinity, 0.0 };

  while ( true ) {
    nc_assert(k >= 1.0 );
    const double oldEmaxMinusEminDivEmin_mult_k = oldEmaxMinusEminDivEmin * k;
    const double m = std::floor(oldEmaxMinusEminDivEmin_mult_k);
    if ( m < 1.0 ) {
      if ( extra_verbose )
        NCRYSTAL_MSG("regulariseVDOSGrid skip k since m<1");
      k += 1.0;
      continue;
    }
    const double binwidth = emin / k;
    double eps = oldEmaxMinusEmin - (m*binwidth);
    nc_assert( eps >= -oldEmaxMinusEmin*1e-10 );
    eps = ncmax(0.0, eps );

    const double kbest = best.second;

    if ( extra_verbose )
      NCRYSTAL_MSG("regulariseVDOSGrid trying k="<<k<<" (m="
                   <<m<<", eps="<<eps<<")");

    //We want to check if "eps < best.first", but the calculation of eps
    //involves a subtraction which is numerically unstable. Thus we do it like
    //this instead:
    if ( eps == 0.0 || kbest == 0.0
         || ( kbest*std::floor(oldEmaxMinusEminDivEmin * k)
              >  k*std::floor(oldEmaxMinusEminDivEmin * kbest) ) ) {
      //Beats it
      if ( extra_verbose )
        NCRYSTAL_MSG("regulariseVDOSGrid NEW BEST k="<<k
                     <<" (reduces epsilon by factor "<<eps/best.first
                     <<" = 1-"<<(1.0-eps/best.first)<<")");
      best = { eps, k };
    }

    //Check if best result so far is acceptable. We lower our requirement as we
    //go along and get more and more desperate:
    double tol = 1e-6;
    if ( m > 5000 ) {
      tol = 1e-5;
      if ( m > 10000 ) {
        tol = 1e-4;
        if ( m > 15000 )
          tol = m > 19000 ? 1e-2 : 1e-3;
      }
    }
    if ( extra_verbose )
      NCRYSTAL_MSG("regulariseVDOSGrid Checking best tol (from k="
                   <<best.second<<"):  "<<best.first<<" versus "
                   <<oldEmaxMinusEmin*tol<<" = "<<oldEmaxMinusEmin<<" * "<<tol);
    if ( best.first < oldEmaxMinusEmin*tol )
      break;
    if ( m >= 20000 )
      NCRYSTAL_THROW(BadInput,"Could not regularise input energy grid."
                     " Are the energy ranges highly unusual?");
    k += 1.0;
  }

  double new_binwidth = emin / best.second;
  double mm = std::floor( oldEmaxMinusEminDivEmin * best.second );
  nc_assert_always ( mm <= 20000 );
  unsigned new_npts = static_cast<unsigned>( mm + 0.5 ) + 1;
  double new_emax = emin + new_binwidth * ( new_npts - 1 );

  if ( new_emax < oldEmax ) {
    //Add one extra point to make sure new range encompasses old range.
    ++new_npts;
    new_emax = emin + new_binwidth * ( new_npts - 1 );
  }
  nc_assert( new_emax >= oldEmax);
  if ( extra_verbose )
    NCRYSTAL_MSG("regulariseVDOSGrid new binwidth="
                 <<fmt(new_binwidth,"%.14g")
                 <<", emin="<<fmt(emin,"%.14g")
                 <<", emax="<<fmt(new_emax,"%.14g")
                 <<", npts="<<new_npts);

  //Ok, we need to regularise! We do this by putting a high number of points
  //linearly spaced between [emin,emax-eps] where eps is determined as the
  //smallest possible number which will make the result regular. The
  //corresponding new density values are found by interpolation in the
  //original values. Note that npts=nbins+1.

  //Calculate new density values by interpolation in old ones:

  VectD newEgrid({emin,new_emax});
  VectD newDensity;
  newDensity.reserve(new_npts);
  double orig_egrid_inverse_binwidth = 0.0;//0.0 means N/A
  if ( orig_egrid.size() == 2 )//todo: also use isLinearlySpacedGrid??
    orig_egrid_inverse_binwidth = ( ( orig_density.size() - 1.0 )
                                    / (orig_egrid.back()-orig_egrid.front()) );

  VectD orig_egrid_expanded = ( orig_egrid.size() == 2
                                ? linspace(orig_egrid.front(),
                                           orig_egrid.back(),
                                           orig_density.size())
                                : orig_egrid );
  nc_assert_always( orig_egrid_expanded.size() == orig_density.size() );
  auto it = orig_egrid_expanded.begin();
  auto itBegin = orig_egrid_expanded.begin();
  auto itLast = std::prev(orig_egrid_expanded.end());

  //Loop like this rather than calling linspace(..), since new_binwidth is not
  //exactly the same as (emax-emin)/nbins-1 due to FP inaccuracies:
  for ( auto ieval : ncrange(new_npts) ) {
    const double eval = ( ieval + 1 == new_npts
                          ? new_emax
                          : emin + new_binwidth * ieval++ );
    //Increment position in old grid if needed:
    while ( it != itLast && eval >= *std::next(it) )
      ++it;
    if ( eval == *it ) {
      //Exactly at old grid point, just transfer value:
      newDensity.push_back( vectAt(orig_density,std::distance(itBegin,it)) );
      continue;
    }
    if ( it == itLast ) {
      //eval is beyond old grid (can only happen for the very last point in the
      //new grid):
      nc_assert( eval == newEgrid.back() );
      nc_assert( eval > *itLast );
      newDensity.push_back(  eval > *it ? 0.0 : orig_density.back() );
    } else {
      //Interpolate:
      auto itNext = std::next(it);
      double y0 = vectAt(orig_density,std::distance(itBegin,it));
      double y1 = vectAt(orig_density,std::distance(itBegin,itNext));
      double x0 = *it;
      double x1 = *itNext;
#if 0
      //Fast, but not quite stable near r=1:
      double r = ( eval - x0 ) / ( x1 - x0 );
      nc_assert(r>=0.0 && r<=1.0);
      newDensity.push_back( y0 * (1.0 - r ) + y1 * r );
#else
      //Stable:
      StableSum sum;
      sum.add(y0*x1);
      sum.add(- y0*eval);
      sum.add(eval*y1);
      sum.add(- x0*y1);
      if ( orig_egrid_inverse_binwidth )
        newDensity.push_back( sum.sum() * orig_egrid_inverse_binwidth );
      else
        newDensity.push_back( sum.sum()/(x1-x0) );
#endif
    }
  }
  nc_assert( newDensity.size() == new_npts );

  if ( s_verbose_vdoseval )
    NCRYSTAL_MSG("regulariseVDOSGrid Grid was regularised using "
                 << newDensity.size()
                 << " equidistant points on interval [" << newEgrid.front()
                 << ", " << newEgrid.back() << "]");

  nc_assert(newEgrid.size()==2);
  return { newEgrid, newDensity };
}

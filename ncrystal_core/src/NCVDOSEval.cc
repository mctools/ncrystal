////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2022 NCrystal developers                                   //
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

#include "NCrystal/internal/NCVDOSEval.hh"
#include "NCrystal/internal/NCMath.hh"
#include "NCrystal/internal/NCIter.hh"
#include <iostream>
namespace NC=NCrystal;

namespace NCrystal {
  static std::atomic<bool> s_verbose_vdoseval( getenv("NCRYSTAL_DEBUG_PHONON")!=nullptr );
}

void NC::VDOSEval::enableVerboseOutput(bool status)
{
  s_verbose_vdoseval = status;
}

bool NC::VDOSEval::verboseOutputEnabled()
{
  return s_verbose_vdoseval;
}


NC::VDOSEval::~VDOSEval() = default;

NC::VDOSEval::VDOSEval(const VDOSData& vd)
  : m_density(vd.vdos_density()),
    m_emin(vd.vdos_egrid().first),
    m_emax(vd.vdos_egrid().second),
    m_kT(constant_boltzmann*vd.temperature().get()),
    m_temperature(vd.temperature()),
    m_elementMassAMU(vd.elementMassAMU())
{
  if ( s_verbose_vdoseval )
    std::cout << "NCrystal::VDOSEval constructed ("<<m_density.size()<<" density pts on egrid spanning ["<<m_emin<<", "<<m_emax<<"]"<<std::endl;

  nc_assert( m_elementMassAMU.dbl()>0.5 && m_elementMassAMU.dbl()<2000.0 );
  nc_assert( m_temperature.get()>=1.0&&m_temperature.get()<1e5 );
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
    std::cout << "NCrystal::VDOSEval Correcting emax slightly for completely regular grid: " << m_emax << " -> "
              << emax_corrected << " (relative change: " << ((emax_corrected-m_emax)/m_emax) << ")" << std::endl;
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

  //Normalise:
  m_originalIntegral = integrateWithFunction([](double){return 1.0;},
                                             [](double e){return e*e;});
  nc_assert_always(m_originalIntegral>0.0);
  double scalefact = 1.0/m_originalIntegral;
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

  //We integrate the density with the function coth(E/2kT)/E = 1/E*tanh(E/2kT)
  //
  //For low values of u=E/2kT, we should expand E^2*coth(E/2kT)/E=2kT*u*coth(u)
  //in u. The expansion will be called only when E<1e-5, so u<1e-5/2kT =
  //0.058/T[K]. To support T=1K, we thus need to ensure that the expansion is
  //precise (at double precision level) to u=0.06. This is satisfied using an
  //8th order taylor expansion.

  const double twokT = 2.0*m_kT;
  const double inv2kT = 1.0 / twokT;
  constexpr double c2 = 1./3.;
  constexpr double c4 = -1./45.;
  constexpr double c6 = 2./945.;
  constexpr double c8 = -1./4725.;
  return integrateWithFunction([inv2kT](double e){ return 1.0/(e*std::tanh(e*inv2kT)); },
                               [c2,c4,c6,c8,twokT,inv2kT](double e)
                               {
                                 (void)c2;//avoid annoying misleading clang warning
                                 (void)c4;//avoid annoying misleading clang warning
                                 (void)c6;//avoid annoying misleading clang warning
                                 (void)c8;//avoid annoying misleading clang warning
                                 const double u=e*inv2kT;
                                 const double u2=u*u;
                                 return (1.0+u2*(c2+u2*(c4+u2*(c6+u2*c8))))*twokT;
                               })*m_emax;

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
  // Teff = (1/2K) * integral( rho(E) * E / tanh(E/2kT) )
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

  const double twokT=2.0*m_kT;
  const double inv2kT = 1.0 / twokT;

  //For Taylor expansion: F(E)*E^2 = E^3 / tanh(E/2kT) = (2kT)^3 * u^3/tanh(u),
  //u=E/2kT. We keep 10th order of Taylor expansion to be valid beyond u=0.06,
  //for the same reasons given in the calcGamma0 function (i.e. to work up to
  //E=1e-5eV and all the way down to T=1K).

  const double kkk = twokT*twokT*twokT;
  const double c4 = 1./3.;
  const double c6 = -1./45.;
  const double c8 = 2./945.;
  const double c10 = -1./4725.;
  constexpr double ccc = 0.5/constant_boltzmann;

  return ccc*integrateWithFunction( [inv2kT](double e){ return e / std::tanh(e*inv2kT); },
                                    [c4,c6,c8,c10,kkk,inv2kT](double e)
                                    {
                                      const double u = inv2kT*e;
                                      const double u2 = u*u;
                                      return kkk * u2 * (1.0+u2*(c4+u2*(c6+u2*(c8+u2*c10))));
                                    });
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

  constexpr double convfact = 0.5*constant_hbar*constant_hbar*constant_c*constant_c / constant_dalton2eVc2;
  return convfact * gamma0 / ( m_elementMassAMU.dbl() * m_emax );
}

double NC::VDOSEval::evalG1Asymmetric( double energy, double gamma0 ) const
{
  double G1sym = evalG1Symmetric(ncabs(energy),gamma0);
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
    const double fact = m_k*m_kT*m_emax / gamma0;
    if (u<0.07) {
      //taylor expand u/sinh(u)
      constexpr double c2 = -1./6.;
      constexpr double c4 = 7./360.;
      constexpr double c6 = -31./15120.;
      constexpr double c8 = 127./604800.;
      const double u2 = u*u;
      return fact * ( 1.0+u2*(c2+u2*(c4+u2*(c6+u2*c8))) );
    } else {
      return fact * u / std::sinh(u);
    }
  } else {
    return eval(energy) * m_emax / ( energy * 2.0 * gamma0 * std::sinh(u) );
  }
}

template <class Fct, class FctEsqTaylor>
double NC::VDOSEval::integrateWithFunction(Fct f, FctEsqTaylor fesqtaylor ) const
{
  //With VDOS=rho(E), integrate rho(E)*f(E) over [0,infinity].
  //
  //F might diverge for small E, but since for small E we have rho(E)=k*E^2, the
  //combined integrand might be stable (and is for the present cases). For that
  //reason, the caller should provide fesqtaylor with an approximation for
  //f(E)*E^2 at small E (i.e. a Taylor expansion of F(E)*E^2). It should contain
  //enough terms to be valid until E=1e-5.

  constexpr double etiny = 0.9e-5;
  nc_assert( m_emin>=1e-5 && etiny < m_emin*0.91 );

  //Sanity check that provided functions are compatible:
  nc_assert( floateq( f(etiny)*etiny*etiny, fesqtaylor(etiny),1e-14,1e-99) );

  //NB: If profiling shows that we spend a significant portion of time here, we
  //could use integrateRomberg17 instead.

  StableSum sum;

  //Contribution in [0,etiny]:
  sum.add(m_k*integrateRomberg33(fesqtaylor, 0.0, etiny));

  //Contribution in [etiny,m_emin]:
  sum.add(m_k*integrateRomberg33([&f](double e){return f(e)*e*e;}, etiny, m_emin));

  //Contribution in [m_emin,m_emax], considering individual bins (since the
  //integrand will then be smooth each time we invoke the integration algorithm):
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
    double bincontrib = integrateRomberg33([&f,A,B](double e){return f(e)*(A*e+B);}, e0, e1);
    sum.add(bincontrib);
  }

  return sum.sum();
}


namespace NCrystal {
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
  nc_assert_always( orig_density.size() > 2) ;
  nc_assert_always( orig_density.size() < 4000000000) ;
  nc_assert_always( orig_egrid.size()==2 || orig_egrid.size() == orig_density.size() );
  nc_assert_always( nc_is_grid(orig_egrid) );
  nc_assert_always( orig_egrid.front() >= 0.0 );

  if ( orig_egrid.front() < 1e-5 )
    NCRYSTAL_THROW(BadInput,"VDOS energy range can not be specified for values less than 1e-5eV = 0.01meV");

  const double tolerance = 1e-6;
  double emax_corrected = NC::checkIsRegularVDOSGrid( orig_egrid, orig_density );
  if ( emax_corrected ) {
      //Is already OK within tolerance! Return regularised egrid, with potential
      //slight correction to emax to correct for numerical imprecision within
      //the allowed tolerance:
    if ( s_verbose_vdoseval ) {
      std::cout<<"NCrystal::regulariseVDOSGrid Grid was already regular within tolerance of "<<tolerance;
      if ( orig_egrid.back() != emax_corrected )
        std::cout<<" (corrected emax slightly "<<orig_egrid.back()<<" -> "<<emax_corrected<<", a relative change of "<<(emax_corrected/orig_egrid.back()-1.0)<<")";
      std::cout<<std::endl;
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
  const double oldEmax = orig_egrid.back();
  const double oldEmaxMinusEmin = oldEmax-emin;
  const double oldEmaxMinusEminDivEmin = oldEmaxMinusEmin/emin;

  double kBegin = 2000.0/oldEmaxMinusEminDivEmin;//at least 2000 bins inside [emin,emax]
  double k = ncmax(1.0,std::round(kBegin));//k is double, to avoid conversions below.
  PairDD best = { kInfinity, 0.0 };
  while ( true ) {
    double m = std::floor(oldEmaxMinusEminDivEmin * k);
    double binwidth = emin / k;
    double eps = oldEmaxMinusEmin - (m*binwidth);
    nc_assert_always(eps>=0.0);
    if ( eps < best.first )
      best = { eps, k };
    //Check if best result so far is acceptable. We lower our requirement as we go along:
    double tol = 1e-6;
    if ( m > 5000 ) {
      tol = 1e-5;
      if ( m > 10000 ) {
        tol = 1e-4;
        if ( m > 15000 )
          tol = m > 19000 ? 1e-2 : 1e-3;
      }
    }
    if ( best.first < oldEmaxMinusEmin*tol )
      break;
    if ( m >= 20000 )
      NCRYSTAL_THROW(BadInput,"Could not regularise input energy grid. Are the energy ranges highly unusual?");
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
  nc_assert ( new_emax >= oldEmax);

  //Ok, we need to regularise! We do this by putting a high number of points
  //linearly spaced between [emin,emax-eps] where eps is determined as the
  //smallest possible number which will make the result regular. The
  //corresponding new density values are found by interpolation in the
  //original values. Note that npts=nbins+1.

  //Calculate new density values by interpolation in old ones:

  VectD newEgrid({emin,new_emax});
  VectD newDensity;
  newDensity.reserve(new_npts);
  VectD orig_egrid_expanded = ( orig_egrid.size() == 2
                                ? linspace(orig_egrid.front(),orig_egrid.back(),orig_density.size())
                                : orig_egrid );
  nc_assert_always( orig_egrid_expanded.size() == orig_density.size() );
  auto it = orig_egrid_expanded.begin();
  auto itBegin = orig_egrid_expanded.begin();
  auto itLast = std::prev(orig_egrid_expanded.end());

  for ( auto eval : linspace(newEgrid.front(),newEgrid.back(),new_npts) ) {
    //Increment position in old grid if needed:
    while ( it != itLast && eval >= *std::next(it) )
      ++it;
    if ( eval == *it ) {
      //Exactly at old grid point, just transfer value:
      newDensity.push_back( vectAt(orig_density,std::distance(itBegin,it)) );
      continue;
    }
    if ( it == itLast ) {
      //eval is beyond old grid (can only happen for the very last point in the new grid):
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
      double r = ( eval - x0 ) / ( x1 - x0 );
      newDensity.push_back( y0 * (1.0 - r ) + y1 * r );
    }
  }
  nc_assert( newDensity.size() == new_npts );

  if ( s_verbose_vdoseval )
    std::cout << "NCrystal::regulariseVDOSGrid Grid was regularised using " << newDensity.size()
              << " equidistant points on interval [" << newEgrid.front() << ", " << newEgrid.back() << "]" << std::endl;

  nc_assert(newEgrid.size()==2);
  return { newEgrid, newDensity };
}

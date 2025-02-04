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

#include "NCrystal/internal/phys_utils/NCGaussOnSphere.hh"
#include "NCrystal/internal/utils/NCRomberg.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCrystal/internal/utils/NCString.hh"
namespace NC = NCrystal;

namespace NCRYSTAL_NAMESPACE {
  class CosExpansionRadiusFct : public Fct1D {
    double m_prec;
  public:
    CosExpansionRadiusFct(double prec) : m_prec(prec) {};
    virtual ~CosExpansionRadiusFct(){}
    virtual double eval(double x) const
    {
      return (1.0-(1.0-0.5*x*x)/cos_mpi2pi2(x)) - m_prec;
    }
  };
  double gos_cosexpansionradius(double target_precision) {
    //Returns largest x for which (cos(x)-(1-0.5*x^2))/cos(x) is smaller than the
    //target_precision (but clips the target precision to at worst 0.999999
    //corresponding to a |x| smaller than roughly 81degrees.
    nc_assert_always(target_precision>0);
    CosExpansionRadiusFct f(ncmin(0.999999,target_precision));
    return findRoot(&f,0.0,0.999999*kPiHalf,1e-13);
  }
}

namespace NCRYSTAL_NAMESPACE {

  class GOSNormInt : public Romberg {
  public:
    GOSNormInt(double mos): Romberg(), m_c(-0.5/(mos*mos)) {}
    virtual ~GOSNormInt(){}
    double evalFunc(double x) const { return ncmax(0.0, std::sin( x ) * std::exp( m_c * x * x)); }
    virtual bool accept(unsigned level, double prev_estimate, double estimate,double,double) const
    {
      return level>12 || ncabs(prev_estimate-estimate)<1e-12*ncabs(estimate);
    }
  private:
    const double m_c;
  };

  class GOSCircleInt : public Romberg {
  public:
    GOSCircleInt(const GaussOnSphere*gos, double sinalpha_times_singamma, double cosalpha_times_cosgamma, double acc, bool statcollect)
      : Romberg(),
        m_gos(gos),
        m_sasg(sinalpha_times_singamma),
        m_cacg(cosalpha_times_cosgamma),
        m_acc(acc),
        m_nevals(statcollect?1:0)
    {}
    virtual ~GOSCircleInt(){}

    virtual double evalFunc(double) const {
      //Must be implemented, but should never be called since we provide
      //evalFuncMany+evalFuncManySum.
      nc_assert_always(false);
    }

    virtual void evalFuncMany(double* fvals, unsigned n, double offset, double delta) const
    {
      if (m_nevals)
        m_nevals += n;
      nc_assert(offset>=0&&offset<kPi*1.00001);
      nc_assert(delta>0&&delta*(n-1)<=kPi*1.00001);
      CosSinGridGen grid(n,offset,delta);
      unsigned i=0;
      do {
        double cb = m_sasg * grid.current_cosval() + m_cacg;
        nc_assert(NC::ncabs(cb)<1.000000001);
        fvals[i++] = m_gos->evalCosXInRange(cb);
      } while (grid.step());
    }

    virtual double evalFuncManySum(unsigned n, double offset, double delta) const
    {
      if (m_nevals)
        m_nevals += n;
      nc_assert(offset>=0&&offset<kPi*1.00001);
      nc_assert(delta>0&&delta*n<=kPi*1.00001);
      CosSinGridGen grid(n,offset,delta);
      double sum(0.);
      do {
        double cb = m_sasg * grid.current_cosval() + m_cacg;
        nc_assert(NC::ncabs(cb)<1.000000001);
        sum += m_gos->evalCosXInRange(cb);
      } while (grid.step());
      return sum;
    }

    virtual bool accept(unsigned level, double prev_estimate, double estimate, double a, double b) const
    {
      if (ncabs(prev_estimate-estimate)<=m_acc*ncabs(estimate))
        return true;
      if (level<11)
        return false;

      //Damn. Almost time to emit warning, but if the result is almost at our
      //level, we just accept it quietly if it is still very precise.
      if (m_acc<5e-4&&ncabs(prev_estimate-estimate)<=10.0*m_acc*ncabs(estimate))
        return true;

      static bool first = true;
      if (first) {
        first = false;
        unsigned twotolevelm1 = 1<<(level-1);//2^(level-1)
        NCRYSTAL_WARN("Problems during numerical integration of Gaussian density on sphere. Romberg integration"
                      " did not converge after "<<2*twotolevelm1+1<<" function evaluations (requested acc="<<m_acc<<", got acc="<<
                      ncabs(prev_estimate-estimate)/(ncmax(1e-300,ncabs(estimate)))<<"). Dumping integrand to"
                      " ncrystal_goscircleintegral_fct.txt for debugging. Further warnings"
                      " of this type will not be emitted.");
        writeFctToFile("ncrystal_goscircleintegral_fct.txt", a, b,twotolevelm1);
      }
      return true;
    }

    unsigned nEvals() const { nc_assert(m_nevals>0); return m_nevals-1; }

  private:
    const GaussOnSphere* m_gos;
    const double m_sasg;
    const double m_cacg;
    const double m_acc;
    mutable unsigned m_nevals;//mutable here is MT-safe, since we are only using local GOSCircleInt instances.
  };

  class SLTFct_EvalCosX : public Fct1D {
    double m_norm, m_expval;
  public:
    SLTFct_EvalCosX(double norm, double expval) : m_norm(norm), m_expval(expval) {}
    virtual ~SLTFct_EvalCosX(){}
    virtual double eval(double cosx) const
    {
      nc_assert(ncabs(cosx)<1.000001);
      double x = std::acos(ncmin(1.0,ncmax(-1.0,cosx)));
      return m_norm * std::exp(m_expval*x*x);
    }
  };
  class SLTFct_SofCosD : public Fct1D {
    double m_k, m_expfact, m_tasq;
  public:
    SLTFct_SofCosD(double norm, double expfact, double sigma, double truncangle)
      : m_k(norm*2.506628274631000502415765284811045253006986740609938316629923576*sigma),//2.5066... is sqrt(2pi)
        m_expfact(expfact), m_tasq(truncangle*truncangle) {}
    virtual ~SLTFct_SofCosD(){}
    virtual double eval(double cosd) const
    {
      nc_assert(ncabs(cosd)<1.000001);
      double d = std::acos(ncmin(1.0,ncmax(-1.0,cosd)));
      double d2 = d*d;
      return m_k * std::exp(m_expfact*d2) * std::erf(std::sqrt(ncmax(0.0,-m_expfact*(m_tasq-d2))));
    }
  };



}

double NC::GaussOnSphere::calcNormFactor( double sigma, double trunc_angle )
{
  //Find by integration, but limit to at most 20*sigma for numerical
  //stability (contributions outside will likely not affect the result
  //anyway):
  GOSNormInt normcalc(sigma);
  double integral = k2Pi*normcalc.integrate(0.0,ncmin(trunc_angle,20.0*sigma));
  nc_assert(integral>0.0);
  return 1.0/integral;
}


NC::GaussOnSphere::GaussOnSphere()
  : m_cta(-1.0),
    m_circleint_k1(-1.0),
    m_circleint_k2(-1.0),
    m_norm(-1.0),
    m_expfact(-1.0),
    m_truncangle(-1.0),
    m_sigma(-1.0),
    m_numint_accuracy(-1.0),
    m_prec(-1),
    m_sta(-1.0)
{
  nc_assert(!isValid());
}

NC::GaussOnSphere::GaussOnSphere( double sigma, double trunc_angle, double prec  )
  : m_cta(-1.0),
    m_circleint_k1(-1.0),
    m_circleint_k2(-1.0),
    m_norm(-1.0),
    m_expfact(-1.0),
    m_truncangle(-1.0),
    m_sigma(-1.0),
    m_numint_accuracy(-1.0),
    m_prec(-1),
    m_sta(-1.0)
{
  nc_assert(!isValid());
  set(sigma,trunc_angle,prec);
  nc_assert(isValid());
}

NC::GaussOnSphere::~GaussOnSphere()
{
#ifndef NDEBUG
  if (m_stats.genpointworst)
    produceStatReport("destructed");
#endif
}

void NC::GaussOnSphere::produceStatReport(const char * callpt)
{
  markused(callpt);
#ifndef NDEBUG
  uint64_t worst = (m_stats.genpointcalled?(uint64_t)m_stats.genpointworst:0);
  NCRYSTAL_MSG("GaussOnSphere(sigma="<<m_sigma<<", truncangle="
               <<m_truncangle/m_sigma<<"sigma, prec="<<m_prec<<") "
               <<callpt<<". Used "<<m_stats.genpointtries
               <<" tries to generate "<<m_stats.genpointcalled
               <<" pts on circles (acceptance rate: "
               <<(m_stats.genpointtries?m_stats.genpointcalled*100.0/m_stats.genpointtries:0.0)
               <<"%). Worst case used "<<worst<<" tries."
               << " Performed "<<m_stats.circleintnumber
               <<" numerical circle integrations using an average of "
               <<(m_stats.circleintnumber?double(m_stats.circleintevals)/m_stats.circleintnumber:0.0)
               << " function evaluations each time (worst case used "
               <<m_stats.circleintworst<<" evaluations).");
#else
  markused(callpt);
#endif
}

void NC::GaussOnSphere::set(double sigma, double trunc_angle, double prec ) {
  nc_assert_always(sigma>0);
  nc_assert_always(trunc_angle>0);
  nc_assert_always(trunc_angle<kPiHalf);
  if ( ! (valueInInterval(0.9999e-7,0.10000001,prec) || valueInInterval(1.0,10000.0,prec)) )
    NCRYSTAL_THROW(BadInput,"prec must either be in the range [1e-7,1e-1] or in the range [1,10000].");
  if (prec<=1)
    prec = ncclamp(prec,1e-7,0.1);
  else
    prec = ncclamp(prec,1.0,10000.0);

  if (m_truncangle==trunc_angle&&m_sigma==sigma&&m_prec == prec)
    return;

#ifndef NDEBUG
  if (m_stats.genpointworst) {
    produceStatReport("settings changed.");
    m_stats.genpointworst = m_stats.genpointtries = m_stats.genpointcalled = 0;
    m_stats.circleintworst = m_stats.circleintnumber = m_stats.circleintevals = 0;
  }
  if (ncgetenv_bool("DEBUG_GAUSSONSPHERE")) {
    m_stats.genpointworst = 1;//enable stat collection
    m_stats.genpointtries = m_stats.genpointcalled = 0;
    m_stats.circleintworst = m_stats.circleintnumber = m_stats.circleintevals = 0;
  }
#endif
  unsigned nlt = 0;
  m_truncangle=trunc_angle;
  m_sigma=sigma;
  m_prec=prec;
  sincos(trunc_angle,m_cta,m_sta);
  m_expfact = -0.5/(sigma*sigma);
  m_norm = calcNormFactor(sigma,trunc_angle);

  //analyse prec and trunc_angle parameters to establish some overall settings:
  nc_assert(prec>0);
  if (prec<1.0) {
    double circleint_approx_maxangle = gos_cosexpansionradius(prec*0.5);//0.5 is safety
    nc_assert_always(circleint_approx_maxangle<kPiHalf);
    if (m_truncangle>circleint_approx_maxangle||circleint_approx_maxangle<=1e-10) {
      //approximation formula never valid at requested precision:
      m_circleint_k1 = kInfinity;
      m_circleint_k2 = 0.0;
    } else {
      m_circleint_k1 = cos_mpi2pi2(circleint_approx_maxangle);
      m_circleint_k2 = m_cta - 1e-5;
      //NB: 1e-5 above is for numerical stability, so we are sure to always fall
      //back to numerical integration when gamma~=0
    }
    m_numint_accuracy = ncclamp(prec * 0.1,1e-7,1e-4);
    nlt = std::max<unsigned>(40,std::min<unsigned>(10000,unsigned(10.0/std::sqrt(prec)+0.5)));
  } else {
    nc_assert_always(prec>=1.0);
    //prec>1: always use approximation formula, use prec pts for lookup tables (but at least 20)
    nlt = std::max<unsigned>(20,unsigned(prec+0.5));
    m_numint_accuracy = 1e-4;//won't actually be used
    m_circleint_k1 = 0.0;
    m_circleint_k2 = kInfinity;
  }
  nc_assert_always(nlt>=20&&nlt<=10000);
  //For extreme (sub-arcsecond mosaicity) crystals, the numerical integration
  //suffers from instability issues to such a degree that we get convergence
  //issues if we request extreme accuracy. Since such crystals are more of an
  //edge case, we simply reduce the requested accuracy in these scenarios a bit
  //(it will still be highly accurate of course):
  if ( m_sigma< 10*kArcSec) {
    m_numint_accuracy = ncmax(m_numint_accuracy,1e-6);
    nlt = ncmin(nlt,1000);
    if ( m_sigma< kArcSec ) {
      nlt = ncmin(nlt,500);
      m_numint_accuracy = ncmax(m_numint_accuracy,1e-5);
      if ( m_sigma< 0.1*kArcSec ) {
        nlt = ncmin(nlt,200);
        m_numint_accuracy = ncmax(m_numint_accuracy,1e-4);
      }
    }
  }

  //Prepare lookup table for N*exp(-k*delta^2) by cos(delta) value, thus
  //avoiding calls to both acos and exp. The lookup table essentially
  //covers the interval from ltmin = cos(trunc_angle*1.01 to 1.0.

  double ltmin(0.), ltinvdelta(0.);
  SplinedLookupTable lt_evalcosx;
  {
    double nltm1 = nlt-1.0;
    ltmin = cos_mpi2pi2(ncmin(kPiHalf,trunc_angle));
    ltinvdelta = nltm1/(1.0-ltmin);

    SLTFct_EvalCosX sltfct(m_norm,m_expfact);
    //determine h, "epsilon" parameter for derivative estimation. It can not be
    //too small, or numerical impression results.
    double h = 0.01/ltinvdelta;
    double minh = ncmin(1e-4,(1.0-ltmin)*0.1);
    h=ncmax(h,minh);
    double fprime_b = estimateSingleSidedDerivative(&sltfct, 1.0,-h);
    double fprime_a = estimateSingleSidedDerivative(&sltfct, ltmin,h);
    lt_evalcosx.set(&sltfct,ltmin,1.0,fprime_a,fprime_b,nlt,
                    "gosdensityfromcosx",
                    "Gaussian density on sphere as a function of cos(x)");
  }

  SplinedLookupTable lt_sofcosd;
  {
    const double cta = cos_mpi2pi2(trunc_angle);
    SLTFct_SofCosD fct_sofcosd(m_norm,m_expfact,m_sigma,m_truncangle);
    //determine h, "epsilon" parameter for derivative estimation. It can not be
    //too small, or numerical impression results.
    double h = 0.01*(1.0-cta)/nlt;
    double minh = ncmin(1e-4,(1.0-cta)*0.1);
    h=ncmax(h,minh);
    double fprime_b = estimateSingleSidedDerivative(&fct_sofcosd, 1.0,-h);
    //Actual derivative at cta is infinite, so to avoid
    //oscillatory spline behaviour, it is better if we simply force the curve at
    //the first point to have a derivate pointing at the second point:
    double fprime_a = fct_sofcosd.eval(cta+(1.0-cta)/nlt)/((1.0-cta)/nlt);
    lt_sofcosd.set(&fct_sofcosd,cta,1.0,fprime_a,fprime_b,nlt,
                   "gosapproxcircleintegral",
                   "Approximate value of circle integral as a function of cos(alpha-gamma). Must be multiplied by sqrt(sin(alpha)/sin(gamma)).");
  }

  m_lt_sofcosd.swap(lt_sofcosd);
  m_lt_evalcosx.swap(lt_evalcosx);
}

double NC::GaussOnSphere::circleIntegralSlow( double cg, double sg, double ca, double sa ) const
{
  const double sasg = sa*sg; nc_assert(sasg>=0.0);
  const double cacg = ca*cg;
  const double cd = cacg+sasg;

  if (cd<=m_cta)
    return 0.0;

  if (sasg<1e-14) {
    //degenerate case, either sa~=0 and the result should vanish since the
    //circle has vanishing circumference, or sg~=0 and the density will be
    //roughly identical on all points of the circle and we can perform the
    //integral by evaluating the Gaussian at x=alpha and multiplying by the
    //circle's circumference.
    return k2Pi*sa*evalCosX(ca);
  }

  //full numerical integration required:
  nc_assert(sasg>0);
  double cos_tmax =  (m_cta-cacg)/sasg;
  double tmax = ( cos_tmax<=-1.0 ? kPi : std::acos(NC::ncmin(1.0,cos_tmax)) );
  if ( tmax<=1e-12 )
    return 0.0;

  //Increase target accuracy for edge cases (similar to what was done based on
  //m_sigma in ::set(..)):
  double intacc = m_numint_accuracy;
  if ( tmax< 10*kArcSec) {
    intacc = ncmax(intacc,1e-6);
    if ( tmax< kArcSec ) {
      intacc = ncmax(intacc,1e-5);
    if ( tmax< 0.1*kArcSec )
      intacc = ncmax(intacc,1e-4);
    }
  }

#ifndef NDEBUG
  const bool statcollect = (m_stats.genpointworst>0);
#else
  const bool statcollect = false;
#endif
  GOSCircleInt gosci(this,sasg,cacg,intacc,statcollect);
  const double res = 2.0*sa*gosci.integrate(0,tmax);//sa is radius of curve, so comes from Jacobian.
  if (!statcollect)
    return res;
#ifndef NDEBUG
  m_stats.circleintworst = std::max<uint64_t>(m_stats.circleintworst,gosci.nEvals());
  ++m_stats.circleintnumber;
  m_stats.circleintevals += gosci.nEvals();
#endif
  return res;

}

bool NC::GaussOnSphere::genPointOnCircle( RNG& rng, double cg, double sg, double ca, double sa,
                                          double& ct, double& st ) const
{
  nc_assert(isValid());
  nc_assert(ncabs(cg*cg+sg*sg-1.0)<1e-6);
  nc_assert(ncabs(ca*ca+sa*sa-1.0)<1e-6);
  const double sasg = sa*sg; nc_assert(sasg>=0.0);
  const double cacg = ca*cg;
  const double cd = cacg+sasg;

  if (cd<=m_cta)
    return false;

  if (sasg<1e-14) {
    //Degenerate case, either sa~=0 and the density integrated over the circle
    //should vanish since the circle has vanishing circumference (so return
    //false), or sg~=0 and the density will be roughly identical on all points
    //of the circle, so pick random t:
    if (sa<1e-7)
      return false;
    std::tie(ct,st) = randPointOnUnitCircle(rng);
    return true;
  }

  //Find valid range of t values:
  nc_assert(sasg>0);
  double cos_tmax =  (m_cta-cacg)/sasg;
  if (cos_tmax>=1.0)//vanishing length of circle inside truncation zone
    return false;
  double tmax = ( cos_tmax<=-1.0 ? kPi : std::acos(cos_tmax) );

  //The highest contribution is at t=0, at which cos(delta) = cd. Generate t via MC-rejection.
  double densitymax = evalCosXInRange(cd)*1.00000001;//1.00000001 is overlay safety
  const int maxtriesplus1(1001);//we should usually use *much* fewer tries than this (averaging around 3-6 depending on parameters).
  int triesleft = maxtriesplus1;
  while (--triesleft) {
    ct = cos_mpipi( rng.generate()*tmax );//generate t uniformly in allowed range
    double cd_at_t = sasg*ct+cacg;
    double density_at_t = evalCosXInRange(cd_at_t);
    if ( density_at_t > densitymax ) {
      static bool first = true;
      if (first) {
        first = false;
        NCRYSTAL_WARN("Problems sampling with rejection method during GaussOnSphere::genPointOnCircle "
                      "invocation. Overlay value was not larger than actual cross-section value at sampled point "
                      "(overshot by factor of "<<(densitymax?density_at_t/densitymax:kInfinity)<<"). Further warnings"
                      " of this type will not be emitted.");
      }
    }
    if ( density_at_t > densitymax * rng.generate())
      break;
  }
#ifndef NDEBUG
  if (m_stats.genpointworst) {
    ++m_stats.genpointcalled;
    uint64_t triesused = maxtriesplus1-triesleft;
    m_stats.genpointtries += triesused;;
    m_stats.genpointworst = std::max<uint64_t>(m_stats.genpointworst,triesused);
  }
#endif
  if (triesleft<=0) {
    static bool first = true;
    if (first) {
      first = false;
      NCRYSTAL_WARN("Problems sampling with rejection method during GaussOnSphere::genPointOnCircle "
                    "invocation. Did not accept sampled value after "<<maxtriesplus1-1<<" attempts. Further warnings"
                    " of this type will not be emitted.");
    }
    return false;
  }
  st = std::sqrt(1.0-ct*ct);
  st = (rng.coinflip()?st:-st);//pick t in [-pi,pi], not just in [0,pi]
  return true;
}

double NC::GaussOnSphere::estimateNTruncFromPrec( double prec, double minval, double maxval )
{
  nc_assert(minval>0&&maxval>minval&&prec>=0.0);
  //special cases: reference mode always selects ntrunc=maxval, approx mode always ntrunc=minval
  if (prec==0)
    return maxval;
  if (prec>=1.0)
    return minval;
  nc_assert_always(prec>0.0&&prec<1.0);
#if 0
  //Volume within a radius of N*sigma from the center of a 1D Gaussian in the
  //plane is 1-erf(N/sqrt(2)), which can not be inverted easily, so we use
  //resort to root-finding:
  class EstNTruncFct : public Fct1D {
    double m_prec;
  public:
    EstNTruncFct(double prec) : m_prec(prec) {};
    virtual ~EstNTruncFct(){}
    virtual double eval(double x) const
    {
      return std::erfc(x*kInvSqrt2) - m_prec;
    }
  };
  EstNTruncFct f(prec);
  if (f.eval(maxval)>=0)
    return maxval;
  if (f.eval(minval)<=0)
    return minval;
  double n = findRoot(&f,minval,maxval,1e-13);
#else
  //Volume within a radius of N*sigma from the center of a 2D Gaussian in the
  //plane is 1-exp(-N^2/2), which can be use to find the following exact
  //formula:
  double n = std::sqrt(-2.0*std::log(ncmax(1e-300,prec)));
#endif
  //Returned clamped value and with a safety factor of 1.1:
  return ncclamp(1.1*n,minval,maxval);
}

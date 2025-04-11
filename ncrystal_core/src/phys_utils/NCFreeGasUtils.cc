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

#include "NCrystal/internal/phys_utils/NCFreeGasUtils.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRandUtils.hh"
namespace NC=NCrystal;

#define NCRYSTAL_FREEGASUTILS_ENABLEEXTRADEBUGGING 0

NC::FreeGasXSProvider::FreeGasXSProvider( Temperature temp_kelvin,
                                          AtomMass target_mass_amu,
                                          SigmaFree sigma )
{
  //Evaluate cross-sections via sampling paper (10.1016/j.jcp.2018.11.043)
  //eq. 20, with explicit prefactor of sigma_free=sigma_bound*(A/(1+A))^2" (see
  //also discussion after eq. 20 in
  //https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4927245/). And discussion on
  //DGSW-399.
  //
  //Sampling paper eq. 20 does not include A=M/mn in "a", but here we scale a by
  //sqrt(A). Not sure where we got this from, but it was heavily validated using
  //numerical integration of sampling paper eq. 19.
  temp_kelvin.validate();
  target_mass_amu.validate();
  sigma.validate();

  double A = target_mass_amu.relativeToNeutronMass();//mass_nuclide / mass_neutron
  m_sigmaFree = sigma.get();
  nc_assert_always( temp_kelvin.get() > 0.0);
  nc_assert_always( target_mass_amu.get() > 0.0);
  nc_assert_always( m_sigmaFree > 0.0);

  nc_assert( temp_kelvin.kT() > 0.0 );
  m_ca = A / temp_kelvin.kT();//a^2 = A*E/kT = ca*E with ca = A/kT [units 1/eV]:
}

NC::FreeGasXSProvider::FreeGasXSProvider( Temperature temp_kelvin,
                                          AtomMass target_mass_amu,
                                          SigmaBound sb )
  : FreeGasXSProvider( temp_kelvin, target_mass_amu, sb.free(target_mass_amu) )
{
}

NC::FreeGasXSProvider::~FreeGasXSProvider() = default;

double NC::FreeGasXSProvider::evalXSShapeASq(double a_squared)
{
  //a^2 = A*E/kT.
  if (a_squared>36.0)
    return 1.0 + 0.5 / a_squared;
  const double a = std::sqrt(a_squared);
  if (a<0.1) {
    if (a==0.0)
      return kInfinity;
    constexpr double c1 = 2.0/3.0;
    constexpr double c2 = 1.0/15.0;
    constexpr double c3 = 1.0/105.0;
    constexpr double c4 = 1.0/756.0;
    constexpr double c5 = 1.0/5940.0;
    const double a2 = a_squared;
    return kInvSqrtPi * ( 2.0 / a +  a *( c1- a2*(c2-a2*(c3-a2*(c4-a2*c5)))));
  };
  //intermediate region, fall-back to full formula (slow):
  const double inva = 1.0 / a;
  return ( 1.0 + 0.5*inva*inva ) * std::erf(a) + kInvSqrtPi * std::exp(-a_squared)*inva;
}

namespace NCRYSTAL_NAMESPACE {

  //////////////////////////////
  // Private helper utilities //
  //////////////////////////////

  namespace ErfcBounds {

    //////////////////////////////////////////////////////////////////////////////
    // Lookup-table based helper function for quickly getting bounds of erfc(x) //
    //////////////////////////////////////////////////////////////////////////////

    //pick parameters to span [-2,9] (= a bit more than needed), and so we have
    //roughly 1000 bins and a nice memorable binwidth of 0.01:
    constexpr std::size_t lookupTableLength = 1103;//= 3+(9.0-(-2.0))/0.01
    constexpr double lookupTableLowerEdge = -2.0;
    constexpr double lookupTableUpperEdge = 9.0;

    //Binwidth calcs us -3 instead of -1, since we also have two fictional bins
    //on each side going to +-infinity:
    constexpr double lookupTableBinWidth = (lookupTableUpperEdge-lookupTableLowerEdge)/(lookupTableLength-3);
    constexpr double lookupTableInvBinWidth = (lookupTableLength-3)/(lookupTableUpperEdge-lookupTableLowerEdge);

    static VectD initCache() {
      //Fill erfc(x) points for x at all bin edges, plus one at -infinity and one at +infinity:
      VectD v;
      v.reserve(lookupTableLength);
      constexpr int nbinedges = lookupTableLength-2;
      //sanity check: nbins*binwidth takes us from lower to upper edge:
      nc_assert(floateq(lookupTableLowerEdge+(nbinedges-1)*lookupTableBinWidth,lookupTableUpperEdge));
      v.push_back(2.0);//erfc(-infinity)
      for (auto e: linspace(lookupTableLowerEdge,lookupTableUpperEdge,nbinedges))
        v.push_back(std::erfc(e));
      v.push_back(0.0);//erfc(+infinity)
      nc_assert( v.size() == lookupTableLength );
      return v;
    }
    static VectD cache = initCache();

    PairDD erfcQuickBounds(const double x) {
      //returns {LB,UB} so LB <= erfc(x) <= UB
      constexpr int nbins_inclinf = lookupTableLength-1;//Number of bins, including those width edges at +-infinity
      constexpr double verylow = lookupTableLowerEdge-0.5*lookupTableBinWidth;//half a bin below lookupTableLowerEdge
      constexpr double veryhigh = lookupTableUpperEdge+0.5*lookupTableBinWidth;//half a bin above lookupTableLowerEdge
      const double x_safe = ncclamp(x,verylow,veryhigh);//clamp, to avoid numerical or overflow issues for huge arguments
      const double x_relpos = (x_safe-lookupTableLowerEdge);
      //Find index of bin (because clamping above, x_relpos/binwidth will be in
      //range [-0.5,lookupTableLength-2.5]. Once we add 1.0 and cast to int, we
      //will get a bin index in the range [0,lookupTableLength-2].
      const int binidx = std::max<int>(0,std::min<int>(nbins_inclinf,static_cast<int>(1.0+x_relpos * lookupTableInvBinWidth)));
      //Evaluating the cache at binidx will give us the lowerbound, and evaluating at binidx+1 will give us the upper bound.
      PairDD bounds;
      bounds.first = cache[binidx+1]*0.99999999;//lowerbound
      bounds.second = cache[binidx]*1.00000001;//upperbound
#if NCRYSTAL_FREEGASUTILS_ENABLEEXTRADEBUGGING
      double test_erfcx = std::erfc(x);
      nc_assert_always(test_erfcx>=bounds.first);
      nc_assert_always(test_erfcx<=bounds.second);
#endif
      return bounds;
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  // Helper class for evaluating f(beta|A,E). Precisely or with quicker bounds //
  ///////////////////////////////////////////////////////////////////////////////

  class FGEvalBetaDistHelper final {
    double m_beta, m_normfact, m_expmbeta;
    double m_k11, m_k12, m_k21, m_k22;
  public:
    FGEvalBetaDistHelper(double c, double invA, double sqrtAc, double beta, double normfact)
      : m_beta(beta),
        m_normfact(normfact),// = 1/(2*erf(sqrt(c/A))) to make f(beta=0)=1
        m_expmbeta(-1.0)
    {
      //distribution function of beta, normalised to unity at beta=0
      nc_assert(beta>=-c);
      nc_assert(c>0.0);
      const double eps = beta/c;
      const double sqrt1pluseps = std::sqrt(1+eps);
      const double S = (beta<0.0?-1.0:1.0);
      //epsprime = 1 + min(0,eps)
      const double sqrtepsprime = ( eps >= 0.0 ? 1.0 : sqrt1pluseps );
      const double sqrtgammaplus = std::sqrt(2.0+eps+2.0*sqrt1pluseps);
      const double SP = 0.5*(S+invA);
      const double SM = 0.5*(S-invA);
      const double invA_sqrtepsprime = invA*sqrtepsprime;
      const double mS_sqrtepsprime = -S*sqrtepsprime;
      const double SPsgp = sqrtgammaplus*SP;
      const double SMsgp = sqrtgammaplus*SM;
      const double j11 = -invA_sqrtepsprime + SPsgp;
      const double j12 = mS_sqrtepsprime + SPsgp;
      const double j21 = mS_sqrtepsprime + SMsgp;
      const double j22 =  invA_sqrtepsprime + SMsgp;
      m_k11 = sqrtAc*j11;
      m_k12 = sqrtAc*j12;
      m_k21 = sqrtAc*j21;
      m_k22 = sqrtAc*j22;
    }
    void evalExpMBeta()
    {
      //beta<-700 check is to ignore exp(largeval) overflow:
      if (m_expmbeta<0)
        m_expmbeta = m_beta<-700.0 ? 0.0 : std::exp(-m_beta);
    }

    double evalExact()
    {
      double t1 = erfcdiff(m_k11,m_k12);
      evalExpMBeta();
      if (!m_expmbeta)
        return m_normfact*t1;
      double t2 = erfcdiff(m_k21,m_k22);
      return m_normfact*(t1+t2*m_expmbeta);
    }

    PairDD evalQuickBounds()
    {
      //Unfortunately, the exp(-beta)*(erfc(k11)-erfc(k12)) term is unstable
      //for beta<0, so for now we abandon the idea of using the bounds technique for those:
      //Comment out, since testBetaDistEval might call with any beta: nc_assert(m_beta>0.0);

      //value pairs are {.first,.second} = {lowerbound,upperbound}
      auto erfck11 =  ErfcBounds::erfcQuickBounds(m_k11);
      auto erfck12 =  ErfcBounds::erfcQuickBounds(m_k12);
      PairDD t1( erfck11.first  - erfck12.second,//LB
                                   erfck11.second - erfck12.first );//UB


      auto erfck21 =  ErfcBounds::erfcQuickBounds(m_k21);
      auto erfck22 =  ErfcBounds::erfcQuickBounds(m_k22);

      PairDD t2( erfck21.first  - erfck22.second,//LB
                                   erfck21.second - erfck22.first );//UB

      if (t2.second>0.0) {
        evalExpMBeta();
        return std::make_pair<double,double>(m_normfact*(t1.first+t2.first*m_expmbeta),
                                             m_normfact*(t1.second+t2.second*m_expmbeta));

      }
      return std::make_pair<double,double>(m_normfact*t1.first,m_normfact*t1.second);
    }
  };

  ///////////////////////////////////////////////////////////////////
  // Utility for sampling f(x)=exp(-1/x-c*x)/sqrt(x) over [xm,xp]  //
  // (could go to NCRandUtils if useful outside free gas sampling) //
  ///////////////////////////////////////////////////////////////////

  double randExpMInvXMCXDivSqrtX( RNG& rng, double c, double xm, double xp )
  {
    //sample f(x)=exp(-1/x-c*x)/sqrt(x) over [xm,xp]
    //
    //This is rather challenging since the transformation method does not yield a
    //closed form expression (and the integral involves expensive and numerically
    //problematic error functions). Furthermore, the function drops off on both
    //sides of a central peak, while there are three free parameters, c, xm and
    //xp.

    nc_assert(xm>=0.0);
    nc_assert(xp>=xm);
    nc_assert(c>0.0);//ok, since sampleAlpha(..) takes care of c~=0.
    if (xp==xm)
      return xm;

    const double sqrtc=std::sqrt(c);
    nc_assert(sqrtc>0.0);
    const double invsqrtc = 1/sqrtc;

    //Peak position (full formula is numerically unstable at small c, so evaluate
    //with taylor expansion there):
    const double xpeak = ( c > 1e-5
                           ? ( c > 1e200 ? invsqrtc : (std::sqrt(16.0*c+1.0)-1.0)/(4.0*c) )
                           : ( 2.0-c*(8.0-c*(64.0-c*(640.0-c*7168.0))) ) );
    nc_assert(xpeak>=0&&xpeak<=2);
    if (xpeak==0.0) {
      //c is infinite or so large that xpeak is calculated as 0. Fall-back to
      //smallest non-zero allowed value.
      return xm > 0.0 ? xm : ncmin(std::numeric_limits<double>::min(),xp);
    }

    //The x-value with highest f(x) on [xm,xp]:
    const double xmax = (xm>xpeak?xm:ncmin(xp,xpeak));
    if (!(xmax>0.0))
      return xm;//could happen if xpeak and xm both evaluate to 0.0

    nc_assert(xmax>0.0);

    //NB: -1/x-c*x = -(c*x^2+1)/x. So if c*xm^2 >> 1, the exp(-1/x) effect is
    //negligible. On the other hand, it dominates when c*xp^2<<1. We pick two
    //thresholds for that effect:
    double xlarge = ncmax(5.0/sqrt(c),2*xpeak);
    double xsmall = ncmin(0.2/sqrt(c),0.5*xpeak);

    //Now truncate [xm,xp] to get rid of edge-areas where the sampled function is
    //completely negligible. This increases efficiency and numerical stability.
    if (xp > xlarge) {
      //Above xlarge exp(-c*x) dominate. So moving K/c up above xlarge makes
      //further contributions negligible when exp(-K) is neglible (we choose
      //K=15).
      xp = ncmin( xp, ncmax(xm,xlarge) + 15.0/c );
    }
    if (xm < xsmall) {
      //Below xsmall, -1/x dominates. We move down so the factor -1/x becomes K
      //smaller (replace xsmall with xp if xp<xsmall):
      //
      //    -1/xm = -1/xsmall - K <=> xm = xsmall / ( 1 + K*xsmall )
      //
      //We choose K=30 to be safe, since the ignored 1/sqrt(x) is more important
      //for small x.
      double xsm = ncmin(xp,xsmall);
      xm = ncmax( xm, xsm / ( 1 + 30.0*xsm ) );
    }

    //guard against zero-division:
    //fpmin is smallest positive number with full precision (non-denormal)
    constexpr double fpmin = std::numeric_limits<double>::min();
    if ( (xm = ncmax(fpmin,ncmax(fpmin/xp,xm))) >= xp )
      return xp;

    //cutoff regions where f(x) is less than this (for efficiency).
    constexpr double fcutoff_limit = 1e-9;

    //Setup f(x) normalised so f(xmax)=1:
    auto f_eval = [xmax,c](double x) {
                    const double exparg = (x-xmax)/(x*xmax)-c*(x-xmax);
                    if (exparg>=706.0)
                      return 1.0;//hopefully this only happens near function maximum
                    const double fval = exparg < -745.1 ? 0.0 : std::exp(exparg)*std::sqrt(xmax/x);
                    nc_assert( fval < 1.00001 );
                    return fval;
                  };

    //speed up some edge-cases by very aggressively checking if xm can be moved
    //towards xp. If not, the cost is at most one function evaluation.
    nc_assert(floateq(f_eval(xmax),1.0));
    if ( xp < xpeak ) {
      while (true) {
        double xm_new = xp - 0.01*(xp-xm);
        double fval = f_eval(xm_new);
        if ( fval>=fcutoff_limit )
          break;
        xm = xm_new;
      }
    }

    //Now, try to determine which of the two chosen sampling methods is optimal for the
    //values of c, xm, and xp:
    //
    // Method A ("flat")     : Sample with flat overlay function.
    // Method B ("largex")   : For large x the exp(-1/x) behaviour is unimportant, and the
    //                         optimal sampling is performed with an overlay function of
    //                         f_overlay(x) = exp(-c*x)/sqrt(x).
    // Method C ("combined") : Neither A nor B are effective, so we combine them.
    //
    // In both cases, we improve efficiency dramatically by dynamically shrinking
    // [xm,xp] whenever encountering sampled points of the overlay function with
    // abysmal acceptance rates (the current code is expected to be called with
    // different parameters each time, so such a narrowing has to be on-the-fly
    // and cheap rather than carefully and expensively precalculated). In addition
    // we also avoid many function evaluations by reusing knowledge obtained in
    // earlier samplings about when the function drops below certain threshold
    // values.

    double probability_flat(-1.0);//=1 => method A, =0 => method B, in between => method C.
    double xswitch(-1.0);//upper point of "flat" sampling, lower point of "largex" sampling.
    double area_right(-1.0);//overlay function area to the right of
    //xswitch. Calculate exactly once and reuse, as it is
    //expensive.

    //Now analyse situation regarding overlay functions:

    if ( xm >= xlarge ) {
      //we are far above the peak, out on the tail.
      probability_flat = 0.0;//
      xswitch = xm;
    } else if ( c > 25 || xp <= xlarge ) {
      //peak drops off sharply on the right side (c large) or we dont have to
      //sample the tails (xp<=xlarge), so a simple global flat overlay should work
      //fine.
      probability_flat = 1.0;
      xswitch = xp;
    } else {
      //Combined mode!
      nc_assert(xm<xlarge&&xp>xlarge);
      xswitch = xlarge;
      nc_assert(xmax<xswitch);
      //We must carefully calculate overlay function areas to each side of
      //xswitch, so we can randomly pick the two sides with the right probability
      //later.

      //We normalise f(x) to have a value of 1 at xmax. The flat overlay to the
      //left of xswitch then has a height of 1 and an area of:
      const double area_left = (xswitch-xm);
      //The non-flat overlay to the right of xswitch
      // is: f1(x)=N*exp(-c*x)/sqrt(x), with N chosen
      // so that f1(xp)==f(xp). This can be shown to
      // result in an area to the right of xswitch
      // of:
      //    k1*(erfc(sqrt(C*xswitch))-erfc(sqrt(C*xp)))
      // with:
      //    k1=exp(c*xmax+1/xmax-1/xp)*sqrt(pi*xmax/c)
      const double B=c*xmax+1/xmax-1/xp;
      area_right=(erfc_rescaled(sqrtc*std::sqrt(xswitch), B)
                  -erfc_rescaled(sqrtc*std::sqrt(xp), B))*sqrt(kPi*(xmax/c));
      nc_assert(area_left+area_right>0.0);
      probability_flat = area_left/(area_left+area_right);
    }

    //Generate!
    auto eval_probability = [fcutoff_limit](double pleft)
                            {
                              (void)fcutoff_limit;//avoid annoying misleading clang warning
                              nc_assert(pleft>=0.0&&pleft<=1.0);
                              bool always_left = (pleft>1.0-fcutoff_limit);
                              bool always_right = (pleft<fcutoff_limit);
                              bool single_side = (always_left||always_right);
                              return std::make_pair<double>(single_side,always_left);
                            };

    bool single_side, always_left;
    std::tie(single_side,always_left) = eval_probability(probability_flat);

    //Edge-case - for consistency we must revert combined mode to flat-only mode
    //if f_eval(xswitch) is already below cut-off (1.1 factor for numerical
    //safety):
    if ( ! single_side && f_eval(xswitch) < fcutoff_limit*1.1 ) {
      probability_flat = 1.0;
      area_right = 0.0;
      xp = xswitch;
      always_left = true;
      single_side = false;
    }

    //Start sampling loop:

    while (true) {
      bool do_flat(single_side?always_left:(rng.generate()<probability_flat));
      if (do_flat) {
        //==> Sampling with flat overlay in [xm,xswitch]
        nc_assert(single_side==always_left);
        double dx = xswitch-xm;
        nc_assert(valueInInterval(xm,xswitch,xmax));
        double xthr_up(xswitch), xthr_low(xm);
        constexpr double fthreshold(0.05);
        double xgen = xm + rng.generate()*dx;
        double Raccept = rng.generate();
        //Idea, since it is cheap to modify area_left. Divide into three flat regions... and adjust area_left etc.
        if ( !valueInInterval(xthr_low,xthr_up,xgen) && Raccept>fthreshold )
          continue;//reject cheaply using previously obtained knowledge

        //Expensive function evaluation:
        const double fval = f_eval(xgen);

        if (fval<fthreshold) {

          //low fval, record this knowledge for cheap rejection of future values:
          if (xgen < xmax)
            xthr_low = xgen;//postponed_fixme: Never read, modifying a local loop varibale which was not intended! (same for xthr_up)
          else
            xthr_up = xgen;//F i x m e: Never read, modifying a local loop varibale which was not intended!

          if (fval < fcutoff_limit) {
            //generated value is in negligible region. Make sure we don't generate here again.
            if ( xgen < xmax ) {
              xm = xgen;
            } else {
              nc_assert( always_left );
              xswitch = xgen;
            }
            dx = xswitch-xm;
            if (!single_side) {
              //Recalculate left/right probabilities:
              const double area_left = dx;
              nc_assert( area_right>=0.0 && area_left+area_right>0.0 );
              probability_flat = area_left/(area_left+area_right);
              std::tie(single_side,always_left) = eval_probability(probability_flat);
            }
            continue;//Reject value.
          }
        }
        if ( Raccept <= fval )
          return xgen;//Accept generated value
      } else {
        //==> Tail sampling in [xswitch,xp]
        //
        //Sample with the rejection method, using exp(-c*x)/sqrt(x) as an overlay
        //function, and reject if R>exp(1/xp-1/x) (which is 1 at x=xp, in
        //accordance with the normalisation assumed when calculating
        //probability_flat). Unlike for the flat case, we do not look for or apply
        //fcutoff_limit (it is not needed, and it is expensive to recalculate
        //area_right each time).
        nc_assert(!always_left);
        double xgen = randExpDivSqrt( rng, c, xswitch, xp );

        nc_assert(xgen>=xswitch&&xgen<=xp);
        if ( rng.generate() < std::exp((xgen-xp)/(xgen*xp)))
          return xgen;
      }
    }
  }

}

NC::FreeGasSampler::FreeGasSampler(NeutronEnergy ekin, Temperature temp_kelvin, AtomMass target_mass_amu)
  : m_c(ncmin(1e14,ncmax(1e-10,ekin.get()/(temp_kelvin.kT())))),
    //nb: we constrain m_c for numerical safety (1.16e13 is 1GeV neutron on 1Kelvin material)
    m_kT(temp_kelvin.kT()),
    m_sqrtAc(std::sqrt(target_mass_amu.get()*m_c/const_neutron_atomic_mass)),
    m_invA(1.0/target_mass_amu.relativeToNeutronMass()),
    m_Adiv4(0.25*target_mass_amu.relativeToNeutronMass()),
    m_normfact(0.5/std::erf(std::sqrt(m_c*m_invA))),
    m_c_real(ekin.get()/(temp_kelvin.kT()))//unconstrainted version of m_c
{
#ifndef NDEBUG
  nc_assert(m_c>0);
  nc_assert(m_c_real>=0);
  nc_assert(ekin.get()>=0);
  nc_assert(temp_kelvin.get()>0);
  nc_assert(m_kT>0);
  nc_assert(target_mass_amu.get()>0);
  nc_assert(m_invA>0);
  nc_assert(m_sqrtAc>0.0);
  ekin.validate();
  temp_kelvin.validate();
  target_mass_amu.validate();
#endif
}

NC::FreeGasSampler::~FreeGasSampler() = default;

void NC::FreeGasSampler::testBetaDistEval (double beta, double & f_exact, double & f_lb, double& f_ub )
{
  if ( beta <= -m_c ) {
    f_lb = f_ub = f_exact = 0.0;
    return;
  }
  FGEvalBetaDistHelper eval_helper(m_c, m_invA, m_sqrtAc, beta, m_normfact);
  std::tie(f_lb, f_ub) = eval_helper.evalQuickBounds();
  f_exact = eval_helper.evalExact();
}

double NC::FreeGasSampler::sampleBeta( RNG& rng ) const
{
  if (m_c_real>1e4) {
    //At extremely high energy the neutron will to a very good approximation
    //experience an energy-loss which is uniformly distributed from deltaE=0 to
    //deltaE=-E' with E' being the maximally kinematically allowed limit. The
    //exact threshold of what constitutes "very high' energy depends on the
    //value of A. We use an ad-hoc determination of the threshold, found from
    //inspection of beta-distributions:
    const double A = 1.0 / m_invA;
    const double A2 = A*A;
    const double c_highe_threshold = 1e4*ncmin(1000.0*A,A2*A2*A2);
    if ( m_c_real > c_highe_threshold ) {
      double am1_div_ap1 = (1.0-m_invA)/(1.0+m_invA);// (A-1)/(A+1);
      double elossmax = m_c_real * (1.0-am1_div_ap1*am1_div_ap1);
      return -elossmax*rng.generate();
    }
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  //The idea is to sample f(beta) with the rejection method, using an overlay
  //function, f_overlay, which is f_overlay=0 outside [a,b], f_overlay=1 on
  //[a,0] and f_overlay=exp(-beta) on [0,b]. For a<0<b, and with a=-m_c and
  //b=infinity. However, simply pursuing this strategy blindly results in
  //abysmal computation speed. Therefore, several optimisation strategies are
  //employed, outlined in the following:
  //
  //--------->
  //
  //OPTIMISATION #1) Using something similar to, but faster than,
  //f_overlay=exp(-beta) for beta>0. Otherwise, a lot of computational time
  //would be wasted on sampling at beta>0, both in the std::log call used for
  //sampling and the subsequent std::exp(-beta) call to find the value of the
  //overlay function (NB: This would not be so bad if we had a very fast
  //ziggurat-like exponential sampling method implemented!). So for the interval
  //[0,2] (or [0,b] if b<2), we use a polynomial approximation to exp(-x), T(x)
  //which is cheap to evaluate and can itself be sampled efficiently (acceptance
  //rate from 44%@b~=2 to 100%@b~=0) with the rejection method and a flat
  //overlay. For T(x) we use the 6th order Taylor expansion of exp(-x):
  //
  //  T(x) = 1/720*x^6 - 1/120*x^5 + 1/24*x^4 - 1/6*x^3 + 1/2*x^2 - x + 1
  //
  //Whose integral from 0 to b is:
  //
  //  1/5040*b^7 - 1/720*b^6 + 1/120*b^5 - 1/24*b^4 + 1/6*b^3 - 1/2*b^2 + b
  //
  //And we also take advantage of the fact that T(x) is itself bounded on [0,2]
  //by the cheap function 1-(19/45)*x.
  //
  //--------->
  //
  //OPTIMISATION #2) Sampling f_overlay on [a,b]=[-m_c,infinity] would give
  //abysmal acceptance rates in most scenarios (the notable exception would be
  //for A=1). We take advantage of the fact that f(beta) drops off monotonically
  //away from beta=0, and employ an fcutoff_limit=1e-6. In principle we can then
  //narrow [a,b] to [a',b'] where f(a')=f(b')=fcutoff_limit. The only problem
  //here is that there is no analytical formula for a' and b', and finding them
  //with numerical methods (like root finding) would be very costly itself, and
  //thus defeat the purpose. Instead we perform a very aggressive and crude
  //pre-sampling search to narrow [a,b], exponentially narrowing [a,b] at each
  //step. That handles the most extreme cases (i.e. A>>1, c>>1) only, so once we
  //proceed to actual sampling, we will narrow [a,b] dynamically whenever we are
  //so unlucky to sample a beta value where we subsequently find
  //f(beta)<fcutoff_limit. That ensures an exponentially fast convergence, and
  //is mathematically sound - technically it simply means that once in a while,
  //after rejecting a sampled value, we abandon further sampling attempts with
  //the current overlay function, and continue after redefining f_overlay to
  //another valid form.
  //
  //--------->
  //
  //OPTIMISATION #3) Related to the previous optimisation (#2), but instead of
  //cropping [a,b] and avoiding sampling beta points where f(beta) is
  //negligible, it merely avoids evaluation of f(beta) at many of the points
  //within [a,b]. Specifically, based on the f(beta) evaluations that anyway
  //have to be performed during the course of rejected samplings, we keep track
  //of [afthreshold,bfthreshold], the interval outside of which we are certain
  //f(beta)<fthreshold, with fthreshold=0.1. That allows us to reject future
  //values outside this interval 90% of the time without an expensive function
  //evaluation. In a way, this method can be thought of as a dynamic ziggurat
  //method with 1 step on the ziggurath.
  //
  //--------->
  //
  //OPTIMISATION #4) The evaluation of f(beta) is costly, mostly since it
  //requires 4 calls to std::erfc. Before evaluating f(beta), we employ a lookup
  //table of erfc values to cheaply construct upper and lower limits of
  //f(beta). This allows us to accept/reject many points without actually
  //evaluating f(beta) exactly. Unfortunately, this method turns out to only
  //work for beta>0, as the bounds for beta<0 are often too imprecise to be of
  //any use.
  //
  //--------->
  //
  //A final note is that the evaluation of f(beta) is complicated by the fact
  //that it contains two erfc(x1)-erfc(x2) terms, which must be evaluated very
  //carefully in order to always produce correct results.
  //
  ////////////////////////////////////////////////////////////////////////////////

  //Constants concerning T(x):
  constexpr double Tlim = 2.0;//NB: Don't change this value without changing the next two constants!
  constexpr double Tlim_k1 = 0.135335283236612691893999494972484403407;//std::exp(-Tlim) for Tlim=2.
  constexpr double Tlim_k2 = 274./315.;//integral of T(x) from 0 to Tlim=2.

  //The fcutoff_limit:
  const double fcutoff_limit = 1e-6;

  //Start by sampling in [-E/kT,bmax] where bmax is chosen so exp(-bmax) is
  //negligible. And use -m_c_real instead of -m_c as lower sampling limit for
  //extremely low-energy neutrons to not violate kinematic constraints:
  double aa(ncmax(-m_c_real,-m_c)), bb(13.815510557964274);//bb=-log(fcutoff_limit)
  nc_assert(aa>=-m_c_real);
  nc_assert(floateq(std::exp(-bb),fcutoff_limit));

  //At high energies and A values, f is more narrowly peaked around beta=0, so
  //we gain efficiency by narrowing the initial interval more strongly than the
  //dynamic narrowing below would do:

  if ( m_invA <= 1.0/10.0 ) {//only for A>=10
    if (m_c>10.1) {
      //lower limit - don't bother for low-energy particles.
      while (true) {
        double aa_new = aa*0.2;//stepsize tuned
        if (aa_new>-1e-99)
          break;
        double fval = FGEvalBetaDistHelper(m_c, m_invA, m_sqrtAc, aa_new, m_normfact).evalExact();
        if ( fval > fcutoff_limit )
          break;
        aa = aa_new;
      }
    }
    //Same for bb (but can use quick upper bound for beta>0):
    {
      while (true) {
        double bb_new = bb*0.25;//stepsize tuned
        if ( bb_new<1e-99 )
          break;
        double fval_upperbound = FGEvalBetaDistHelper(m_c, m_invA, m_sqrtAc, bb_new, m_normfact).evalQuickBounds().second;
        if ( fval_upperbound > fcutoff_limit )
          break;
        bb = bb_new;
      }
    }
  }
  if (!(bb>aa))
    return aa;

  //Variables defining overlay (initialised with setAB):
  double a, b, prob_downscat, prob_notclosetail;
  RandExpIntervalSampler expsampler;

  auto setAB = [&a,&b,&prob_downscat,&prob_notclosetail,&expsampler,Tlim,Tlim_k1,Tlim_k2](double aaa,double bbb) {
                 (void)Tlim;//avoid annoying misleading clang warning
                 (void)Tlim_k1;//avoid annoying misleading clang warning
                 (void)Tlim_k2;//avoid annoying misleading clang warning
                 a = aaa;
                 b = bbb;
                 double area_downscat(-aaa),area_fartail,area_closetail;
                 if (bbb<=Tlim) {
                   area_fartail = 0.0;
                   constexpr double c2 = -1./2.;
                   constexpr double c3 =  1./6.;
                   constexpr double c4 = -1./24.;
                   constexpr double c5 =  1./120.;
                   constexpr double c6 = -1./720.;
                   constexpr double c7 =  1./5040.;
                   area_closetail =  b * (1.0+b*(c2+b*(c3+b*(c4+b*(c5+b*(c6+b*c7))))));
                   nc_assert(area_closetail>=0.0&&area_closetail<=Tlim_k2);
                 } else {
                   area_fartail = Tlim_k1 - std::exp(-bbb);
                   nc_assert(area_fartail>=0.0);
                   area_closetail = Tlim_k2;
                 }
                 nc_assert(area_fartail<=area_closetail);
                 nc_assert(area_downscat>=0.0&&area_closetail>=0.0&&area_fartail>=0.0);
                 double invareatot = 1.0 / ( area_downscat + area_closetail + area_fartail );
                 prob_downscat = area_downscat * invareatot;
                 prob_notclosetail = ( area_fartail + area_downscat ) * invareatot;
                 expsampler.invalidate();
                 nc_assert( prob_downscat>=0.0&&prob_downscat<=1.0 );
                 nc_assert( prob_notclosetail>=0.0&&prob_notclosetail<=1.0 );
                 nc_assert(a<=0&&b>=0);
               };
  //initialise:
  setAB(aa,bb);
  constexpr double fthreshold = 0.1;
  double afthreshold(a), bfthreshold(b);

  //Sampling loop:

  while (true) {
    double beta, foverlay;

    //////////////////////////////
    // Sample overlay function: //
    //////////////////////////////

    double R_selectregion = rng.generate();
    if (R_selectregion<prob_downscat) {
      //downscat, f_overlay(beta)=1
      beta = rng.generate()*a;
      nc_assert(beta>=a&&beta<=0.0);
      foverlay = 1.0;
    } else {
      if ( R_selectregion < prob_notclosetail ) {
        //far tail, f_overlay(beta)=exp(-beta)
        if (!expsampler.isValid())
          expsampler.set(Tlim,b,1.0);
        beta = expsampler.sample(rng);
        nc_assert(beta>=Tlim&&beta<=b*1.000001);
        foverlay = std::exp(-beta);
      } else {
        //close tail, f_overlay(beta)=T(beta).
        //Generate with rejection method, using flat overlay.
        double bmax=ncmin(b,Tlim);
        while (true) {
          beta = rng.generate()*bmax;
          double Raccept0 = rng.generate();
          constexpr double kcheap = 19./45.;
          if (Raccept0 > 1.0 - kcheap*beta)
            continue;//cheaply reject 75% to 100% of values which would be rejected below.
          //Must evaluate T(beta) (which becomes foverlay):
          constexpr double c1 = -1.;
          constexpr double c2 = 1./2.;
          constexpr double c3 = -1./6.;
          constexpr double c4 = 1./24.;
          constexpr double c5 = -1./120.;
          constexpr double c6 = 1./720.;
          foverlay =  1.0+beta*(c1+beta*(c2+beta*(c3+beta*(c4+beta*(c5+beta*c6)))));
          nc_assert(foverlay>=std::exp(-beta)*(1.0-1e-15));
          nc_assert(foverlay<=1.0);
          if (Raccept0<foverlay)
            break;
        }
      }

    }

    nc_assert(beta>=a&&beta<=b*1.000001);
    nc_assert(beta>-m_c_real);

    ///////////////////////////////////////////////////////////////////////////////
    // Deal with the sampled value - accept, reject, update cutoffs & thresholds //
    ///////////////////////////////////////////////////////////////////////////////

    //Generate random value uniformly in [0,foverlay] for acceptance test:
    const double faccept = rng.generate()*foverlay;

    if ( faccept > fthreshold && !valueInInterval(afthreshold,bfthreshold,beta) )
      continue;//we can reject here without function evaluation

    //In the following we attempt to carry out cheap estimation of f(beta)
    //bounds, in order to perform our rejection/acceptance without a more
    //expensive evaluation. Unfortunately, the bounds calculations are more
    //expensive and less useful at beta<0, since they rely on careful
    //cancellation of erfc() results. For that reason, using the bounds at
    //beta<0 makes the fcutoff method unreliable, and in the end the result is
    //not an improvement.

    FGEvalBetaDistHelper eval_helper(m_c, m_invA, m_sqrtAc, beta, m_normfact);
    bool need_exact(true);
    double fval;
    if (beta>0) {
      //At beta>0, we can get cheap and normally useful bounds on f(beta):

      double fval_lowerbound, fval_upperbound;
      std::tie(fval_lowerbound, fval_upperbound) = eval_helper.evalQuickBounds();

#if NCRYSTAL_FREEGASUTILS_ENABLEEXTRADEBUGGING
      double testfval = eval_helper.evalExact();
      nc_assert_always(testfval>=fval_lowerbound);
      nc_assert_always(testfval<=fval_upperbound);
#endif

      if ( faccept <= fval_lowerbound )
        return beta;//Accepted (cheaply)!

      //For further tests and updating of cutoff/thresholds, we need to do a full
      //expensive evaluation of f(beta) only when faccept falls in the (hopefully
      //narrow) band [fval_lowerbound,fval_upperbound]. Otherwise we can proceed
      //using fval_upperbound:
      fval = fval_upperbound;
      if (faccept > fval_upperbound)
        need_exact = false;
    }

    if ( need_exact ) {
      fval = eval_helper.evalExact();//Expensive evaluation!
      nc_assert( -1e-20 <= fval&&fval<=1.0000001 );
      nc_assert( fval<=foverlay*1.0000001 );
      if ( faccept < fval )
        return beta;//Accepted (expensive)!
    }

    //Rejected (cheaply or expensively). Check fval to see if we should update
    //any parameters for more efficient sampling in subsequent tries.

    if (fval<fcutoff_limit) {
      //Negligible contribution here. Narrow [a,b] to avoid wasting time here
      //again (this is technically/mathematically OK, it simply corresponds to
      //switching to a new overlay function for future tries)!
      if (beta<0)
        setAB(beta,b);
      else
        setAB(a,beta);
      continue;
    }

    if (fval<fthreshold) {
      //cache new-found knowledge of where function is < fthreshold:
      if ( beta<0 )
        afthreshold = ncmax(afthreshold,beta);
      else
        bfthreshold = ncmin(bfthreshold,beta);
    }

  }
}

double NC::FreeGasSampler::sampleAlpha( double beta, RNG& rng ) const
{
  // We must sample alpha in [alpha-,alpha+] from the free-gas S(alpha,beta),
  // which is proportional to (A=M/m_neutron):
  //
  // exp( -beta^2/(4alpha/A) - alpha/4A ) / sqrt(alpha) = exp(-t/alpha - s*alpha)/sqrt(alpha)
  //
  // with t=beta^2*A/4 and s=1/4A.
  //
  // Now, for beta!=0 we can put x=alpha/t and sample x instead. This means we
  //  should sample x in [alpha-/t,alpha+/t] from the distribution:
  //
  //  exp(-1/x - (s*t)*x)/sqrt(x) =  exp(-1/x - c*x)/sqrt(x)
  //
  // with c = s*t = beta^2/16  [NOTE: this is a different c than m_c which is E/kT]
  //
  // After sampling x, we of course get alpha as t*x
  //
  // For the case where beta ~= 0, ...

  if ( m_c_real < m_c || muIsotropicAtBeta(beta,m_c) ) {
    nc_assert( beta >= -m_c_real*1.001 );
    //Either we are close to the kinematical end-point, where [alpha-,alpha+] is
    //so tiny that S can be assumed to be constant over it, or the initial
    //energy is so extremely tiny that the neutron is basically at rest with
    //undefined direction. It seems appropriate to deal with both cases with a
    //simple isotropic scattering.
    auto alim = getAlphaLimits( m_c_real, beta );
    double alpha = alim.first+rng.generate()*(alim.second-alim.first);
    return ncclamp( alpha,alim.first, alim.second );
  }

  nc_assert( beta >= -m_c*1.001 );
  beta = ncmax( -m_c, beta );
  auto alims = getAlphaLimits( m_c, beta );
  double am = alims.first;
  double ap = alims.second;
  if (am==ap)
    return am;

  const double betasq = beta*beta;
  const double t = betasq * m_Adiv4;
  const double c = 0.0625 * betasq;//NB: this is different c than m_c

  if ( ncmin(t,c) < 1e-5 ) {
    //special case, beta~=0.0. Sample for x=alpha/4A according to
    //exp(-x)/sqrt(x).
    const double fourA = m_Adiv4*16.0;
    const double inv4A = 1.0/fourA;
    const double xxm(am*inv4A), xxp(ap*inv4A);
    NCRYSTAL_DEBUGONLY(unsigned iloop(0));
    while (true) {
      //sample xx from exp(-x)/sqrt(x)
      const double xx = randExpDivSqrt( rng, 1.0, xxm, xxp );
      double alpha = xx * fourA;
      if (alpha<am||alpha>ap)
        continue;

      //Since beta might be slightly different than 0, we must now use the
      //rejection method to correct for the fact that we neglected a factor of
      //exp(-beta^2/(4alpha/A)) = exp(-t/alpha), but normalised to attain unity
      //at alpha=ap, so exp(t/ap-t/alpha).
      //
      // Accept if: R <= std::exp(t/ap-t/alpha)
      //          <=>  logR <= t/ap-t/alpha
      //          <=>  -logR >= -t * (alpha-ap)/(ap*alpha)
      //          <=> (ap*alpha)*(-logR) >= t*(ap-alpha)

      if ( alpha*ap*randExp(rng) >= t*(ap-alpha) )
        return alpha;
      nc_assert(iloop++<100);
    }
  } else {
    //General case, t!=0, c~=0.
    double invt(1.0/t);
    double xm(am*invt), xp(ap*invt);
    double x = randExpMInvXMCXDivSqrtX( rng, c, xm, xp );
    nc_assert(x>=xm*0.999);
    nc_assert(x<=xp*1.001);
    return ncclamp(x*t,am,ap);
  }
}

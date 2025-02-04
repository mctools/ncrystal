#ifndef NCrystal_GaussOnSphere_hh
#define NCrystal_GaussOnSphere_hh

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

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCSpline.hh"

namespace NCRYSTAL_NAMESPACE {

  class GaussOnSphere {
  public:

    //Class providing evaluations, normalisations and integrations of a
    //truncated 2D Gaussian density on a unit sphere (naturally intended to help
    //modelling Bragg diffraction in Gaussian mosaic crystals). Input parameters
    //are the sigma (in radians) of the Gaussian, the angle (in radians) at
    //which the truncation is implemented, and a precision parameter. The
    //density more than the trancation angle away from the Gaussian center will
    //be 0, and the distribution will be properly normalised over the entire
    //sphere.
    //
    //Evaluations can be done either directly via the angles, but also via the
    //cosine of the angles (which often are given directly from e.g. a dot
    //product of unit vectors). In either case, approximations and splined
    //lookup tables and utilised in order to avoid expensive std::exp and
    //std::acos calls.
    //
    //Finally it is possible to integrate the density along a circle on the
    //sphere, and to sample a point on a circle according to the density. The
    //circle integrals will be evaluated either via an efficient Romberg
    //integration scheme, or - where appropriate - via an appropriate
    //approximation formula. The latter once again employ splined lookup tables
    //to avoid expensive calls (this time replacing exp+acos+sqrt+erf calls),
    //and the former uses the CosSinGridGen to avoid expensive cosine
    //calls. Sampling points on circles is implemented precisely as rejection
    //sampling.
    //
    //The precision parameter is used to control the rough precision of the
    //approximations used: When and how to carry out full numerical integration
    //rather than using the approximation formulas, and the number of points in
    //lookup tables. Values can be in the range [1e-7,0.1], with 1e-3 or 1e-4
    //being a sensible default which still provides rather high precision. As a
    //special case for developers, setting the precision parameter to a value in
    //[1,10000] will disable the full numerical integration in circle integrals
    //and will cause lookup tables to contain int(prec+0.5) points (but at least
    //20).
    //
    //If running with the NCRYSTAL_DEBUG_GAUSSONSPHERE environment variable set,
    //the GaussOnSphere destructor will print out various statistics of interest
    //to developers.

    GaussOnSphere();//Constructs invalid object, must call set() before using.
    GaussOnSphere( double sigma, double trunc_angle, double prec = 1e-3 );
    void set(double sigma, double trunc_angle, double prec = 1e-3 );

    ~GaussOnSphere();

    double getNormFactor() const { nc_assert(isValid()); return m_norm; }
    double getSigma() const { nc_assert(isValid()); return m_sigma; }
    double getTruncangle() const { nc_assert(isValid()); return m_truncangle; }
    double getCosTruncangle() const { nc_assert(isValid()); return m_cta; }
    double getSinTruncangle() const { nc_assert(isValid()); return m_sta; }
    double getPrecisionParameter() const { nc_assert(isValid()); return m_prec; }

    //Evaluate the density at x, the angular distance to the center of the Gaussian:
    double evalX(double x) const;//safe everywhere, fast approximate values.
    double evalXInRange(double x) const;//Like evalX, but faster since it incorrectly returns non-zero values
                                        //outside the truncation limit. User code should precheck the range.
    double evalXExact(double x) const;//safe and precise values everywhere (the slowest).

    //For efficiency, it is also possible to evaluate the density directly via
    //the cosine of the angle x, in case the cosine is available (e.g. as result
    //of a dot product of unit vectors or another relatively inexpensive
    //operation). It is undefined behaviour to call these functions without
    //having enabled the lookup table by setting nlt_cosx to a positive value.

    double evalCosX(double cosx) const;//safe for any cosine value,  might be slow.
    double evalCosXInRange(double cosx) const;//like evalXCosX but slightly faster since it assumes cosx>=cos(truncangle) (within
                                              //numerical precision) or something bad will happen (an assert or a bad memory access).

    //Integrated density along a circle which is the result of intersecting the
    //sphere with a cone located an angle of gamma away from the Gaussian center
    //and with an opening angle of alpha.
    double circleIntegral( double cosgamma, double singamma, double cosalpha, double sinalpha ) const;

    //Generate random point on circle, according to the density there. Returns
    //false in case of vanishing density everywhere on circle. The ct=cos(t) and
    //st=sin(t) values can be used to construct the coordinate of the chosen
    //point on the circle, in a frame where the center of the circle is at
    //(0,0,1) and the center of the Gaussian density is at
    //(singamma,0,cosgamma), as (x,y,z)=(sinalpha*ct,sinalpha*st,cosalpha);
    bool genPointOnCircle( RNG&, double cosgamma, double singamma, double cosalpha, double sinalpha,
                           double& cost, double& sint ) const;

    static double calcNormFactor( double sigma, double trunc_angle );

    //Estimate ntrunc by finding the N at which a 2D *planar* Gaussian would
    //contain 1-prec of its volume within a radius of N from it's center:
    static double estimateNTruncFromPrec( double prec, double minval = 3.0, double maxval = 8.0 );

    bool isValid() const { return m_norm>0.0; }

    GaussOnSphere(const GaussOnSphere&) = delete;
    void operator=(const GaussOnSphere&) = delete;

  private:
    double circleIntegralSlow( double cacg, double sasg, double ca, double sa ) const;
    double m_cta;
    double m_circleint_k1;
    double m_circleint_k2;
    double m_norm;
    double m_expfact;
    double m_truncangle;
    double m_sigma;
    double m_numint_accuracy;
    SplinedLookupTable m_lt_sofcosd;
    SplinedLookupTable m_lt_evalcosx;
    double m_prec;//for reference
    double m_sta;//for reference
    //Enable sampling efficiency report when NCRYSTAL_DEBUG_GAUSSONSPHERE is set (only in debug builds):
    void produceStatReport(const char *);
#ifndef NDEBUG
    mutable struct Stats {
      std::atomic<uint64_t> genpointworst = {0};
      std::atomic<uint64_t> genpointcalled = {0};
      std::atomic<uint64_t> genpointtries = {0};
      std::atomic<uint64_t> circleintworst = {0};
      std::atomic<uint64_t> circleintnumber = {0};
      std::atomic<uint64_t> circleintevals = {0};
    } m_stats;
#endif
  };

}


////////////////////////////
// Inline implementations //
////////////////////////////


inline double NCrystal::GaussOnSphere::evalXInRange(double x) const
{
  nc_assert(isValid());
  double res = m_norm * exp_negarg_approx( m_expfact * x * x );
  return res;
}

inline double NCrystal::GaussOnSphere::evalXExact(double x) const
{
  nc_assert(isValid());
  return ncabs(x) <= m_truncangle ? m_norm * std::exp( m_expfact * x * x ) : 0.0;
}

inline double NCrystal::GaussOnSphere::evalX(double x) const {
  nc_assert(isValid());
  return ncabs(x) <= m_truncangle ? evalXInRange(x) : 0.0;
}

inline double NCrystal::GaussOnSphere::evalCosX(double cosx) const {
  nc_assert(isValid());
  if (cosx>=m_lt_evalcosx.getLower())
    return evalCosXInRange(cosx);
  return 0.0;
}

inline double NCrystal::GaussOnSphere::evalCosXInRange(double cosx) const {
  nc_assert(isValid());
  nc_assert(cosx>=m_lt_evalcosx.getLower()-1e-7&&"do not use this method outside the truncation area");
  nc_assert(cosx<=1.0000001);
  double res = ncmax(0.0,m_lt_evalcosx.eval(cosx));
  return res;
}

inline double NCrystal::GaussOnSphere::circleIntegral( double cg, double sg, double ca, double sa ) const
{
  nc_assert(isValid());
  nc_assert(ncabs(cg*cg+sg*sg-1.0)<1e-6);
  nc_assert(ncabs(ca*ca+sa*sa-1.0)<1e-6);
  const double sasg = sa*sg; nc_assert(sasg>=0.0);
  const double cacg = ca*cg;
  const double cd = cacg+sasg;

  if (cd>m_cta&&sasg>=1e-14&&m_circleint_k2 > m_circleint_k1*sasg+cacg) {
    //closed-form approximation gives accurate result here:
    return m_lt_sofcosd.eval(cd)*std::sqrt(sa/sg);
  }

  //special case, needs more careful treatment:
  return circleIntegralSlow( cg,sg,ca,sa );
  //NB: Passing in sasc and cacg to circleIntegralSlow seems like it might save
  //a few multiplications, but it actually slows down in certain benchmarks
  //(because it breaks tail calls or inlining somehow?).
}

#endif

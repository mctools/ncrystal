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

#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCRotMatrix.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include <sstream>
#include <list>
#include <queue>

namespace NC = NCrystal;

bool NC::nc_is_grid( Span<const double> v )
{
  if ( v.size() < 2 )
    return false;
  double last = v.front();
  if ( !std::isfinite(last) )
    return false;
  for ( auto e : v.subspan(1) ) {
    if ( !(std::isfinite(e) && e>last ) )
      return false;
    last = e;
  }
  return true;
}

NC::VectD NC::logspace(double start, double stop, unsigned num)
{
  nc_assert(num>1);
  nc_assert(stop>start);
  VectD vec(num) ;
  double interval = (stop-start)/(num-1);
  for(VectD::iterator it=vec.begin();it!=vec.end();++it)  {
    *it = std::pow(10.0,start);
    start += interval;
  }
  vec.back() = std::pow(10.0,stop);
  return vec;
}

NC::VectD NC::geomspace(double start, double stop, unsigned num)
{
  nc_assert(num>1);
  nc_assert(stop>start);
  auto v = logspace(std::log10(start),std::log10(stop),num);
  v.front() = start;
  v.back() = stop;
  return v;
}

NC::VectD NC::linspace(double start, double stop, unsigned num)
{
  nc_assert(num>1);
  nc_assert(stop>start);
  VectD v;
  v.reserve(num) ;
  unsigned num_minus_1 = num - 1;
  const double interval = (stop-start)/num_minus_1;
  for (unsigned i = 0; i<num_minus_1;++i)
    v.push_back(std::fma(static_cast<double>(i),interval,start));
  v.push_back( stop );
  return v;
}

NC::VectD NC::powspace(double a, double b, unsigned num, double p )
{
  nc_assert(num >= 2);
  nc_assert(a > 0.0);
  nc_assert(b > a);
  nc_assert(p > 0.0);
  nc_assert(std::isfinite(a));
  nc_assert(std::isfinite(b));
  nc_assert(std::isfinite(p));
  nc_assert(num < 1000000000ULL);

  const double step = 1.0 / static_cast<double>(num - 1);
  const double delta = b - a;
  const double nm1 = static_cast<double>(num - 1);

  VectD res;
  if ( p==1.0 ) {
    res = linspace(a, b, num);
    return res;
  }

  res.reserve(num);
  res.push_back(a);

  //NB: 2.0, 1.5, 4.0 are used in our code. We can add others as needed.
  if (p == 2.0) {
    for (double i = 1.0; i < nm1; i += 1.0) {
      const double s = i * step;
      res.push_back(std::fma(delta,ncsquare(s),a));
    }
  } else if (p == 3.0) {
    for (double i = 1.0; i < nm1; i += 1.0) {
      const double s = i * step;
      res.push_back(std::fma(delta*ncsquare(s),s,a));
    }
  } else if (p == 4.0) {
    for (double i = 1.0; i < nm1; i += 1.0) {
      const double s = i * step;
      res.push_back(std::fma(delta,ncsquare(ncsquare(s)),a));
    }
  } else if (p == 1.5) {
    for (double i = 1.0; i < nm1; i += 1.0) {
      const double s = i * step;
      res.push_back(std::fma(delta,s * std::sqrt(s),a));
    }
  } else if (p == 0.5) {
    for (double i = 1.0; i < nm1; i += 1.0) {
      const double s = i * step;
      res.push_back(std::fma(delta,std::sqrt(s),a));
    }
  } else {
    for (double i = 1.0; i < nm1; i += 1.0) {
      const double s = i * step;
      res.push_back(std::fma(delta,std::pow(s, p),a));
    }
  }

  res.push_back(b);
  nc_assert( res.size() == num );
  return res;
}

bool NC::isPrime(unsigned n) {
  if (n>3) {
    if ( !(n%2) || !(n%3) )
      return false;//fast precheck for factors of 2 and 3
    for (unsigned i = 5; i*i <= n; i += 6 ) {
      if ( !(n % i) || !(n%(i + 2)) )
        return false;
    }
    return true;
  }
  return n>1;//2+3:prime, 0+1: not prime
}

void NC::sincos_mpi2pi2(double A, double& cosA, double& sinA) {
  nc_assert(ncabs(A)<=kPiHalf);
  //Evaluate at A/2 via Taylor expansions and get final results via
  //double-angle formula.
  double x = 0.5*A;
  double mx2 = -x*x;
  double s2 =  x*(1.0 + mx2 * ( 1.66666666666666666666666666666666666666666667e-1 // - x^3 / 3!
                      + mx2 * ( 8.33333333333333333333333333333333333333333333e-3 // + x^5 / 5!
                      + mx2 * ( 1.98412698412698412698412698412698412698412698e-4 // - x^7 / 7!
                      + mx2 * ( 2.75573192239858906525573192239858906525573192e-6 // + x^9 / 9!
                      + mx2 * ( 2.50521083854417187750521083854417187750521084e-8 // - x^11 / 11!
                      + mx2 * ( 1.60590438368216145993923771701549479327257105e-10 // + x^13 / 13!
                      + mx2 * ( 7.64716373181981647590113198578807044415510024e-13 // - x^15 / 15!
                                ))))))));
  double c2m1 = mx2 * ( 0.5 // - x^2 / 2!
              + mx2 * ( 4.16666666666666666666666666666666666666666667e-2 // + x^4 / 4!
              + mx2 * ( 1.38888888888888888888888888888888888888888889e-3 // - x^6 / 6!
              + mx2 * ( 2.48015873015873015873015873015873015873015873e-5 // + x^8 / 8!
              + mx2 * ( 2.75573192239858906525573192239858906525573192e-7 // - x^10 / 10!
              + mx2 * ( 2.08767569878680989792100903212014323125434237e-9 // + x^12 / 12!
              + mx2 * ( 1.14707455977297247138516979786821056662326504e-11 // - x^14 / 14!
              + mx2 * ( 4.77947733238738529743820749111754402759693765e-14 // + x^16 / 16!
                        ))))))));
  double k = 2.0*c2m1;
  sinA = (k+2.0)*s2;
  cosA = 1.0+k*(c2m1+2.0);
}

void NC::sincos_mpi8pi8(double A, double& cosA, double& sinA) {
  nc_assert(ncabs(A)<=kPi*0.12500001);
  //Evaluate at A/2 via two 6 terms Taylor expansions and get final
  //results via double-angle formula.
  double x = 0.5*A;
  double mx2 = -x*x;
  double s2 = x*(1.0 + mx2 * ( 1.66666666666666666666666666666666666666666667e-1 // - x^3 / 3!
                     + mx2 * ( 8.33333333333333333333333333333333333333333333e-3 // + x^5 / 5!
                     + mx2 * ( 1.98412698412698412698412698412698412698412698e-4 // - x^7 / 7!
                     + mx2 * ( 2.75573192239858906525573192239858906525573192e-6 // + x^9 / 9!
                     + mx2 * ( 2.50521083854417187750521083854417187750521084e-8 // - x^11 / 11!
                               ))))));
  //Taylor expansion for cosine is short enough that it beats
  //c2=std::sqrt(1-s2*s2). Evaluate without 1.0 term for numerical stability:
  double c2m1 = mx2 * ( 0.5 // - x^2 / 2!
              + mx2 * ( 4.16666666666666666666666666666666666666666667e-2 // + x^4 / 4!
              + mx2 * ( 1.38888888888888888888888888888888888888888889e-3 // - x^6 / 6!
              + mx2 * ( 2.48015873015873015873015873015873015873015873e-5 // + x^8 / 8!
              + mx2 * ( 2.75573192239858906525573192239858906525573192e-7 // - x^10 / 10!
                        )))));
  double k = 2.0*c2m1;
  sinA = (k+2.0)*s2;
  cosA = 1.0+k*(c2m1+2.0);

}

double NC::cos_mpipi(double A)
{
  //Use abs/min/copysign tricks to actually do the evaluation in [-pi/2,pi/2].
  double Aabs = ncabs(A);
  nc_assert(Aabs<=kPi);

  //Taylor expansion to 22nd order
  double x = ncmin(Aabs,kPi-Aabs);
  double mx2 = -x*x;
  double c = 1.0 + mx2 * ( 0.5 // - x^2 / 2!
             + mx2 * ( 4.16666666666666666666666666666666666666666667e-2 // + x^4 / 4!
             + mx2 * ( 1.38888888888888888888888888888888888888888889e-3 // - x^6 / 6!
             + mx2 * ( 2.48015873015873015873015873015873015873015873e-5 // + x^8 / 8!
             + mx2 * ( 2.75573192239858906525573192239858906525573192e-7 // - x^10 / 10!
             + mx2 * ( 2.08767569878680989792100903212014323125434237e-9 // + x^12 / 12!
             + mx2 * ( 1.14707455977297247138516979786821056662326504e-11 // - x^14 / 14!
             + mx2 * ( 4.77947733238738529743820749111754402759693765e-14 // + x^16 / 16!
             + mx2 * ( 1.56192069685862264622163643500573334235194041e-16 // - x^18 / 18!
             + mx2 * ( 4.1103176233121648584779906184361403746103695e-19 // + x^20 / 20!
             + mx2 * ( 8.89679139245057328674889744250246834331248809e-22 // - x^22 / 22!
              )))))))))));
  return std::copysign(c,kPiHalf-Aabs);
}

double NC::cos_mpi2pi2(double x)
{
  nc_assert(ncabs(x)<=kPiHalf);
  //Taylor expansion to 22nd order, precision better than 1.1e-16 over entire range
  double mx2 = -x*x;
  return 1.0 + mx2 * ( 0.5 // - x^2 / 2!
             + mx2 * ( 4.16666666666666666666666666666666666666666667e-2 // + x^4 / 4!
             + mx2 * ( 1.38888888888888888888888888888888888888888889e-3 // - x^6 / 6!
             + mx2 * ( 2.48015873015873015873015873015873015873015873e-5 // + x^8 / 8!
             + mx2 * ( 2.75573192239858906525573192239858906525573192e-7 // - x^10 / 10!
             + mx2 * ( 2.08767569878680989792100903212014323125434237e-9 // + x^12 / 12!
             + mx2 * ( 1.14707455977297247138516979786821056662326504e-11 // - x^14 / 14!
             + mx2 * ( 4.77947733238738529743820749111754402759693765e-14 // + x^16 / 16!
             + mx2 * ( 1.56192069685862264622163643500573334235194041e-16 // - x^18 / 18!
             + mx2 * ( 4.1103176233121648584779906184361403746103695e-19 // + x^20 / 20!
             + mx2 * ( 8.89679139245057328674889744250246834331248809e-22 // - x^22 / 22!
              )))))))))));
}

double NC::cos_mpi8pi8(double x)
{
  nc_assert(ncabs(x)<=0.125000000001*kPi);
  //Taylor expansion to 12th order, precision better than 3.3e-16 over entire range
  double mx2 = -x*x;
  return 1.0 + mx2 * ( 5.0e-1 // - x^2 / 2!
             + mx2 * ( 4.16666666666666666666666666666666666666666667e-2 // + x^4 / 4!
             + mx2 * ( 1.38888888888888888888888888888888888888888889e-3 // - x^6 / 6!
             + mx2 * ( 2.48015873015873015873015873015873015873015873e-5 // + x^8 / 8!
             + mx2 * ( 2.75573192239858906525573192239858906525573192e-7 // - x^10 / 10!
             + mx2 * ( 2.08767569878680989792100903212014323125434237e-9 // + x^12 / 12!
              ))))));
}

double NC::sin_mpipi(double A)
{
  //Use abs/min/copysign tricks to actually do the evaluation in [-pi/2,pi/2].
  double Aabs = ncabs(A);
  nc_assert(Aabs<=kPi);

  //Taylor expansion to 19th order
  double x = ncmin(Aabs,kPi-Aabs);
  double mx2 = -x*x;
  double s = x*(1.0 + mx2 * ( 1.66666666666666666666666666666666666666666667e-1 // - x^3 / 3!
                + mx2 * ( 8.33333333333333333333333333333333333333333333e-3 // + x^5 / 5!
                + mx2 * ( 1.98412698412698412698412698412698412698412698e-4 // - x^7 / 7!
                + mx2 * ( 2.75573192239858906525573192239858906525573192e-6 // + x^9 / 9!
                + mx2 * ( 2.50521083854417187750521083854417187750521084e-8 // - x^11 / 11!
                + mx2 * ( 1.60590438368216145993923771701549479327257105e-10 // + x^13 / 13!
                + mx2 * ( 7.64716373181981647590113198578807044415510024e-13 // - x^15 / 15!
                + mx2 * ( 2.81145725434552076319894558301032001623349274e-15 // + x^17 / 17!
                + mx2 * ( 8.22063524662432971695598123687228074922073899e-18 // - x^19 / 19!
              ))))))))));
  return std::copysign(s,A);
}

double NC::sin_mpi2pi2(double x)
{
  nc_assert(ncabs(x)<=kPiHalf);
  //Taylor expansion to 19th order, precision better than 6e-17 over entire range
  double mx2 = -x*x;
  return x*(1.0 + mx2 * ( 1.66666666666666666666666666666666666666666667e-1 // - x^3 / 3!
                + mx2 * ( 8.33333333333333333333333333333333333333333333e-3 // + x^5 / 5!
                + mx2 * ( 1.98412698412698412698412698412698412698412698e-4 // - x^7 / 7!
                + mx2 * ( 2.75573192239858906525573192239858906525573192e-6 // + x^9 / 9!
                + mx2 * ( 2.50521083854417187750521083854417187750521084e-8 // - x^11 / 11!
                + mx2 * ( 1.60590438368216145993923771701549479327257105e-10 // + x^13 / 13!
                + mx2 * ( 7.64716373181981647590113198578807044415510024e-13 // - x^15 / 15!
                + mx2 * ( 2.81145725434552076319894558301032001623349274e-15 // + x^17 / 17!
                + mx2 * ( 8.22063524662432971695598123687228074922073899e-18 // - x^19 / 19!
              ))))))))));
}

double NC::sin_mpi8pi8(double x)
{
  nc_assert(ncabs(x)<=0.125000000001*kPi);
  //Taylor expansion to 13th order, precision better than 6e-17 over entire range
  double mx2 = -x*x;
  return x*(1.0 + mx2 * ( 1.66666666666666666666666666666666666666666667e-1 // - x^3 / 3!
                + mx2 * ( 8.33333333333333333333333333333333333333333333e-3 // + x^5 / 5!
                + mx2 * ( 1.98412698412698412698412698412698412698412698e-4 // - x^7 / 7!
                + mx2 * ( 2.75573192239858906525573192239858906525573192e-6 // + x^9 / 9!
                + mx2 * ( 2.50521083854417187750521083854417187750521084e-8 // - x^11 / 11!
                + mx2 * ( 1.60590438368216145993923771701549479327257105e-10 // + x^13 / 13!
              )))))));
}

double NC::estimateDerivative(const Fct1D* f, double x, double h, unsigned order)
{
  nc_assert(f);
  nc_assert(h>0);
  nc_assert(order==4||order==6);
  if (order==6)
    return (  256.*f->eval(x+0.25*h)-256.*f->eval(x-0.25*h)-32.*f->eval(x+0.5*h)+32.*f->eval(x-0.5*h)
              -8.*f->eval(x+0.5*h)+8.*f->eval(x-0.5*h)+f->eval(x+h)-f->eval(x-h) ) / (90.*h);
  return (-f->eval(x+h)+8*f->eval(x+0.5*h)-8*f->eval(x-0.5*h)+f->eval(x-h))/(6.0*h);
}

double NC::estimateSingleSidedDerivative(const Fct1D* f, double x, double h, unsigned order)
{
  nc_assert(f);
  nc_assert(h!=0);
  nc_assert(order>=1&&order<=4);

  switch(order) {
  case 1:
    return ( f->eval(x+h)-f->eval(x) ) / h;
  case 2:
    return ( -3.0*f->eval(x) +4.0*f->eval(x+0.5*h) -f->eval(x+h)) / h;
  case 3:
    return -2.*(4.5*f->eval(x)-8.*f->eval(x+0.25*h)+4.*f->eval(x+0.5*h)-0.5*f->eval(x+h)) /h;
  case 4:
    return -(2/3.)*(31.5*f->eval(x)-64.*f->eval(x+0.125*h)+40.*f->eval(x+0.25*h)-8.*f->eval(x+0.5*h)+0.5*f->eval(x+h))/h;
  default:
    nc_assert_always(false);
  }
}

double NC::findRoot(const Fct1D*f,double a, double b, double acc)
{
  //Basically a mix between a binary search (bisection) algorithm and a "false
  //position" algorithm, trying to get both the stability of binary search and
  //the improved convergence rate of "false position".
  nc_assert(f);
  double fa = f->eval(a);
  double fb = f->eval(b);
  if (!(b > a))
    NCRYSTAL_THROW(CalcError,"root finding requires b>a.");
  if (fa == 0.0)
    return a;
  if (fb == 0.0)
    return b;
  if ( (fa < 0.0) == (fb < 0.0) )
    NCRYSTAL_THROW(CalcError,"root finding requires f(a) * f(b) <= 0");
  acc *= 0.5;//safety
  unsigned i(60);
  while(--i) {
    //At a point c inside [a,b]. In the pure bisection method this would be
    //b=0.5*(a+b) and in the pure false position method this would be c =
    //(a*fb-b*fa)/(fb-fa). The problem with the latter is that is can sometime
    //be too near the edges of the interval, even though the root is further
    //inside, thus resulting for a very slow convergence for certain functions.
    //
    //Our ad hoc combination here uses c from the false position method, but
    //constrained so that the next step will always split the interval at a
    //point which is at least 15% from the edges.
    double dfba = fb-fa;
    nc_assert(dfba);
    double c = (a*fb-b*fa)/dfba;
    if ( b-a<acc )
      return c;
    double k = 0.15*(b-a);
    c = ncmax(a+k,ncmin(b-k,c));
    double fc = f->eval(c);
    if ( !fc )
      return c;
    if (fa*fc<0) {
      //root must be in [a,c]
      b=c;fb=fc;
    } else {
      //root must be in [c,b]
      a=c;fa=fc;
    }
  }
  NCRYSTAL_THROW(CalcError,"Root search failed to converge!");
}

NC::Fct1D::~Fct1D(){}

double NC::stable_exp( double x )
{
  const double y0 = std::exp(x);
  if ( y0 == 0.0 || !std::isfinite(y0) )
    return y0;//underflow, overflow, or nan-in-nan-out: nothing to refine.
  //One Newton-Raphson step on f(y)=ln(y)-x=0 (f'(y)=1/y):
  //y1 = y0 - (ln(y0)-x)*y0 = y0+y0*(x-ln(y0)). Use std::fma for the final
  //combination (single rounding), rather than y0*(1+correction) (two
  //roundings: forming 1+correction, then the product):
  const double correction = x - std::log(y0);
  return std::fma( y0, correction, y0 );
}

double NC::stable_expm1( double x )
{
  const double y0 = std::expm1(x);
  if ( !std::isfinite(y0) )
    return y0;//overflow, or nan-in-nan-out.
  const double onepy0 = 1.0 + y0;
  if ( onepy0 == 0.0 )
    return y0;//y0=-1 exactly (x very negative): nothing left to refine.
  //One Newton-Raphson step on f(y)=log1p(y)-x=0 (f'(y)=1/(1+y)):
  //y1 = y0 - (log1p(y0)-x)*(1+y0) = y0 + (x-log1p(y0))*(1+y0). Use
  //std::fma for the final combination (single rounding):
  const double correction = x - std::log1p(y0);
  return std::fma( correction, onepy0, y0 );
}

double NC::stable_log( double x )
{
  const double y0 = std::log(x);
  if ( !std::isfinite(y0) )
    return y0;//x<=0 (-inf/nan) or x=+inf: nothing to refine.
  //One Newton-Raphson step on f(y)=exp(y)-x=0 (f'(y)=exp(y)):
  //y1 = y0 - (exp(y0)-x)/exp(y0) = y0 + (x-exp(y0))/exp(y0). The
  //subtraction x-exp(y0) is close to exact (Sterbenz's lemma: exp(y0) is
  //within a handful of ULP of x by construction, well within the 2x
  //factor Sterbenz's lemma requires for an exact FP subtraction), so no
  //separate stabilisation is needed there -- just a plain division and
  //sum (no std::fma: the correction here has unit scale, so a plain sum
  //is already a single rounding, unlike stable_exp's multiply-add):
  const double expy0 = std::exp(y0);
  if ( !( expy0 > 0.0 ) || !std::isfinite(expy0) )
    return y0;//underflow/overflow reconstructing exp(y0): nothing to refine.
  const double correction = (x - expy0)/expy0;
  return y0 + correction;
}

double NC::stable_tanh( double x )
{
  //tanh(x) = (e^(2x)-1)/(e^(2x)+1) = expm1(2x)/(expm1(2x)+2). Well
  //conditioned for all x (no cancellation near x=0, since stable_expm1
  //already handles that), and inherits stable_expm1's refined accuracy
  //directly rather than needing its own Newton-Raphson correction (which
  //would be ill-conditioned for tanh: atanh'(y)->infinity as y->+-1):
  if ( ncisnan(x) )
    return x;
  const double t = stable_expm1(2.0*x);
  if ( t == kInfinity )
    return 1.0;//avoids inf*inf/inf further down for large positive x.
  return t/(t+2.0);
}

double NC::stable_sinh( double x )
{
  //sinh(x) = t*(t+2)/(2*(1+t)) for t=expm1(x)=e^x-1, since
  //t*(t+2)/(1+t) = (e^x-1)(e^x+1)/e^x = e^x-e^-x = 2*sinh(x). Well
  //conditioned for x>=0 (then t>=0, so 1+t>=1), but the (1+t) denominator
  //goes to 0 as x->-infinity (t->-1), amplifying t's own tiny relative
  //error into a large one -- confirmed empirically (e.g. off by ~2e-8
  //relative for x=-20, far more than a rounding-level error). Sidestep
  //this entirely using the odd-function symmetry sinh(-x)=-sinh(x) (an
  //exact identity: negation is a free, lossless sign-bit flip), so the
  //division is only ever evaluated in the well-conditioned x>=0 branch:
  if ( ncisnan(x) )
    return x;
  if ( x < 0.0 )
    return -stable_sinh(-x);
  const double t = stable_expm1(x);
  if ( t == kInfinity )
    return kInfinity;//avoids inf*inf/inf further down for large positive x.
  return 0.5*t*(t+2.0)/(1.0+t);
}

namespace NCRYSTAL_NAMESPACE {
  double erfcdiff_notaylor(double a, double b)
  {
    nc_assert(b>=a);
    if (b<0) {
      //both numbers are negative, more precise to evaluate with both numbers
      //positive where erfc~=0 instead of erfc~=2:
      //We use erfcdiff(a,b)=erfcdiff(-b,-a):
      b=-b;
      a=-a;
      std::swap(a,b);
      //NB: b>=a still holds now
    }
    //ok, b>a and at least one number is positive:
    nc_assert( b>=a && b>=0 );
    //evaluate (remembering std::erfc(x)=0 for x above ~27.3 or so:
    const double erfca = a>27.3 ? 0.0 : std::erfc(a);
    if ( b > a+4.0 && ( a>=4 || ( a < 0.0 && b > 6.0 ) ) ) {
      //erfc(b) contribution negligible at double precision
      return erfca;
    }
    const double erfcb = b>27.3 ? 0.0 : std::erfc(b);
    return erfca - erfcb;
  }
}

double NC::erfcdiff(double a, double b)
{
  if ( ncmax(ncabs(a),ncabs(b)) < 0.32 ) {
    //Both arguments are small (which will happen in free gas scattering for
    //neutrons at low energies). Evaluate via taylor expansions (important
    //also for numerical stability since erfc(x)~=1+O(x) for tiny arguments,
    //and the 1 cancels out in erfc(a)-erfc(b).
    //
    //Use expansion: erfc(x)-1 ~= c1*x+c3*x^3+...+c11*x^11+O(x^13)
    //
    constexpr double c1  = - 2.0 * kInvSqrtPi;
    constexpr double c3  =   2.0 * kInvSqrtPi / 3.0;
    constexpr double c5  = - 0.2 * kInvSqrtPi;
    constexpr double c7  =   kInvSqrtPi / 21.0;
    constexpr double c9  = - kInvSqrtPi / 108.0;
    constexpr double c11 =   kInvSqrtPi / 660.0;
    constexpr double c13 = - kInvSqrtPi / 4680.0;
    constexpr double c15 =   kInvSqrtPi / 37800.0;
    const double a2 = a*a;
    const double b2 = b*b;
    const double a3to11 = a * a2 * ( c3 + a2 * ( c5 + a2 * ( c7 + ( a2 * ( c9 + a2 * ( c11 + a2 * ( c13 + a2 * c15 ) ) ) ) ) ) );
    const double b3to11 = b * b2 * ( c3 + b2 * ( c5 + b2 * ( c7 + ( b2 * ( c9 + b2 * ( c11 + b2 * ( c13 + b2 * c15 ) ) ) ) ) ) );
    return c1*(a-b) + ( a3to11 - b3to11 );
  }

  //Use erfcdiff_notaylor, ordering arguments so b>=a.
  return a>b ? -erfcdiff_notaylor(b,a) : erfcdiff_notaylor(a,b);
}


double NC::erfc_rescaled(double x, double b)
{
  //exp(b)*erfc(x), but faster and more precise.
  nc_assert(x>=0.0);
  if (b<-745.1)
    return 0.0;//erfc <= 1, so exp(b) will always force strictly 0 here.
  if ( ( x<23.0 && ncabs(b)<700 ) || x < 5 ) {
    //standard functions provide full precision here (at least when |b|<700).
    return std::exp(b)*std::erfc(x);
  }
  //large x, employ expansion and combine exp(b)*exp(-x^2)=exp(b-x^2). If the
  //caller picked b appropriately, so |b-x^2|<700, this provides full precision.
  const double bxx=b-x*x;
  if (bxx<-745.1)
    return 0.0;
  const double c3  = -0.5         ;// = 1/2
  const double c5  =  0.75        ;// = 3/4
  const double c7  = -1.875       ;// = 15/8
  const double c9  =  6.5625      ;// = 105/16
  const double c11 = -29.53125    ;// = 945/32
  //const double c13 = 162.421875 ;// = 10395/64
  const double y = 1/x;
  const double y2 = y*y;
  return kInvSqrtPi*std::exp(bxx)*(y+y2*(c3+y2*(c5+y2*(c7+y2*(c9+y2*c11)))));
}

std::pair<NC::VectD, NC::VectD>
NC::reducePtsInDistribution(Span<const double> x,
                            Span<const double> y,
                            std::size_t targetN,
                            const PtReduceCfg& cfg)
{
  // Strategi to ensure good performance: Perform a coarse pre-thinning, removing
  // approximately 25% of points per pass.  The final reduction uses a priority
  // queue and local importance updates.

#ifndef NDEBUG
  {
    const std::size_t nn = x.size();
    nc_assert(nn == y.size());
    nc_assert(nn >= 2);
    nc_assert(targetN >= 2);
    nc_assert(std::isfinite(cfg.equidistant_fraction));
    nc_assert(cfg.equidistant_fraction >= 0.0);
    nc_assert(std::isfinite(cfg.tail_floor));
    nc_assert(cfg.tail_floor > 0.0);
    nc_assert(nc_is_grid(x));
    for (std::size_t i = 0; i < nn; ++i) {
      nc_assert(std::isfinite(y[i]));
      nc_assert(y[i]>=0.0);
    }
  }
#endif

  const std::size_t NONE = std::numeric_limits<std::size_t>::max();

  if (targetN >= x.size())
    return { VectD{x.begin(), x.end()}, VectD{y.begin(), y.end()} };

  const double ymax = *std::max_element(y.begin(), y.end());
  nc_assert(ymax > 0.0);

  const double invYmax = 1.0 / ymax;
  const double gappts = cfg.equidistant_fraction * targetN - 1.0;

  const double maxGap =
    gappts > 0.0 ? (x.back() - x.front()) / gappts : kInfinity;

  struct Candidate {
    double importance;
    std::size_t i;
  };

  VectD px(x.begin(), x.end());
  VectD py(y.begin(), y.end());
  VectD plny;

  auto rebuildLog = [&]() {
    const std::size_t ncur = px.size();

    plny.resize(ncur);

    for (std::size_t i = 0; i < ncur; ++i) {
      vectAt(plny, i) =
        std::log(ncmax(cfg.tail_floor,
                       vectAt(py, i) * invYmax));
    }
  };

  auto calcImportance = [&](std::size_t i0,
                            std::size_t i1,
                            std::size_t i2) -> double {
    const double gap = vectAt(px, i2) - vectAt(px, i0);

    if (gap > maxGap) {
      //Instead of just returning infinity, we return something that scales with
      //the gap size. Just in case weird settings results in a too small
      //max_gap. In that case we at least should prevent bigger gaps before
      //smaller gaps.
      return 1e250 * ncclamp(gap, 1e-305, 1e55);
    }

    //Area of triangle with 3 points at the corners is the area the curve
    //integral will change if we remove the point (x1,y1). Apart from a missing
    //factor of 0.5, the change in area is thus:
    const double dx1 = vectAt(px, i1) - vectAt(px, i0);
    const double dx2 = vectAt(px, i2) - vectAt(px, i0);
    const double area =
      ncabs(dx1 * (vectAt(py, i2) - vectAt(py, i0)) -
            dx2 * (vectAt(py, i1) - vectAt(py, i0)));

    //To preserve also small features in the tails, we find the equivalent
    //change in the logy-curve, and combine the two for the final importance. We
    //did consider other options of combining the two factors, but the one below
    //seemed to work best for our use-cases related to VDOS expansion.

    const double logArea =
      ncabs(dx1 * (vectAt(plny, i2) - vectAt(plny, i0)) -
            dx2 * (vectAt(plny, i1) - vectAt(plny, i0)));
    return area * logArea * logArea;
  };

  //Coarse prefilter: Crucially, we are only considering alternating points,
  //which means importance scores of all considered points are stable even if
  //some are removed.

  std::vector<Candidate> cand;
  while (px.size()/targetN > 8 ) {//fixme: revisit!
    rebuildLog();//FIXME: Can we avoid this?
    const std::size_t ncur = px.size();
    cand.clear();
    cand.reserve((ncur - 2) / 2);
    for (std::size_t i = 1; i + 1 < ncur; i += 2)
      cand.push_back( { calcImportance(i - 1, i, i + 1), i });

    std::sort( cand.begin(), cand.end(),
               [](const Candidate& a, const Candidate& b) {
                 if (a.importance != b.importance)
                   return a.importance < b.importance;
                 return a.i < b.i;
               });

    const std::size_t ndel = cand.size() / 2;
    std::vector<unsigned char> remove(ncur, 0);

    for (std::size_t j = 0; j < ndel; ++j)
      vectAt(remove, cand[j].i) = 1;

    VectD nx, ny;
    nx.reserve(ncur - ndel);
    ny.reserve(ncur - ndel);
    for (std::size_t i = 0; i < ncur; ++i) {
      if (!vectAt(remove, i)) {
        nx.push_back(vectAt(px, i));
        ny.push_back(vectAt(py, i));
      }
    }

    nc_assert( nx.size() == ncur-ndel );
    nc_assert( ny.size() == ncur-ndel );
    px.swap(nx);
    py.swap(ny);
  }
  rebuildLog();

  const std::size_t n = px.size();

  struct Entry {
    double importance;
    std::size_t i, gen;
  };

  struct EntryLess {
    bool operator()(const Entry& a, const Entry& b) const
    {
      if (a.importance != b.importance)
        return a.importance > b.importance;
      return a.i > b.i;
    }
  };

  using Queue = std::priority_queue<Entry, std::vector<Entry>, EntryLess>;
  std::vector<std::size_t> prev(n), next(n), gen(n, 0);

  for (std::size_t i = 0; i < n; ++i) {
    vectAt(prev, i) = i == 0 ? NONE : i - 1;
    vectAt(next, i) = i + 1 == n ? NONE : i + 1;
  }

  Queue q;

  auto enqueue = [&](std::size_t i) {
    if (i == NONE)
      return;
    if (vectAt(prev, i) == NONE ||
        vectAt(next, i) == NONE)
      return;
    ++vectAt(gen, i);
    double importance = calcImportance(vectAt(prev,i), i, vectAt(next,i));
    q.push( { importance, i, vectAt(gen,i) });
  };

  for (std::size_t i = 1; i + 1 < n; ++i)
    enqueue(i);

  std::size_t left = n;
  while (left > targetN) {
    Entry e;
    for (;;) {
      nc_assert(!q.empty());
      e = q.top();
      q.pop();
      const std::size_t i = e.i;
      if (vectAt(prev, i) != NONE &&
          vectAt(next, i) != NONE &&
          vectAt(gen, i) == e.gen)
        break;
    }
    const std::size_t i = e.i;
    const std::size_t i0 = vectAt(prev, i);
    const std::size_t i2 = vectAt(next, i);
    vectAt(next, i0) = i2;
    vectAt(prev, i2) = i0;
    vectAt(prev, i) = NONE;
    vectAt(next, i) = NONE;
    ++vectAt(gen, i);
    --left;
    enqueue(i0);
    enqueue(i2);
  }

  VectD nx, ny;
  nx.reserve(targetN);
  ny.reserve(targetN);
  std::size_t i = 0;
  while (i != NONE) {
    nx.push_back(vectAt(px, i));
    ny.push_back(vectAt(py, i));
    i = vectAt(next, i);
  }

  return {std::move(nx), std::move(ny)};
}

std::pair<NC::VectD, NC::VectD>
NC::reducePtsByEquidistribution(Span<const double> x,
                                Span<const double> y,
                                std::size_t targetN,
                                const PtReduceCfg& cfg)
{
  // Inputs: increasing x, non-negative y (same size >=2), number of points
  // requested (>=2), and cfg (see comments in header). Returns targetN of the
  // input points (or all of them if fewer than targetN).

#ifndef NDEBUG
  {
    nc_assert(x.size() == y.size());
    nc_assert(x.size() >= 2);
    nc_assert(targetN >= 2);
    nc_assert(std::isfinite(cfg.equidistant_fraction));
    nc_assert(cfg.equidistant_fraction >= 0.0);
    nc_assert(cfg.equidistant_fraction <= 1.0);
    nc_assert(std::isfinite(cfg.tail_floor));
    nc_assert(cfg.tail_floor > 0.0);
    nc_assert(nc_is_grid(x));
    for (std::size_t i = 0; i < y.size(); ++i) {
      nc_assert(std::isfinite(y[i]));
      nc_assert(y[i] >= 0.0);
    }
  }
#endif

  const std::size_t n = x.size();
  if (targetN >= n)
    return { VectD{x.begin(), x.end()}, VectD{y.begin(), y.end()} };

  const double ymax = *std::max_element(y.begin(), y.end());

  //Densities are piecewise constant on the n-1 intervals between the input
  //points, so we need n-1 interval widths:
  VectD width(n - 1);
  StableSumKahan xrange;
  for (std::size_t i = 0; i + 1 < n; ++i) {
    vectAt(width, i) = x[i + 1] - x[i];
    xrange.add(vectAt(width, i));
  }
  const double invL = 1.0 / xrange.sum();

  //Density based on |f''| of the function values in f (which are of order 1
  //or below). First we estimate |f''| at each point from the change of slope,
  //and then set the density in each interval as the square root of the average
  //of the |f''| at its endpoints. The result is normalised to have a total
  //integral of 1, unless it is (essentially) 0, in which case the returned
  //vector is empty.
  //
  //Since the normalisation would otherwise turn rounding errors (for a
  //function with no curvature, like a straight line) into a density with full
  //weight, curvatures which are not significantly above what values with a
  //relative accuracy of 1e-12 can resolve at the given spacing (an f'' of
  //1e-12/width^2) are considered to be 0.
  constexpr double relative_noise_level = 1e-12;
  VectD kappa(n), dens(n - 1);
  auto sqrtCurvatureDensity = [&](const VectD& f) -> VectD
  {
    for (std::size_t i = 1; i + 1 < n; ++i) {
      const double s0 = (vectAt(f, i) - vectAt(f, i - 1)) / vectAt(width, i - 1);
      const double s1 = (vectAt(f, i + 1) - vectAt(f, i)) / vectAt(width, i);
      const double wmean = 0.5 * (vectAt(width, i - 1) + vectAt(width, i));
      const double k = ncabs(s1 - s0) / wmean;
      const double noise = relative_noise_level / (wmean * wmean);
      vectAt(kappa, i) = ( k > noise ? k - noise : 0.0 );
    }
    kappa.front() = vectAt(kappa, 1);
    kappa.back() = vectAt(kappa, n - 2);
    StableSumKahan total;
    for (std::size_t i = 0; i + 1 < n; ++i) {
      vectAt(dens, i) = std::sqrt(0.5 * (vectAt(kappa, i)
                                         + vectAt(kappa, i + 1)));
      total.add(vectAt(dens, i) * vectAt(width, i));
    }
    const double t = total.sum();
    if (!(t > 0.0) || !std::isfinite(t))
      return {};
    const double invt = 1.0 / t;
    VectD res(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i)
      vectAt(res, i) = vectAt(dens, i) * invt;
    return res;
  };

  VectD dlin, dlog;
  if (ymax > 0.0) {
    VectD f(n);
    const double invymax = 1.0 / ymax;
    for (std::size_t i = 0; i < n; ++i)
      vectAt(f, i) = y[i] * invymax;
    dlin = sqrtCurvatureDensity(f);
    for (std::size_t i = 0; i < n; ++i)
      vectAt(f, i) = std::log(ncmax(cfg.tail_floor, vectAt(f, i)));
    dlog = sqrtCurvatureDensity(f);
  }

  //Combine the densities. Weights are shares of the total integral. Terms which
  //are absent (functions with no curvature) are replaced by the constant term.
  double w_uniform = cfg.equidistant_fraction;
  double w_lin = 0.5 * (1.0 - cfg.equidistant_fraction);
  double w_log = w_lin;
  if (dlin.empty()) {
    w_uniform += w_lin;
    w_lin = 0.0;
  }
  if (dlog.empty()) {
    w_uniform += w_log;
    w_log = 0.0;
  }

  //Cumulative integral at each input point:
  VectD cum(n);
  cum.front() = 0.0;
  StableSumKahan csum;
  for (std::size_t i = 0; i + 1 < n; ++i) {
    double d = w_uniform * invL;
    if (w_lin > 0.0)
      d += w_lin * vectAt(dlin, i);
    if (w_log > 0.0)
      d += w_log * vectAt(dlog, i);
    csum.add(d * vectAt(width, i));
    vectAt(cum, i + 1) = csum.sum();
  }
  const double mtot = cum.back();
  nc_assert(mtot > 0.0);

  //Place points at equal steps of the cumulative integral, always using the
  //closest input point, but ensure that we get exactly targetN distinct points
  //in increasing order (so also leave room for the points still to come):
  std::vector<std::size_t> sel;
  sel.reserve(targetN);
  sel.push_back(0);
  for (std::size_t k = 1; k + 1 < targetN; ++k) {
    const double q = mtot * (static_cast<double>(k)
                             / static_cast<double>(targetN - 1));
    //first input point with cum >= q:
    std::size_t hi = static_cast<std::size_t>
      (std::lower_bound(cum.begin(), cum.end(), q) - cum.begin());
    std::size_t idx = hi;
    if (hi > 0 && (q - vectAt(cum, hi - 1)) <= (vectAt(cum, hi) - q))
      idx = hi - 1;
    const std::size_t lo_allowed = sel.back() + 1;
    const std::size_t hi_allowed = n - targetN + k;
    idx = ncmax(lo_allowed, ncmin(hi_allowed, idx));
    sel.push_back(idx);
  }
  sel.push_back(n - 1);
  nc_assert(sel.size() == targetN);

  VectD nx, ny;
  nx.reserve(targetN);
  ny.reserve(targetN);
  for (auto i : sel) {
    nx.push_back(x[i]);
    ny.push_back(y[i]);
  }
  return { std::move(nx), std::move(ny) };
}

NC::VectD::const_iterator NC::findClosestValInSortedVector(const VectD& v, double value)
{
  nc_assert(!v.empty());
  nc_assert(!ncisnan(value));
  auto it = std::lower_bound(v.begin(),v.end(),value);
  //it is the first element with *it >= value.
  if (it == v.begin())
    return it;
  if (it == v.end())
    return std::prev(v.end());
  //either it or std::prev(it), depending on which is closer:
  return ncabs(*it-value) < ncabs(*std::prev(it)-value) ? it : std::prev(it);
}

double NC::integrate01_kpowx( double k, const Optional<double>& opt_lnk )
{
  nc_assert(k>0.0&&std::isfinite(k));
  const double u = k - 1.0;
  if ( ncabs(u) < 0.1 ) {
    //direct Taylor, not even needing log:
    // For the record we simply got the Taylor coefficients with sagemath:
    // > sage: u,x=var('u,x');f=( ((1+u)**x).integrate(x,0,1) )
    //
    // First investigating number of orders needed for u=+-0.1 with command
    // (with a few orders added for safety):
    //
    // > sage: ((f-f.taylor(u,0,13))(u=1/10)).n()
    // > sage: ((f-f.taylor(u,0,13))(u=-1/10)).n()
    // Then generate the coefficients with:
    // > sage: print( '\n'.join(('constexpr double c%i = %s;'%(c[1],str(c[0])))
    // > ....:   .replace('/','./') for c in (f.taylor(u,0,13)).coefficients()))
    //
    // In this case, one can also simply evaluate the 14th term at u=0.1, giving
    // 4.21e-17, so the taylor expansion is indeed giving full precision.

    constexpr double c1 = 1.0/2.0;
    constexpr double c2 = -1.0/12.0;
    constexpr double c3 = 1.0/24.0;
    constexpr double c4 = -19.0/720.0;
    constexpr double c5 = 3.0/160.0;
    constexpr double c6 = -863.0/60480.0;
    constexpr double c7 = 275.0/24192.0;
    constexpr double c8 = -33953.0/3628800.0;
    constexpr double c9 = 8183.0/1036800.0;
    constexpr double c10 = -3250433.0/479001600.0;
    constexpr double c11 = 4671.0/788480.0;
    constexpr double c12 = -13695779093.0/2615348736000.0;
    constexpr double c13 = 2224234463.0/475517952000.0;
    return ( 1.0+u*(c1+u*(c2+u*(c3+u*(c4+u*(c5+u*(c6+u*(c7
             +u*(c8+u*(c9+u*(c10+u*(c11+u*(c12+u*c13)))))))))))) );
  }

  //Full formula, needs lnk=log(k):
  double lnk;
  if ( opt_lnk.has_value() ) {
    lnk = opt_lnk.value();
    nc_assert(std::isfinite(lnk));
    nc_assert(floateq(std::log(k),lnk));
  } else {
    lnk = std::log(k);
    nc_assert(std::isfinite(lnk));
  }
  return ( k - 1.0 ) / lnk;
}

double NC::integrate01_xkpowx( double k, const Optional<double>& opt_lnk )
{
  nc_assert(k>0.0&&std::isfinite(k));
  const double u = k - 1.0;
  if ( ncabs(u) < 0.22 ) {
    //direct Taylor, not even needing log:

    // For the record we simply got the Taylor coefficients with sagemath:
    // > sage: u,x=var('u,x');f=( (x*(1+u)**x).integrate(x,0,1) )
    // > sage: print( '\n'.join(('constexpr double c%i = %s;'%(c[1],str(c[0])))
    // > ....:   .replace('/','./') for c in (f.taylor(u,0,18)).coefficients()))
    //
    // The final order and thresholds were determined manually after checking
    // with mpmath, since the analytical (large u) formula below also has
    // numerical issues. To truly improve on this, one would likely need several
    // taylor domains - but for now the worst error seen by the entire function
    // was 2e-15 which is OK.

    constexpr double c0 = 1.0/2.0;
    constexpr double c1 = 1.0/3.0;
    constexpr double c2 = -1.0/24.0;
    constexpr double c3 = 7.0/360.0;
    constexpr double c4 = -17.0/1440.0;
    constexpr double c5 = 41.0/5040.0;
    constexpr double c6 = -731.0/120960.0;
    constexpr double c7 = 8563.0/1814400.0;
    constexpr double c8 = -27719.0/7257600.0;
    constexpr double c9 = 190073.0/59875200.0;
    constexpr double c10 = -516149.0/191600640.0;
    constexpr double c11 = 1013143139.0/435891456000.0;
    constexpr double c12 = -1519024289.0/747242496000.0;
    constexpr double c13 = 14108351869.0/7846046208000.0;
    constexpr double c14 = -14399405173.0/8966909952000.0;
    constexpr double c15 = 23142912688967.0/16005934264320000.0;
    constexpr double c16 = -83945247395407.0/64023737057280000.0;
    constexpr double c17 = 84894728616107.0/70959641905152000.0;
    constexpr double c18 = -3204549982389941.0/2919482409811968000.0;

    return (c0+u*(c1+u*(c2+u*(c3+u*(c4+u*(c5+u*(c6+u*(c7+u*(c8+u*(c9+u*(c10
   +u*(c11+u*(c12+u*(c13+u*(c14+u*(c15+u*(c16+u*(c17+u*c18))))))))))))))))));
  }
  //Full formula, needs lnk=log(k):
  double lnk;
  if ( opt_lnk.has_value() ) {
    lnk = opt_lnk.value();
    nc_assert(std::isfinite(lnk));
    nc_assert(floateq(std::log(k),lnk));
  } else {
    lnk = std::log(k);
    nc_assert(std::isfinite(lnk));
  }
  StableSum ss;
  ss.add(lnk);
  ss.add(-1.0);
  ss.mult(k);
  ss.add(1.0);
  return ss.sum() / ncsquare(lnk);
}

static_assert( NC::ncconstexpr_ispow2( 1 ), "" );
static_assert( NC::ncconstexpr_ispow2( 2 ), "" );
static_assert( NC::ncconstexpr_ispow2( 4 ), "" );
static_assert( NC::ncconstexpr_ispow2( 64 ), "" );
static_assert( NC::ncconstexpr_ispow2( 1024 ), "" );
static_assert( !NC::ncconstexpr_ispow2( 0 ), "" );
static_assert( !NC::ncconstexpr_ispow2( 3 ), "" );
static_assert( !NC::ncconstexpr_ispow2( 5 ), "" );
static_assert( !NC::ncconstexpr_ispow2( 48 ), "" );
static_assert( NC::ncconstexpr_roundupnextpow2(1) == 1, "" );
static_assert( NC::ncconstexpr_roundupnextpow2(2) == 2, "" );
static_assert( NC::ncconstexpr_roundupnextpow2(3) == 4, "" );
static_assert( NC::ncconstexpr_roundupnextpow2(4) == 4, "" );
static_assert( NC::ncconstexpr_roundupnextpow2(17) == 32, "" );
static_assert( NC::ncconstexpr_roundupnextpow2(24) == 32, "" );
static_assert( std::numeric_limits<double>::is_iec559,
               "NCrystal requires IEEE 754 floating point numbers" );

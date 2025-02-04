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
#include "NCrystal/internal/utils/NCRotMatrix.hh"
#include "NCrystal/internal/utils/NCIter.hh"
#include <sstream>
#include <list>

namespace NC = NCrystal;

bool NC::nc_is_grid(NC::Span<const double> v)
{
  if ( v.empty() )
    return false;
  double last = v.front();
  if ( ncisnan(last) || ncisinf(last) )
    return false;
  for ( auto e : Span<const double>(std::next(v.begin()),v.end()) ) {
    if ( !(e>last) || ncisnan(e) || ncisinf(e) )
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
  //Like this for highest numerical precision:
  for (unsigned i = 0; i<num_minus_1;++i)
    v.push_back(start+i*interval);
  v.push_back( stop );
  return v;
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
  if (!(b>a)||!fa*fb<0.0)
    NCRYSTAL_THROW(CalcError,"root finding requires b>a and f(a)*f(b)<0.");
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

std::pair<NC::VectD,NC::VectD> NC::reducePtsInDistribution( const NC::VectD& x,
                                                            const NC::VectD& y,
                                                            std::size_t targetN )
{
  nc_assert_always(x.size()==y.size());
  nc_assert_always(x.size()>=targetN);
  nc_assert_always(targetN>=2);
  if ( targetN >= x.size()  )
    return { x, y };
  const double ymax = *std::max_element(y.begin(),y.end());
  nc_assert_always(ymax>0.0);
  const double inv_ymax = 1.0/ymax;

  //Put point-data (including caches of expensive std::log results) into
  //doubly linked list:
  struct PtData;
  typedef std::multimap<double,std::list<PtData>::iterator> ImpMap_t;
  struct PtData {
    double x,y,lny;
    ImpMap_t::iterator impMapIter;
    PtData(double x_, double y_, double lny_, ImpMap_t::iterator it_)
      :x(x_), y(y_), lny(lny_), impMapIter(it_) {}
  };

  std::list<PtData> pts;
  ImpMap_t importanceMap;

  for ( auto&& e: enumerate(x) ) {
    double xval(e.val), yval(y.at(e.idx));
    pts.emplace_back(xval, yval, std::log(std::max<double>(1e-20,yval*inv_ymax)), importanceMap.end() );
  }

  //Definition of importance score:
  auto importanceCalc = [](const decltype(pts)::iterator it)
                        {
                          const auto& pt1 = *it;
                          const auto& pt0 = *std::prev(it);//previous neighbour
                          const auto& pt2 = *std::next(it);//next neighbour
                          //Area of triangle with 3 points at the corners is
                          //the area the curve integral will change if we
                          //remove the point (x1,y1). Apart from a missing
                          //factor of 0.5, the change in area is thus:
                          double area_change = ncabs(pt0.x*(pt1.y-pt2.y)+pt1.x*(pt2.y-pt0.y)+pt2.x*(pt0.y-pt1.y));
                          //To preserve also small features in the tails, we
                          //find the equivalent change in the logy-curve, and
                          //combine the two for the final importance (one can
                          //change the power of each factor to change how much
                          //focus should be on the features in tails and how
                          //much focus should be on strong peaks):
                          double logarea_change = ncabs(pt0.x*(pt1.lny-pt2.lny)+pt1.x*(pt2.lny-pt0.lny)+pt2.x*(pt0.lny-pt1.lny));
                          return area_change * logarea_change * logarea_change;
                        };

  //We keep multimap of importance -> PtData-iterator. Once initialised, the
  //first element corresponds to the least important remaining point on the
  //curve. It can be removed while also updating the importance scores of its
  //neighbours - repeat until done.
  auto itLast = std::prev(pts.end());
  for (auto it = std::next(pts.begin()); it!=itLast; ++it)
    it->impMapIter = importanceMap.emplace(importanceCalc(it),it);

  auto updateImportance = [&importanceMap,&importanceCalc](decltype(pts)::iterator it)
                          {
                            if (it->impMapIter == importanceMap.end())
                              return;//first or last point, nothing to update
                            importanceMap.erase(it->impMapIter);
                            it->impMapIter = importanceMap.emplace(importanceCalc(it),it);
                          };
  while( pts.size() > targetN ) {
    //Remove least important point and update importance scores of neighbours:
    nc_assert(!importanceMap.empty());
    auto impMapIter = importanceMap.begin();
    auto itPt = impMapIter->second;
    nc_assert(impMapIter==itPt->impMapIter);
    //Remove point:
    auto itPtN1 = std::prev(itPt);
    auto itPtN2 = std::next(itPt);
    pts.erase(itPt);
    importanceMap.erase(impMapIter);
    //Update neighbours:
    nc_assert(itPtN2!=pts.end());
    updateImportance(itPtN1);
    updateImportance(itPtN2);
  }
  //Done, put in vectors and return:
  VectD newx, newy;
  newx.reserve(pts.size());
  newy.reserve(pts.size());
  for (auto& pt: pts) {
    newx.emplace_back(pt.x);
    newy.emplace_back(pt.y);
  }
  return { newx, newy };
}

NC::VectD::const_iterator NC::findClosestValInSortedVector(const NC::VectD& v, double value)
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
static_assert( std::numeric_limits<double>::is_iec559, "NCrystal requires IEEE 754 floating point numbers" );

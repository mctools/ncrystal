////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2017 NCrystal developers                                   //
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

#include "NCMath.hh"
#include "NCRotMatrix.hh"
#include <sstream>
#include <algorithm>

namespace NCrystal {
  void evalY(double rand, double y, double& f_of_y, double& fprime_of_y)
  {
    // With a=sqrt(kT/m), the maxwell velocity distribution is:
    //
    //   f'(v) = sqrt(2/pi)/ * (1/a) * (v/a)^2 * exp[ - 0.5*(v/a)^2]
    //
    // And the commulative density function:
    //
    //   f(v) = erf[(v/a)*(1/sqrt(2))] - sqrt(2/pi) * (v/a) * exp[- 0.5*(v/a)^2]
    //
    // Making a variable change, y=v/a, we get (df/dv) / (1/a) = df/dy, and thus:
    //
    //   f'(y) = sqrt(2/pi)/ * y^2 * exp[ - 0.5*y^2]
    //   f(y) = erf[y/sqrt(2)] - sqrt(2/pi) * y * exp[- 0.5*y^2]
    const double y2 = y*y;
    const double k = 0.79788456080286541 * exp( - 0.5 * y2 ) * y;// 0.79788... = sqrt(2/pi)
    fprime_of_y =  k * y;
    f_of_y = ncerf(0.70710678118654757*y) - k - rand;//0.707... = 1/sqrt(2)
  }
}

double NCrystal::genThermalY(double rand)
{
  if (rand<=0.0)
    return 0.0;
  //First, bracket and estimate result using lookup table:
  const double mwlookup_delta = 0.01;
  const double mwlookup_maxr (1.0-1e-13);
  const double mwlookup_maxres (7.97482162522363);
  static const double mwlookup_table[] = {
    0, 0.3388684138410028, 0.4299207124908707, 0.4950744476442889,
    0.5478607658860285, 0.5931663491376359, 0.6333824737081255, 0.6698796170409779,
    0.7035250967074413, 0.7349067919985254, 0.7644438332246415, 0.7924468910441906,
    0.8191534750352306, 0.8447498018866729, 0.8693849737267244, 0.8931805216442653,
    0.9162370345787935, 0.9386388886086633,  0.960457699785782, 0.9817548962985233,
    1.002583668853801,   1.02299047305457,  1.043016203109471,   1.06269712050892,
    1.08206559735863,  1.101150717679314,  1.119978768567496,   1.13857364502852,
    1.156957186481244,  1.175149458698043,  1.193168991817706,  1.211032982734667,
    1.228757468406866,  1.246357475282379,  1.263847149013001,   1.28123986782289,
    1.298548342275297,  1.315784703688883,   1.33296058306604,  1.350087182014055,
    1.367175337470917,   1.38423557982608,  1.401278187914843,  1.418313239178709,
    1.435350657229599,  1.452400256951501,  1.469471787788336,  1.486574975746351,
    1.503719564604272,  1.520915356800695,  1.538172254455052,  1.555500300975186,
    1.572909723710827,  1.590410978128382,  1.608014794008489,  1.625732224204914,
    1.643574696552397,  1.661554069573791,  1.679682692715337,  1.697973471935873,
    1.716439941594798,   1.73509634372912,  1.753957715987648,  1.773039989707631,
    1.792360099885219,  1.811936109117952,  1.831787348000512,  1.851934574954198,
    1.872400159092117,  1.893208290500268,  1.914385223215159,  1.935959558060747,
    1.957962571859416,  1.980428605961377,   2.00339552479314,   2.02690526064548,
    2.051004468132225,  2.075745312458949,  2.101186428702039,  2.127394098201242,
    2.15444370455286,  2.182421553645298,  2.211427173624066,   2.24157625625212,
    2.273004468923244,  2.305872467704384,  2.340372599464152,  2.376738032044664,
    2.41525545510403,  2.456283185097282,  2.500277710809408,  2.547833926268537,
    2.599748568553671,   2.65712516969768,  2.721558385025618,  2.795483482915113,
    2.882910146158972,  2.991201681411585,  3.136464460374558,  3.368214175143001 };

  double xlow,xhigh,rts,f,df;

  if (rand>=mwlookup_maxr) {
    xlow = mwlookup_maxres-1.0e-13;
    xhigh = 10*mwlookup_maxres;
    evalY(rand,xhigh,f,df);
    if (f<=0)
      return xhigh;//we prefer not to return infinity, even when called with rand=1
    rts = 0.5*(xlow+xhigh);
  } else {
    int ltidx = int(rand/mwlookup_delta);
    const int nltm1 = sizeof(mwlookup_table)/sizeof(double)-1;
    nc_assert(ltidx>=0);
    if (ltidx<nltm1) {
      xlow = mwlookup_table[ltidx];
      xhigh = mwlookup_table[ltidx+1];
      rts = xlow + (rand/mwlookup_delta-ltidx)*(xhigh-xlow);
    } else {
      xlow = mwlookup_table[nltm1];
      xhigh = mwlookup_maxres;
      rts = 0.5 * (xlow+xhigh);
    }
  }

  //Push limits out slightly, to guard against rounding errors:
  xlow -= 1.0e-9;
  if (xlow<0)
    xlow=0.0;
  xhigh += 1.0e-9;
  nc_assert(rts>xlow&&rts<xhigh);

  //Now, apply Newton-Raphson root-finding to solve F(y)==rand for y
  nc_assert(xhigh>xlow);
  double dxold(ncabs(xhigh-xlow));
  double dx(dxold), temp;
  evalY(rand,rts,f,df);
  for (int j=1;j<=100;++j) {
    if ((((rts-xhigh)*df-f)*((rts-xlow)*df-f) >= 0.0) || (ncabs(2.0*f) > ncabs(dxold*df))) {
      dxold=dx;
      dx=0.5*(xhigh-xlow);
      rts=xlow+dx;
      if (xlow == rts)
        return rts;
    } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts -= dx;
      if (temp == rts)
        return rts;
    }
    if (ncabs(dx) < 1.0e-10 )
      return rts;
    evalY(rand,rts,f,df);
    if (f < 0.0)
      xlow=rts;
    else
      xhigh=rts;
  }
  NCRYSTAL_THROW(CalcError,"genThermalY: Max iterations exceeded");
  return 1.0;
}

bool NCrystal::ncis_sorted(std::vector<double>::const_iterator itb, std::vector<double>::const_iterator ite)
{
#if __cplusplus >= 201103L
  return std::is_sorted(itb,ite);
#else
  nc_assert_always(itb <= ite);
  if ( itb==ite )
    return true;
  ++itb;
  for(; itb != ite; ++itb )
    if( *itb < *(itb-1) )
      return false;
  return true;
#endif
}

double NCrystal::trapz(const std::vector<double> &y,const std::vector<double> &x)
{
  nc_assert( x.size()==y.size() );
  double sum=0.;
  for(unsigned i=1;i<x.size();i++) {
    sum += (x[i]-x[i-1])*(y[i]+y[i-1])/2.;
  }
  return sum;
}

void NCrystal::flip(const std::vector<double> & arr, std::vector<double> & arr_dst, bool opposite_sign)
{
  arr_dst.resize(arr.size());
  std::reverse_copy(arr.begin(),arr.end(),arr_dst.begin());

  if(opposite_sign) {
    std::vector<double>::iterator it(arr_dst.begin()), itE(arr_dst.end());
    for (;it!=itE;++it)
      *it = -(*it);
  }
}

std::vector<double> NCrystal::logspace(double start, double stop, unsigned num)
{
  std::vector<double> vec(num) ;
  double interval = (stop-start)/(num-1);
  for(std::vector<double>::iterator it=vec.begin();it!=vec.end();++it)  {
    *it =  pow(10,start);
    start += interval;
  }
  return vec;
}

std::vector<double> NCrystal::linspace(double start, double stop, unsigned num)
{
  nc_assert(num>1);
  std::vector<double> vec(num) ;
  double interval = (stop-start)/(num-1);
  for(std::vector<double>::iterator it=vec.begin();it!=vec.end();++it)  {
    *it =  start;
    start += interval;
  }
  return vec;
}

void NCrystal::concatenate(std::vector<double> & arr,const std::vector<double> & arr_back, unsigned skip_pos)
{
  arr.insert(arr.end(), arr_back.begin()+skip_pos, arr_back.end());
}

void NCrystal::gauleg_10_ord(const double x1, const double x2, std::vector<double>& x, std::vector<double>& w)
{
  static const double gau_legendre_x_10[] = { -9.73906528517171743431e-01,
                                              -8.65063366688984536346e-01,
                                              -6.79409568299024435589e-01,
                                              -4.33395394129247157888e-01,
                                              -1.48874338981631215706e-01,
                                              1.48874338981631215706e-01,
                                              4.33395394129247157888e-01,
                                              6.79409568299024435589e-01,
                                              8.65063366688984536346e-01,
                                              9.73906528517171743431e-01 };

  static const double gau_legendre_w_10[] = { 6.66713443086837109774e-02,
                                              1.49451349150580531377e-01,
                                              2.19086362515982180366e-01,
                                              2.69266719309991742204e-01,
                                              2.95524224714752925536e-01,
                                              2.95524224714752925536e-01,
                                              2.69266719309991742204e-01,
                                              2.19086362515982180366e-01,
                                              1.49451349150580531377e-01,
                                              6.66713443086837109774e-02 };

  nc_assert_always(x.size()==10 && w.size() == 10);

  for(unsigned i=0;i<10;i++)
    {
      w[i] = gau_legendre_w_10[i]*((x2-x1)/2.);
      x[i] = (gau_legendre_x_10[i]*(x2-x1) + (x2+x1)) /2.;
    }
}


double NCrystal::gaulobatto_grid(const std::vector<double>& old_beta, std::vector<double>& grid, std::vector<double>& totwgt)
{
  grid=old_beta;
  unsigned grid_size=grid.size();

  std::vector<double> betavec;
  betavec.reserve( grid.size()*10 );

  totwgt.resize(0);
  totwgt.reserve(grid.size()*10);

  std::vector<double> x(10), w(10);
  for(unsigned i=0;i<grid_size-1;i++)
    {
      gauloba_10_ord(grid[i], grid[i+1], x, w );
      if(betavec.empty())
        {
          betavec.insert(betavec.end(), x.begin(), x.end());
          totwgt.insert(totwgt.end(), w.begin(), w.end());
        }
      else
        {
          betavec.insert(betavec.end(), ++x.begin(), x.end());
          *(--totwgt.end()) += *w.begin(); //the last edge is the same as the first one
          totwgt.insert(totwgt.end(), ++w.begin(), w.end());
        }
    }

  if(true) //reflect
    {
      std::vector<double> neg_beta;
      flip(betavec,neg_beta,true);
      concatenate(neg_beta,betavec,1);
      std::swap(neg_beta,betavec);
      grid=betavec;

      totwgt[0]*=2;
      std::vector<double> neg_w;
      flip(totwgt,neg_w,false);
      concatenate(neg_w,totwgt,1);
      std::swap(neg_w,totwgt);
    }

  double bw =0.;
  for(std::vector<double>::const_iterator it = totwgt.begin(); it!=totwgt.end();++it)
    {
      bw+=*it;
    }
  return bw;
}

double NCrystal::gaulobatto_grid( unsigned grid_size, double lower_limit, double upper_limit,
                                  std::vector<double>& grid, std::vector<double>& totwgt)
{

  grid = logspace(log10(lower_limit),log10(upper_limit), grid_size-1);
  grid.insert(grid.begin(),0.);

  std::vector<double> betavec;
  betavec.reserve( grid.size()*10 );

  totwgt.resize(0);
  totwgt.reserve(grid.size()*10);

  std::vector<double> x(10), w(10);
  for(unsigned i=0;i<grid_size-1;i++)
    {
      gauloba_10_ord(grid[i], grid[i+1], x, w );
      if(betavec.empty())
        {
          betavec.insert(betavec.end(), x.begin(), x.end());
          totwgt.insert(totwgt.end(), w.begin(), w.end());
        }
      else
        {
          betavec.insert(betavec.end(), ++x.begin(), x.end());
          *(--totwgt.end()) += *w.begin(); //the last edge is the same as the first one
          totwgt.insert(totwgt.end(), ++w.begin(), w.end());
        }
    }

  if(true) //reflect
    {
      std::vector<double> neg_beta;
      flip(betavec,neg_beta,true);
      concatenate(neg_beta,betavec,1);
      std::swap(neg_beta,betavec);
      grid=betavec;

      totwgt[0]*=2;
      std::vector<double> neg_w;
      flip(totwgt,neg_w,false);
      concatenate(neg_w,totwgt,1);
      std::swap(neg_w,totwgt);
    }

  double bw =0.;
  for(std::vector<double>::const_iterator it = totwgt.begin(); it!=totwgt.end();++it)
    {
      bw += *it;
    }
  return bw;

}

void NCrystal::gauloba_10_ord(const double x1, const double x2, std::vector<double>& x, std::vector<double>& w)
{
  static const double gau_lobatto_x_10[] = { -1,
                                             -0.919533908166458814,
                                             -0.738773865105505075,
                                             -0.477924949810444496,
                                             -0.165278957666387025,
                                             0.165278957666387025,
                                             0.477924949810444496,
                                             0.738773865105505075,
                                             0.919533908166458814,
                                             1 };

  static const double gau_lobatto_w_10[] = { 0.0222222222222222222,
                                             0.13330599085107011,
                                             0.22488934206312645,
                                             0.292042683679683758,
                                             0.327539761183897457,
                                             0.327539761183897457,
                                             0.292042683679683758,
                                             0.22488934206312645,
                                             0.133305990851070111,
                                             0.0222222222222222222 };

  nc_assert_always(x.size()==10 && w.size() == 10);

  for(unsigned i=0;i<10;i++)
    {
      w[i] = gau_lobatto_w_10[i]*((x2-x1)/2.);
      x[i] = (gau_lobatto_x_10[i]*(x2-x1) + (x2+x1)) /2.;
    }
}

double NCrystal::simpsons_irregular(const std::vector<double> &y, const std::vector<double> &x)
{
  nc_assert_always(x.size()==y.size());

  double sum=0.;
  double onethird = 1./3;
  const int vecsize = x.size();

  std::vector<double> matele (9.,0);

  for(int i=0;i<vecsize-2;i+=2)
    {
      double x0 = x[i], x1=x[i+1], x2=x[i+2];
      double xx0 = x0*x0, xx1 = x1*x1, xx2 = x2*x2;
      matele[0]=xx0;
      matele[1]=x0;
      matele[2]=1.;

      matele[3]=xx1;
      matele[4]=x1;
      matele[5]=1.;

      matele[6]=xx2;
      matele[7]=x2;
      matele[8]=1.;

      RotMatrix mat (matele);
      Vector yvec (y[i],y[i+1],y[i+2]);
      mat.inv();
      Vector ceo = mat*yvec;
      sum += ceo.x()*onethird*(xx2*x2-xx0*x0) +  ceo.y()*.5*(xx2-xx0) +  ceo.z()*(x2-x0);
    }

  if(!(vecsize%2))
    sum += (y[vecsize-1]+y[vecsize-2])*.5*(x[vecsize-1]-x[vecsize-2]);

  return sum;
}

bool NCrystal::isPrime(unsigned n) {
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


void NCrystal::sincos_mpi2pi2(double A, double& cosA, double& sinA) {
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

void NCrystal::sincos_mpi8pi8(double A, double& cosA, double& sinA) {
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

double NCrystal::cos_mpipi(double A)
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
  return nccopysign(c,kPiHalf-Aabs);
}

double NCrystal::cos_mpi2pi2(double x)
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

double NCrystal::cos_mpi8pi8(double x)
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

double NCrystal::sin_mpipi(double A)
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
  return nccopysign(s,A);
}

double NCrystal::sin_mpi2pi2(double x)
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

double NCrystal::sin_mpi8pi8(double x)
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

double NCrystal::estimateDerivative(const Fct1D* f, double x, double h, unsigned order)
{
  nc_assert(f);
  nc_assert(h>0);
  nc_assert(order==4||order==6);
  if (order==6)
    return (  256.*f->eval(x+0.25*h)-256.*f->eval(x-0.25*h)-32.*f->eval(x+0.5*h)+32.*f->eval(x-0.5*h)
              -8.*f->eval(x+0.5*h)+8.*f->eval(x-0.5*h)+f->eval(x+h)-f->eval(x-h) ) / (90.*h);
  return (-f->eval(x+h)+8*f->eval(x+0.5*h)-8*f->eval(x-0.5*h)+f->eval(x-h))/(6.0*h);
}

double NCrystal::estimateSingleSidedDerivative(const Fct1D* f, double x, double h, unsigned order)
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

double NCrystal::findRoot(const Fct1D*f,double a, double b, double acc)
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
    double c = (a*fb-b*fa)/(fb-fa);
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

NCrystal::Fct1D::~Fct1D(){}

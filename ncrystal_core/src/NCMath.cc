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
#include "NCrystal/NCException.hh"
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

double NCrystal::trapz(const std::vector<double> &y,const std::vector<double> &x)
{
  if(x.size()!= y.size())
    NCRYSTAL_THROW(BadInput,"trapz: wrong input");

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

  if(opposite_sign)
    std::transform(arr_dst.begin(), arr_dst.end(), arr_dst.begin(),
        std::bind1st(std::multiplies<double>(),-1.));
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
  std::vector<double> vec(num) ;
  double interval = (stop-start)/(num-1);
  for(std::vector<double>::iterator it=vec.begin();it!=vec.end();++it)  {
    *it =  start;
    start += interval;
  }
  return vec;
}

static const double gau_legendre_x_10 []={-9.73906528517171743431e-01,
  -8.65063366688984536346e-01,
  -6.79409568299024435589e-01,
  -4.33395394129247157888e-01,
  -1.48874338981631215706e-01,
  1.48874338981631215706e-01,
  4.33395394129247157888e-01,
  6.79409568299024435589e-01,
  8.65063366688984536346e-01,
  9.73906528517171743431e-01};
static const double gau_legendre_w_10 []={6.66713443086837109774e-02,
  1.49451349150580531377e-01,
  2.19086362515982180366e-01,
  2.69266719309991742204e-01,
  2.95524224714752925536e-01,
  2.95524224714752925536e-01,
  2.69266719309991742204e-01,
  2.19086362515982180366e-01,
  1.49451349150580531377e-01,
  6.66713443086837109774e-02};

void NCrystal::gauleg_10_ord(const double x1, const double x2, std::vector<double>& x, std::vector<double>& w)
{
  if(x.size()!=10 || w.size()!=10)
    NCRYSTAL_THROW(BadInput,"gauleg_10_ord: wrong array size");


  for(unsigned i=0;i<10;i++)
  {
    w[i] = gau_legendre_w_10[i]*((x2-x1)/2.);
    x[i] = (gau_legendre_x_10[i]*(x2-x1) + (x2+x1)) /2.;
  }
}

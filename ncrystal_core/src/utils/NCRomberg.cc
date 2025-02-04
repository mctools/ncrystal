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

#include "NCrystal/internal/utils/NCRomberg.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCMsg.hh"

void NCrystal::Romberg::evalFuncMany(double* fvals, unsigned n, double offset, double delta) const
{
  double * it = &fvals[0];
  double nn = n;//only cast once
  for ( double i = 0; i < nn; ++i )
    *(it++) = evalFunc( offset + delta * i );
}

double NCrystal::Romberg::evalFuncManySum(unsigned n, double offset, double delta) const
{
  double sum = 0.0;
  double nn = n;//only cast once
  for ( double i = 0; i < nn; ++i )
    sum += evalFunc( offset + delta * i );
  return sum;
}

bool NCrystal::Romberg::accept(unsigned, double prev_estimate, double estimate,double,double) const
{
  return ncabs(estimate-prev_estimate)<1e-8;
}

void NCrystal::Romberg::convergenceError(double a, double b) const
{
  NCRYSTAL_RAWOUT("NCrystal ERROR: Romberg integration did not converge. Will"
                  " attempt to write function curve to ncrystal_romberg.txt"
                  " for potential debugging purposes.\n");
  writeFctToFile("ncrystal_romberg.txt", a, b,16384);//2^(maxlevel-2)+1,
                                                     //i.e. last amount of pts
                                                     //sampled at once (not
                                                     //exactly at those precise
                                                     //points though).
  NCRYSTAL_THROW(CalcError,"Romberg integration did not converge. Wrote"
                 " function curve to ncrystal_romberg.txt for potential"
                 " debugging purposes.");
}

double NCrystal::Romberg::integrate(double a, double b) const
{

  double h = (b-a);
  double fvals[17];//R(4,4) needs 17 equally spaced evaluations, we do them in one go:
  evalFuncMany(&fvals[0], 17, a, h*0.0625);

  //To reduce overhead, we unroll the calculations for R(n,k) up to R(5,5),
  //since they are anyway short enough to carry out before entering the main
  //loop:

  h *= 0.5;
  const double R00 = (fvals[0] + fvals[16])*h;
  const double R10 = h*fvals[8] + 0.5*R00;
  const double R11 = (4./3.)*R10 + (-1./3.)*R00;
  h *= 0.5;
  const double R20 = h*(fvals[4]+fvals[12]) + 0.5*R10;
  const double R21 = (4./3.) * R20 + (-1./3.)* R10;
  const double R22 = (16./15.) * R21 + (-1./15.) * R11;
  h *= 0.5;
  const double R30 = h*((fvals[2]+fvals[6])+(fvals[10]+fvals[14])) + 0.5*R20;
  const double R31 = (4./3.) * R30 + (-1./3.)* R20;
  const double R32 = (16./15.) * R31 + (-1./15.) * R21;
  const double R33 = (64./63.) * R32 + (-1./63.) * R22;
  h *= 0.5;
  const double R40 = h*(((fvals[1]+fvals[3])+(fvals[5]+fvals[7]))+((fvals[9]+fvals[11])+(fvals[13]+fvals[15]))) + 0.5*R30;
  const double R41 = (4./3.) * R40 + (-1./3.)* R30;
  const double R42 = (16./15.) * R41 + (-1./15.) * R31;
  const double R43 = (64./63.) * R42 + (-1./63.) * R32;
  const double R44 = (256./255.) * R43 + (-1./255.) * R33;

  if (accept(4,R33,R44,a,b))
    return R44;

  //R(4,4) was not enough, try R(5,5):
  const double c5 = evalFuncManySum(16, a+h*0.5, h);
  h *= 0.5;
  const double R50 = h*c5 + 0.5*R40;
  const double R51 = (4./3.) * R50 + (-1./3.)* R40;
  const double R52 = (16./15.) * R51 + (-1./15.) * R41;
  const double R53 = (64./63.) * R52 + (-1./63.) * R42;
  const double R54 = (256./255.) * R53 + (-1./255.) * R43;
  const double R55 = (1024./1023.) * R54 + (-1./1023.) * R44;

  if (accept(5,R44,R55,a,b))
    return R55;

  //Still not accepted. Use generic loop for R(6,6) or higher.

  //Set up cache arrays to keep row data of current and previous rows:
  const unsigned maxlevel = 16;
  double cache1[maxlevel], cache2[maxlevel];
  double *row_prev = &cache1[0], *row = &cache2[0];

  row_prev[0] = R50;
  row_prev[1] = R51;
  row_prev[2] = R52;
  row_prev[3] = R53;
  row_prev[4] = R54;
  row_prev[5] = R55;

  unsigned nj = 16;
  for(unsigned i = 6; i < maxlevel; ++i){
    double hh = h;
    h *= 0.5;
    nj *= 2;
    double c = evalFuncManySum(nj, a+h, hh);

    row[0] = h*c + 0.5*row_prev[0]; //R(i,0)

    double n_k = 1.;
    for(unsigned j = 0; j < i; ++j) {
      n_k *= 4.0;
      row[j+1] = ( n_k * row[j] - row_prev[j] ) / (n_k-1.0); //extrapolate value for R(i,j)
    }

    if (accept(i,row_prev[i-1],row[i],a,b))
      return row[i];

    std::swap(row_prev,row);
  }

  //Did not converge:
  convergenceError(a,b);

  return row_prev[maxlevel-1];//convergenceError() did not throw or otherwise die, so return best estimate.
}

#include "NCrystal/internal/utils/NCFileUtils.hh"
#include <fstream>
#include <iomanip>
void NCrystal::Romberg::writeFctToFile(const std::string& filename, double a, double b, unsigned n) const
{
  nc_assert_always(b>a);
  if (file_exists(filename)) {
    NCRYSTAL_WARN("Aborting writing of "<<filename<<" since it already exists");
    return;
  }
  std::ofstream ofs (filename.c_str(), std::ofstream::out);
  ofs << std::setprecision(20);
  ofs << "#ncrystal_xycurve\n";
  ofs << "#colnames = evalFuncManySum(n=1)xN;evalFuncMany(n=N);reldiff\n";
  VectD y;
  y.resize(n);
  double delta = (b-a)/(n-1);
  evalFuncMany(&y[0], n, a, delta);
  for (unsigned i = 0; i<n; ++i) {
    double x = (i+1==n?b:a+i*delta);
    //evalFunc might throw exception if user implemented evalFuncManySum, so
    //call the latter with n=1 instead of evalFunc:
    double y0 = evalFuncManySum(1,x,1e-10/*delta will be unused*/);
    ofs << x<<" "<<y0<<" "<<y.at(i)<<" "<< ncabs(y.at(i)-y0)/(ncmax(1e-300,ncabs(y0)))<<"\n";
  }
  NCRYSTAL_MSG("Wrote "<<filename);
}

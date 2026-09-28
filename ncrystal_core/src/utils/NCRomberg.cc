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

#include "NCrystal/internal/utils/NCRomberg.hh"
#include "NCrystal/internal/utils/NCMath.hh"
#include "NCrystal/internal/utils/NCMsg.hh"
#include "NCRomberg_FMA.hh"

void NCrystal::Romberg::evalFuncMany(double* fvals, unsigned n, double offset, double delta) const
{
  double * it = &fvals[0];
  double nn = n;//only cast once
  for ( double i = 0; i < nn; ++i )
    *(it++) = evalFunc( std::fma( delta, i, offset ) );
}

double NCrystal::Romberg::evalFuncManySum(unsigned n, double offset, double delta) const
{
  double sum = 0.0;
  double nn = n;//only cast once
  for ( double i = 0; i < nn; ++i )
    sum += evalFunc( std::fma( delta, i, offset ) );
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

//fixme: need to unit test the following new functions + use in NCMath.hh

double NCrystal::Romberg::fixedOrderIntegration9pts( const double* fvals )
{
  StableSum sum;
  fixedOrderIntegration9pts( fvals, sum );
  return sum.sum();
}

double NCrystal::Romberg::fixedOrderIntegration17pts( const double* fvals )
{
  StableSum sum;
  fixedOrderIntegration17pts( fvals, sum );
  return sum.sum();
}

double NCrystal::Romberg::fixedOrderIntegration33pts( const double* fvals )
{
  StableSum sum;
  fixedOrderIntegration33pts( fvals, sum );
  return sum.sum();
}

void NCrystal::Romberg::fixedOrderIntegration5pts( const double* fvals,
                                                   StableSum& tgt )
{
  //coeffs created via sb_ncdev_romberg_print_coeffs script
  constexpr double coeffs5[5] = {
    7./90., 16./45., 2./15., 16./45., 7./90.
  };
  for ( int i = 0; i < 5; ++i )
    tgt.add( coeffs5[i]*fvals[i] );
}

void NCrystal::Romberg::fixedOrderIntegration9pts( const double* fvals,
                                                   StableSum& tgt )
{
  //coeffs created via sb_ncdev_romberg_print_coeffs script
  constexpr double coeffs9[9] = {
    31./810., 512./2835., 176./2835., 512./2835., 218./2835., 512./2835.,
    176./2835., 512./2835., 31./810.
  };
  for ( int i = 0; i < 9; ++i )
    tgt.add( coeffs9[i]*fvals[i] );
}

void NCrystal::Romberg::fixedOrderIntegration17pts( const double* fvals,
                                                    StableSum& tgt )
{
  //coeffs created via sb_ncdev_romberg_print_coeffs script
  constexpr double coeffs17[17] = {
    3937./206550., 65536./722925., 22016./722925., 65536./722925.,
    27728./722925., 65536./722925., 22016./722925., 65536./722925.,
    3062./80325., 65536./722925., 22016./722925., 65536./722925.,
    27728./722925., 65536./722925., 22016./722925., 65536./722925.,
    3937./206550.
  };
  for ( int i = 0; i < 17; ++i )
    tgt.add( coeffs17[i]*fvals[i] );
}

void NCrystal::Romberg::fixedOrderIntegration33pts( const double* fvals,
                                                    StableSum& tgt )
{
  //coeffs created via sb_ncdev_romberg_print_coeffs script
  constexpr double coeffs33[33] = {
    64897./6816150., 33554432./739552275., 1245184./82172475.,
    33554432./739552275., 404992./21130065., 33554432./739552275.,
    1245184./82172475., 33554432./739552275., 14081968./739552275.,
    33554432./739552275., 1245184./82172475., 33554432./739552275.,
    404992./21130065., 33554432./739552275., 1245184./82172475.,
    33554432./739552275., 563306./29582091., 33554432./739552275.,
    1245184./82172475., 33554432./739552275., 404992./21130065.,
    33554432./739552275., 1245184./82172475., 33554432./739552275.,
    14081968./739552275., 33554432./739552275., 1245184./82172475.,
    33554432./739552275., 404992./21130065., 33554432./739552275.,
    1245184./82172475., 33554432./739552275., 64897./6816150.
  };
  for ( int i = 0; i < 33; ++i )
    tgt.add( coeffs33[i]*fvals[i] );
}

void NCrystal::Romberg::fixedOrderIntegration65pts( const double* fvals,
                                                   StableSum& tgt )
{
  //coeffs created via sb_ncdev_romberg_print_coeffs script
  constexpr double coeffs65[65] = {
    18977737./3987447750., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    29018619904./3028466566125., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    9608565248./1009488855375., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    29018619904./3028466566125., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    9609061744./1009488855375., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    29018619904./3028466566125., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    9608565248./1009488855375., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    29018619904./3028466566125., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    355891142./37388476125., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    29018619904./3028466566125., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    9608565248./1009488855375., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    29018619904./3028466566125., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    9609061744./1009488855375., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    29018619904./3028466566125., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    9608565248./1009488855375., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    29018619904./3028466566125., 68719476736./3028466566125.,
    22917677056./3028466566125., 68719476736./3028466566125.,
    18977737./3987447750.
  };
  for ( int i = 0; i < 65; ++i )
    tgt.add( coeffs65[i]*fvals[i] );
}

void NCrystal::Romberg::fixedOrderIntegration129pts( const double* fvals,
                                                     StableSum& tgt )
{
  //coeffs created via sb_ncdev_romberg_print_coeffs script
  constexpr double coeffs129[129] = {
    1223989321./514380759750., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    236111080914944./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    1049437669888./220512745568115., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    236111080914944./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    78707817290384./16538455917608625., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    236111080914944./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    1049437669888./220512745568115., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    236111080914944./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    236123451882074./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    236111080914944./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    1049437669888./220512745568115., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    236111080914944./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    78707817290384./16538455917608625., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    236111080914944./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    1049437669888./220512745568115., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    236111080914944./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    237697616576512./49615367752825875., 562949953421312./49615367752825875.,
    187672890966016./49615367752825875., 562949953421312./49615367752825875.,
    1223989321./514380759750.
  };
  for ( int i = 0; i < 129; ++i )
    tgt.add( coeffs129[i]*fvals[i] );
}

//detail_romberg_integrate: see NCRomberg_FMA.hh (included above).

double NCrystal::Romberg::integrate(double a, double b) const
{
  return NCRYSTAL_APPLY_C_NAMESPACE(detail_romberg_integrate)( this, a, b );
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

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2023 NCrystal developers                                   //
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

#include "NCrystal/NCDefs.hh"
#include "NCrystal/internal/NCFastConvolve.hh"
#include <cstdlib>

namespace NC = NCrystal;

//Temporarily uncomment the following define to test with safer but slower code:
//#define NCRYSTAL_FASTCONVOLVE_EXTRASAFEMATH

NC::FastConvolve::FastConvolve() = default;
NC::FastConvolve::~FastConvolve() = default;


void NC::FastConvolve::initWTable( unsigned n_size )
{
  //round n_size up to next power of 2:
  unsigned n_size_pow2 = 1;
  while (n_size_pow2<n_size)
    n_size_pow2 *= 2;
  unsigned n=0;
  unsigned tmp = 1;
  while ( tmp < n_size_pow2 ) {
    n+=1;
    tmp *= 2;
  }
  nc_assert(tmp==n_size_pow2);

  m_w.clear();
  m_w.reserve(n_size_pow2);
  for( unsigned i = 0; i < n_size_pow2; ++i ) {
    NC::PairDD res = NC::FastConvolve::calcPhase(i, n);
    m_w.emplace_back(std::complex<double>(res.first,res.second));
  }
}

void NC::FastConvolve::fftconv( const NC::VectD& a1, const NC::VectD& a2, NC::VectD& y, double dt )
{
  const int minimum_out_size = a1.size() + a2.size() - 1;

  std::vector<std::complex<double> > b1(a1.begin(),a1.end());
  fftd(b1,FT_forward,minimum_out_size);
  std::vector<std::complex<double> > b2(a2.begin(),a2.end());
  fftd(b2,FT_forward,minimum_out_size);


  nc_assert(b1.size()==b2.size());
  std::vector<std::complex<double> >::iterator itb1(b1.begin()), itb1E(b1.end()), itb2(b2.begin());
  while (itb1!=itb1E)
    *itb1++ *= *itb2++;

  fftd(b1,FT_inverse,minimum_out_size);

  y.resize(minimum_out_size);
  const double k = dt/b1.size();
  nc_assert(b1.size()==b2.size());
  nc_assert(y.size()<=b1.size());
  VectD::iterator ity(y.begin()), ityE(y.end());
  itb1 = b1.begin();
  for(;ity!=ityE;++ity,++itb1) {
#ifdef NCRYSTAL_FASTCONVOLVE_EXTRASAFEMATH
    //use std::abs which calls std::hypot behind the scenes (expensive but can avoid overflows)
    *ity = k * std::abs(*itb1);
#else
    //naive and simple, avoids std::hypot
    double a(itb1->real());
    double b(itb1->imag());
    *ity = std::sqrt(a*a+b*b)*k;
#endif
  }
}

void NC::FastConvolve::fftd( std::vector<std::complex< double> > &data, FastConvolve::caltype ct,
                             unsigned minimum_output_size )
{
  double output_log_size_fp = std::ceil(std::log2(minimum_output_size));
  nc_assert_always(output_log_size_fp<32);
  const int output_log_size = output_log_size_fp;
  const int output_size = ( 1 << output_log_size );//this is now minimum_output_size rounded up to next power of 2

  const unsigned wTableSizeNeeded = std::max<unsigned>(output_size,data.size());
  if ( m_w.size() < wTableSizeNeeded )
    initWTable( wTableSizeNeeded );

  nc_assert_always( data.size() <= (std::size_t)output_size );
  if( data.size() != (size_t)output_size )
    data.resize(output_size,std::complex<double>());

  for(int j=1;j<output_size-1;++j)
    {
      int i=0;
      for(int k=1,tmp=j; k<output_size; i=(i<<1)|(tmp&1),k<<=1,tmp>>=1)
        {
        }
      if(j<i)
        swap(data[i],data[j]);
    }

  const int jump = m_w.size()/output_size;
#ifndef NCRYSTAL_FASTCONVOLVE_EXTRASAFEMATH
  const double convfact = (ct==FT_inverse?-1.0:1.0);
#endif

  for(int i=0;i<output_log_size;++i){
    int z=0;
    const int i1 = (1<<i);
    const int i1m1 = i1-1;
    const int i2 = 1<<(output_log_size-i-1);
    for(int j=0;j<output_size;++j){
      if((j/i1)%2){
        std::complex<double>& data_j = data[j];
        std::complex<double>& data_sympos = data[j-i1];
#ifdef NCRYSTAL_FASTCONVOLVE_EXTRASAFEMATH
        //std::complex<> multiplication is slow since it takes care of proper inf/nan/overflow
        data_j *= (ct==inverse?std::conj(m_w[z*jump]):m_w[z*jump]);
        //and the -=,+= operators seems to carry significant overhead for some reason:
        std::complex<double> temp = data_sympos;
        data_sympos += data_j;
        temp -= data_j;
        data_j = temp;
#else
        //naive and simple is faster:
        const std::complex<double>& w_zjump = m_w[z*jump];
        const double a(data_j.real()), b(data_j.imag()), c(w_zjump.real()), d(convfact * w_zjump.imag());
        const double jr(a*c-b*d);
        const double ji(a*d+b*c);
        const double sr(data_sympos.real());//no by-ref access to real/imag parts
        const double si(data_sympos.imag());
        data_j.real( sr - jr );
        data_j.imag( si - ji );
        data_sympos.real( sr + jr );
        data_sympos.imag( si + ji );
#endif
        z += i2;
      } else {
        z = 0;
        //will be the same result for the next i1-1 loops, so skip ahead:
        j += i1m1;
      }
    }
  }

}

NC::PairDD NC::FastConvolve::calcPhase(unsigned k, unsigned n)
{
  //Calculate exp(i*2*pi*k/n^2) where n must a nonzero number n=1,2,3,4,...
  //and k must be a number in 0...n-1.
  //
  //NB: Using PairDD for (real,imag) parts, rather than std::complex<double> => for
  //efficiency (we don't need the safety checks for obscure cases....
  //
  //Idea, use exp(a*b)=exp(a)*exp(b) and a small cache of
  //exp(i2pi/2^N), N=0,1,2,. when k!=1, one can combine as in these examples:
  //
  // 13/32 = 1/32 + 12/32 = 1/32 + 3/8 = 1/32 + 1/8 + 1/4
  //
  // 17/32 = 1/32 + 1/2
  //
  // 11/32 = 1/32 + 5/16 = 1/32 + 1/16 + 1/4
  //

  //Trivial case:
  if ( k == 0 )
    return { 1.0, 0.0 };

  //Eliminate common factors of 2 in fraction k/n^2:
  while ( k%2==0 ) {
    nc_assert( n>=1 );
    n -= 1;
    k /= 2;
  }

  if ( k == 1 ) {
    //Fundamental form.
    nc_assert(n>=1);//k>0, so it follows from k<=n-1 that n>1

    //Cache high-precision numbers (from Python's mpmath module) (for
    //efficiency, accuracy and reproducability):
    static std::array<double, 20> cosvals = { -1.0, // cos(2pi/2^1)
                                              0.0, // cos(2pi/2^2)
                                              0.707106781186547524401, // cos(2pi/2^3)
                                              0.923879532511286756128, // cos(2pi/2^4)
                                              0.980785280403230449126, // cos(2pi/2^5)
                                              0.995184726672196886245, // cos(2pi/2^6)
                                              0.998795456205172392715, // cos(2pi/2^7)
                                              0.999698818696204220116, // cos(2pi/2^8)
                                              0.999924701839144540922, // cos(2pi/2^9)
                                              0.999981175282601142657, // cos(2pi/2^10)
                                              0.999995293809576171512, // cos(2pi/2^11)
                                              0.999998823451701909929, // cos(2pi/2^12)
                                              0.99999970586288221916, // cos(2pi/2^13)
                                              0.999999926465717851145, // cos(2pi/2^14)
                                              0.999999981616429293808, // cos(2pi/2^15)
                                              0.999999995404107312891, // cos(2pi/2^16)
                                              0.999999998851026827563, // cos(2pi/2^17)
                                              0.999999999712756706849, // cos(2pi/2^18)
                                              0.99999999992818917671, // cos(2pi/2^19)
                                              0.999999999982047294177 }; // cos(2pi/2^20)

    static std::array<double, 20> sinvals = { 0.0, // sin(2pi/2^1)
                                              1.0, // sin(2pi/2^2)
                                              0.707106781186547524401, // sin(2pi/2^3)
                                              0.382683432365089771728, // sin(2pi/2^4)
                                              0.195090322016128267848, // sin(2pi/2^5)
                                              0.0980171403295606019942, // sin(2pi/2^6)
                                              0.049067674327418014255, // sin(2pi/2^7)
                                              0.0245412285229122880317, // sin(2pi/2^8)
                                              0.0122715382857199260794, // sin(2pi/2^9)
                                              0.00613588464915447535964, // sin(2pi/2^10)
                                              0.00306795676296597627015, // sin(2pi/2^11)
                                              0.0015339801862847656123, // sin(2pi/2^12)
                                              0.000766990318742704526939, // sin(2pi/2^13)
                                              0.000383495187571395589072, // sin(2pi/2^14)
                                              0.00019174759731070330744, // sin(2pi/2^15)
                                              0.0000958737990959773458705, // sin(2pi/2^16)
                                              0.000047936899603066884549, // sin(2pi/2^17)
                                              0.0000239684498084182187292, // sin(2pi/2^18)
                                              0.0000119842249050697064215, // sin(2pi/2^19)
                                              0.00000599211245264242784288 }; // sin(2pi/2^20)

    double cosval = ( n < cosvals.size() ? cosvals.at(n-1) : std::cos(k2Pi / (double(n)*n)) );
    double sinval = ( n < sinvals.size() ? sinvals.at(n-1) : std::sin(k2Pi / (double(n)*n)) );
    return { cosval, sinval };
  }

  //Non-fundamental form, must combine results from several fundamental forms
  //using multiplication of complex numbers.

  nc_assert(k%2==1);//must be odd at this point
  PairDD factor1 = calcPhase(1, n);
  PairDD factor2 = calcPhase(k-1, n);
  return { (factor1.first*factor2.first-factor1.second*factor2.second),
           (factor1.first*factor2.second+factor1.second*factor2.first) };

}


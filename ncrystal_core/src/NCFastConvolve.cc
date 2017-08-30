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

#include "NCrystal/NCException.hh"
#include "NCFastConvolve.hh"
#include <cstdlib>
#include <cmath>
#include <algorithm>

//Temporarily uncomment the following define to test with safer but slower code:
//#define NCRYSTAL_FASTCONVOLVE_EXTRASAFEMATH

NCrystal::FastConvolve::FastConvolve(unsigned n_size)
{
  if(!n_size)
    NCRYSTAL_THROW(BadInput,"FastConvolve::FastConvolve size of the FFT can not be zero.");

  //round n_size up to next power of 2:
  unsigned n_size_pow2 = 1;
  while (n_size_pow2<n_size)
    n_size_pow2 *= 2;

  m_w.resize(n_size_pow2);
  for( size_t i = 0; i < n_size_pow2; ++i )
    m_w[i] = std::exp( std::complex<double>(0.,i*(M_PI*2.0/n_size_pow2)));
}


NCrystal::FastConvolve::~FastConvolve()
{
}

void NCrystal::FastConvolve::fftconv(const std::vector<double>& a1, const std::vector<double>& a2, std::vector<double>& y, double dt) const
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
  std::vector<double>::iterator ity(y.begin()), ityE(y.end());
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

void NCrystal::FastConvolve::fftd(std::vector<std::complex< double> > &data, FastConvolve::caltype ct,
                                  unsigned minimum_output_size ) const
{
  double output_log_size_fp = std::ceil(log2(minimum_output_size));//log2 not in std:: namespace until C++11
  nc_assert_always(output_log_size_fp<32);
  const int output_log_size = output_log_size_fp;
  const int output_size = ( 1 << output_log_size );//this is now minimum_output_size rounded up to next power of 2

  if( (size_t)output_size > m_w.size() )
    NCRYSTAL_THROW(BadInput,"output size is larger than the pre-computed w table.");
  if( data.size() > m_w.size())
    NCRYSTAL_THROW(BadInput,"input size is larger than the pre-computed w table");

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

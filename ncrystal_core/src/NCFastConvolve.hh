#ifndef NCrystal_FastConvolve_hh
#define NCrystal_FastConvolve_hh

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

#include "NCrystal/NCDefs.hh"
#include <vector>
#include <complex>

namespace NCrystal {

  class FastConvolve {

    // Class which, using the Fast-Fourier-Transform algorithm (FFT) and its
    // inverse (IFFT), transforms the input arrays a1 and a2 to the complex
    // array given by ('*' here means element-wise multiplication):
    //
    //   b = IFFT ( FFT(a1) * FFT(a2) )
    //
    // And places |b| * constant in the output vector y (of length
    // a1.size()+a2.size()-1) where the constant is given by dt/N where N is
    // y.size() rounded up to the next power of 2.
    //
    // This is intended for convolving different parts of phonon spectra into a
    // single part, thus merging all terms of a phonon expansion into one final
    // spectrum.

  public:
    FastConvolve(unsigned n_size);
    ~FastConvolve();
    void fftconv( const std::vector<double>& a1,
                  const std::vector<double>& a2,
                  std::vector<double>& y,
                  double dt) const;

  private:

    enum  caltype {FT_forward,FT_inverse};
    std::vector< std::complex<double> > m_w;

    //The actual fast-fourier transform algorithm. It reads a vector that has
    //less than m_size (i.e. 2^N) elements, and the output vector of FFT
    //contains m_size elements

    void fftd( std::vector<std::complex<double> > &inout, caltype ct,
               unsigned minimum_output_size ) const;

    //auto-correlation function, may be useful for MD calculation results only
    //void fft_auto(std::vector<std::complex<double> > &inout, caltype ct=forward) const;

  };
}

#endif


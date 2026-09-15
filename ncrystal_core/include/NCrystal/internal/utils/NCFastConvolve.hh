#ifndef NCrystal_FastConvolve_hh
#define NCrystal_FastConvolve_hh

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

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  class FastConvolve : private MoveOnly {

    // Class which performs convolutions via a Fast-Fourier-Transform (FFT)
    // algorithm.
    //
    // If M = a1.size() + a2.size() - 1 and N is the smallest power of
    // two with N >= M. The inputs are zero-padded to length N and transformed
    // using an iterative, radix-2, decimation-in-time complex FFT:
    //
    //        y = constant*Real( IFFT ( FFT(a1) * FFT(a2) ) )
    //
    // The output vector y includes a constant for normalisation.
    //
    // For reference, the alg used is a radix-2, iterative, in-place
    // Cooley–Tukey FFT algorithm, which is used to compute a zero-padded linear
    // convolution. Special attention is given to avoid using trigonometric
    // functions platform irreproducibilities

  public:
    FastConvolve();
    ~FastConvolve();
    FastConvolve( FastConvolve&& ) noexcept;
    FastConvolve& operator=( FastConvolve&& ) noexcept;

    void convolve( const VectD& a1, const VectD& a2, VectD& y, double dt);

    //Legacy version actually use |y| instead if y.real() for output production:
    void convolveLegacy( const VectD& a1, const VectD& a2, VectD& y, double dt);

    //Internal function for calculating exp(i*2pi*k/2^n), exposed for unit
    //testing:
    static PairDD calcPhase(unsigned long k, unsigned long n);

  private:
    struct Impl;
    Pimpl<Impl> m_impl;
  };
}

#endif


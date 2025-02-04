#ifndef NCrystal_FastConvolve_hh
#define NCrystal_FastConvolve_hh

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

#include "NCrystal/core/NCDefs.hh"

namespace NCRYSTAL_NAMESPACE {

  class FastConvolve : private MoveOnly {

    // Class which, using the Fast-Fourier-Transform algorithm (FFT) and its
    // inverse (IFFT), transforms the input arrays a1 and a2 to the complex
    // array given by ('*' here means element-wise multiplication):
    //
    //   b = IFFT ( FFT(a1) * FFT(a2) )
    //
    // And places |b| * constant in the output vector y (of length
    // a1.size()+a2.size()-1) where the constant is given by dt/N where N is
    // y.size() rounded up to the next power of 2.

  public:
    FastConvolve();
    ~FastConvolve();
    FastConvolve( FastConvolve&& ) noexcept;
    FastConvolve& operator=( FastConvolve&& ) noexcept;

    void convolve( const VectD& a1,
                   const VectD& a2,
                   VectD& y,
                   double dt);

    //Internal function for calculating exp(i*2pi*k/2^n), exposed for unit
    //testing:
    static PairDD calcPhase(unsigned k, unsigned n);

  private:
    struct Impl;
    Pimpl<Impl> m_impl;
  };
}

#endif


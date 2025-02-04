#ifndef NCTests_FPE_hh
#define NCTests_FPE_hh

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

#include <stdexcept>

// Used to enable floating-point exceptions for tests written in C/C++.

namespace NCTests {

  struct FloatingPointException : public std::runtime_error {
    using std::runtime_error::runtime_error;//same constructors as base class
  };

  //A call to the following function will install a signal handler for SIGFPE
  //which will print diagnostics, and raise a FloatingPointException. Calling it
  //multiple times has no further effect.
  //
  //This is a no-op on certain platforms (apply/windows):
  void catch_fpe();
}

#endif

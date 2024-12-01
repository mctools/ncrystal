#ifndef NCrystal_TestModUtils_hh
#define NCrystal_TestModUtils_hh

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   //
//                                                                            //
//  Copyright 2015-2024 NCrystal developers                                   //
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

// Macros and inline functions for implementing compiled test modules (to be
// loaded with the Lib class from loadlib.py.

#if defined (_WIN32) || defined (__CYGWIN__) || defined (WIN32)
#  define NCTEST_API __declspec(dllexport)
#elif defined(__GNUC__) || defined(__clang__)
#  define NCTEST_API __attribute__ ((visibility ("default")))
#else
#  define NCTEST_API
#endif
#define NCTEST_CTYPES extern "C" NCTEST_API

#include <iostream>
#include <cstdlib>
#include <stdexcept>

namespace nctest {
  namespace {
    inline void printErrAndExit(const std::exception &e) noexcept(true) {
      const std::runtime_error* stdrte = dynamic_cast<const std::runtime_error*>(&e);
      if (stdrte)
        std::cout<<"NCTest ERROR (std::runtime_error): "<<stdrte->what()<<std::endl;
      else
        std::cout<<"NCTest ERROR (unknown)"<<std::endl;
      //std::exit(1); FIXME: What to do here? Perhaps revisit how we handle the
      //exceptions in NCrystal's own chooks.
    }
  }
}
#define NCCATCH catch (std::exception& e) { nctest::printErrAndExit(e); }

#define NCTEST_CTYPE_DICTIONARY NCTEST_CTYPES const char * nctest_ctypes_dictionary()

#ifndef NCTESTMODUTILS_NO_NCRYSTAL_INCLUDE
#  include "NCrystal/core/NCDefs.hh"
namespace NC = NCrystal;
#endif


#endif

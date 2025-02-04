#ifndef NCrystal_TestModUtils_hh
#define NCrystal_TestModUtils_hh

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
#include <sstream>

namespace nctest {
  namespace {
    inline std::string& detail_lastErrorBuffer()
    {
      static std::string buf;
      return buf;
    }
    inline void printErrAndExit(const std::exception &e) noexcept(true) {
      std::ostringstream ss;
      const std::runtime_error* stdrte = dynamic_cast<const std::runtime_error*>(&e);
      if (stdrte)
        ss << stdrte->what();
      else
        ss << "NCTest ERROR (unknown)";
      detail_lastErrorBuffer() = std::move(ss).str();
    }
  }
}
#define NCCATCH catch (std::exception& e) { nctest::printErrAndExit(e); }

#define NCTEST_CTYPE_DICTIONARY NCTEST_CTYPES const char * nctest_ctypes_dictionary()

NCTEST_CTYPES const char * nctestdetail_get_lasterror()
{
  return nctest::detail_lastErrorBuffer().c_str();
}

NCTEST_CTYPES int nctestdetail_has_lasterror()
{
  return nctest::detail_lastErrorBuffer().empty() ? 0 : 1;
}

NCTEST_CTYPES void nctestdetail_clear_lasterror()
{
  nctest::detail_lastErrorBuffer().clear();
}

#ifndef NCTESTMODUTILS_NO_NCRYSTAL_INCLUDE
#  include "NCrystal/core/NCDefs.hh"
namespace NC = NCrystal;
#endif


#endif

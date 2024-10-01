
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

//////////////////////////////////////////////////////
//FIXME: To common test header:
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
  inline void printErrAndExit(const std::exception &e) noexcept(true) {
    const std::runtime_error* stdrte = dynamic_cast<const std::runtime_error*>(&e);
    if (stdrte)
      std::cout<<"NCTest ERROR (std::runtime_error): "<<stdrte->what()<<std::endl;
    else
      std::cout<<"NCTest ERROR (unknown)"<<std::endl;
    //std::exit(1);
  }
}
#define NCCATCH catch (std::exception& e) { nctest::printErrAndExit(e); }
//////////////////////////////////////////////////////


NCTEST_CTYPES const char * nctest_ctypes_dictionary()
{
  return "void nctest_catch_fpe()";
}

//todo: try to enable on osx/windows?
#if defined (_WIN32) || defined (__CYGWIN__) || defined (WIN32) || defined(__APPLE__)
#  define NCTEST_SKIP_FPE
#endif

#ifdef NCTEST_SKIP_FPE
//do nothing in release build or on osx
NCTEST_CTYPES void nctest_catch_fpe(){}
#else

#  include <cassert>
#  include <cfenv>
#  include <cstring>
#  include <csignal>
#  include <cstdio>
#  include <execinfo.h>
#  include <atomic>

namespace {
  void nctest_custom_sigfpe_handler( int sig, siginfo_t *info, void* ) {

    if (sig!=SIGFPE)
      return;

    //protect against secondary SIGFPE's:
    static bool inHandler = false;//fixme atomic_bool!
    if ( inHandler )
      return;
    inHandler = true;

    //flush pending output:
    std::cout.flush();
    std::cerr.flush();

    //Print diagnostics (TODO: Optionally print the stack-trace):
    int fpe = info->si_code;
    if      (fpe==FPE_INTDIV) printf("ERROR - FPE detected: integer divide by zero\n");
    else if (fpe==FPE_INTOVF) printf("ERROR - FPE detected: integer overflow\n");
    else if (fpe==FPE_FLTDIV) printf("ERROR - FPE detected: floating point divide by zero\n");
    else if (fpe==FPE_FLTOVF) printf("ERROR - FPE detected: floating point overflow\n");
    else if (fpe==FPE_FLTUND) printf("ERROR - FPE detected: floating point underflow\n");
    else if (fpe==FPE_FLTRES) printf("ERROR - FPE detected: floating point inexact result\n");
    else if (fpe==FPE_FLTINV) printf("ERROR - FPE detected: floating point invalid operation\n");
    else if (fpe==FPE_FLTSUB) printf("ERROR - FPE detected: subscript out of range\n");
    else printf("ERROR - FPE detected: unknown fpe\n");

    //Produce backtrace:
    const unsigned traceoffset = 1;//!=0, to avoid printing trace from the FPE
                                   //catcher itself, but better too small than
                                   // missing something.
    const unsigned tracelimit = 30 + traceoffset;
    void *trace_addr[tracelimit];
    std::size_t ntrace = backtrace(trace_addr, tracelimit);
    char ** trace_symbs = backtrace_symbols(trace_addr,ntrace);
    for(std::size_t itrace=traceoffset;itrace<ntrace;++itrace) {
      if (!trace_symbs[0]||!trace_symbs[itrace])
        break;
      printf("Backtrace [%i]: %s\n",int(itrace-traceoffset),trace_symbs[itrace]);
    }
    std::cout.flush();
    throw std::runtime_error("SIGFPE detected - invalid mathematical operation");
  }
}

NCTEST_CTYPES void nctest_catch_fpe()
{
  try {
    //Only on first call:
    static std::atomic<bool> first(true);
    bool btrue(true);
    if ( !first.compare_exchange_strong(btrue,false) )
      return;
    //Do it:
    struct sigaction sigact, initial_sa;
    std::memset (&sigact, 0, sizeof(sigact));
    sigact.sa_sigaction = nctest_custom_sigfpe_handler;
    sigemptyset(&sigact.sa_mask);
    sigact.sa_flags = SA_SIGINFO;
    feenableexcept(FE_DIVBYZERO|FE_INVALID);
    //We do not enable FE_INEXACT and FE_UNDERFLOW or FE_OVERFLOW since they are
    //too easy to trigger (FE_OVERFLOW can be triggered by DBL_MAX*a for a>1).
    int ret = sigaction(SIGFPE, &sigact, &initial_sa);
    if ( ret!=0 ) {
      printf("NCTest::FPE ERROR: Could not install FPE handler.\n");
      std::exit(1);
    }
  } NCCATCH;
}

#endif

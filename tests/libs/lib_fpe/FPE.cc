#include "TestLib_fpe/FPE.hh"

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

#if defined(__APPLE__) || defined(_WIN32) || defined(WIN32)
void NCTests::catch_fpe(){}
#else

#include <cassert>
#include <cfenv>
#include <cstring>
#include <csignal>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <execinfo.h>

namespace NCTests {
  namespace {
    void custom_sigfpe_handler( int sig, siginfo_t *info, void* ) {

      if (sig!=SIGFPE)
        return;

      //protect against secondary SIGFPE's:
      static bool inHandler = false;
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
      const unsigned traceoffset = 1;//to avoid printing trace from the FPE catcher, but better too small than missing something
      const unsigned tracelimit = 30 + traceoffset;
      void *trace_addr[tracelimit];
      size_t ntrace = backtrace(trace_addr, tracelimit);
      char ** trace_symbs = backtrace_symbols (trace_addr,ntrace);
      for(size_t itrace=traceoffset;itrace<ntrace;++itrace) {
        if (!trace_symbs[0])
          break;
        printf("Backtrace [%i]: %s\n",int(itrace-traceoffset),trace_symbs[itrace]);
      }

      std::cout.flush();
      throw FloatingPointException("SIGFPE detected - invalid mathematical operation");
    }

    bool& catch_fpe_installed() {
      static bool value = false;
      return value;
    }
  }
}

void NCTests::catch_fpe() {
  if (catch_fpe_installed())
    return;
  catch_fpe_installed() = true;
  struct sigaction sigact, initial_sa;
  std::memset (&sigact, 0, sizeof(sigact));
  sigact.sa_sigaction = NCTests::custom_sigfpe_handler;
  sigemptyset(&sigact.sa_mask);
  sigact.sa_flags = SA_SIGINFO;
  feenableexcept(FE_DIVBYZERO|FE_INVALID);
  //We do not enable FE_INEXACT and FE_UNDERFLOW or FE_OVERFLOW since they are
  //too easy to trigger (FE_OVERFLOW can be triggered by DBL_MAX*a for a>1).
  int ret = sigaction(SIGFPE, &sigact, &initial_sa);
  if ( ret!=0 )
    printf("NCTests::FPE ERROR: Could not install FPE handler.\n");
}

#endif


#include "TestLib_fpe/FPE.hh"

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
  //  static long disable_fpe_count = 0;
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

//disable/reenable fpe. Note nested functions of python context managers
//applying these might result in call order like:
//disable_catch_fpe->disable_catch_fpe->reenable_catch_fpe->reenable_catch_fpe,
//which is why we need the disable_fpe_count variable to ensure FPEs are not
//inadvertently reenabled before the final call to reenable_catch_fpe.
//
//void NCTests::disable_catch_fpe()
//{
//  if (!catch_fpe_installed())
//    return;
//  ++disable_fpe_count;
//  if (disable_fpe_count==1)
//    fedisableexcept(FE_DIVBYZERO|FE_INVALID);
//}
//
//void NCTests::reenable_catch_fpe()
//{
//  if (!catch_fpe_installed())
//    return;
//  --disable_fpe_count;
//  if (disable_fpe_count==0)
//    feenableexcept(FE_DIVBYZERO|FE_INVALID);
//}

#endif


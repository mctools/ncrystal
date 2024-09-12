#ifndef NCTests_FPE_hh
#define NCTests_FPE_hh

#include <stdexcept>

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

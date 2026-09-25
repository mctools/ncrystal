
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2026 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

include_guard()

# Determines contents of C++ define NCRYSTAL_FMADISPATCH_ATTR and
# NCRYSTAL_FMADISPATCH_ENABLED.
#
# If supported and needed, it will usually be defined as GCC/Clang's
# __attribute__((target_clones("default,fma"))), allowing portable binaries in
# which can use hardware FMA instruction on CPUs which support it, via runtime
# dispatch (an ifunc, resolved once by the loader, not per call). This allows
# usage of std::fma for higher precision and more portable results, while still
# retaining high performance in selected hotspot functions
#
# We can not guess this only from compiler-version preprocessor macros
# (e.g. Apple's Clang and LLVM keep separate version numbering, or musl-libc
# systems which do not support the GNU ifunc mechanism.
#
# IMPORTANT NOTE: If cross-compiling, try_run can not detect the final runtime
# environment, so we conservatively disable the usage!!

function( ncrystal_probe_fmadispatch resvar_defs )
  set( "${resvar_defs}" "" PARENT_SCOPE )

  if ( NOT DEFINED CACHE{NCRYSTAL_FMADISPATCH_SUPPORTED} )
    set( supported OFF )
    if ( CMAKE_CROSSCOMPILING )
      message( WARNING "FMA dispatch probe: skipped (cross-compiling, disabling)" )
    elseif ( NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|AMD64|i[3-6]86|x86)$" )
      # The "fma" target string is x86-specific; nothing to probe elsewhere
      # (and no point: other architectures, e.g. aarch64, either always have
      # hardware FMA already, or this whole technique does not apply):
      message( STATUS "FMA dispatch probe: skipped (non-x86 CMAKE_SYSTEM_PROCESSOR=${CMAKE_SYSTEM_PROCESSOR}, disabling)" )
    else()
      set( testsrc "${CMAKE_CURRENT_BINARY_DIR}/ncrystal_fmadispatch_probe.cc" )
      #NB: probes a plain extern "C" free function -- confirmed (a real
      #basictest.yml CI failure on macOS/Intel) that this is not just a
      #simplification: Apple Clang's target_clones lowering on Mach-O
      #silently produces no linkable definition for a *namespaced or
      #member* (C++-mangled) target_clones target, while the identical
      #attribute on a plain extern "C" function links and runs fine there.
      #This is why any function needing NCRYSTAL_FMADISPATCH_ATTR that
      #cannot live in an anonymous namespace (i.e. needs external/member
      #visibility) must route the actual decorated definition through an
      #extern "C" + NCRYSTAL_APPLY_C_NAMESPACE-wrapped free function
      #instead of decorating itself directly -- see
      #docs/devel_fma_attribute.md and NC::stable_exp/stable_expm1/
      #NCrystal::Romberg::integrate for the established pattern (which
      #this probe exists to validate is safe to rely on):
      file(
        WRITE "${testsrc}"
[[
#include <cmath>
#include <cstdio>
#include <cstddef>

extern "C" __attribute__((target_clones("default,fma")))
void ncrystal_probe_cmul( double* re, double* im,
                          const double* cre, const double* cim,
                          std::size_t n )
{
  for ( std::size_t i = 0; i < n; ++i ) {
    double a = re[i], b = im[i], c = cre[i], d = cim[i];
    re[i] = std::fma( a, c, -(b*d) );
    im[i] = std::fma( a, d, b*c );
  }
}

int main()
{
  const std::size_t n = 64;
  double re[n], im[n], cre[n], cim[n], re_ref[n], im_ref[n];
  for ( std::size_t i = 0; i < n; ++i ) {
    re[i] = re_ref[i] = 1.0 + 0.01*static_cast<double>(i);
    im[i] = im_ref[i] = 0.3*std::sin(0.1*static_cast<double>(i));
    cre[i] = std::cos(0.2*static_cast<double>(i));
    cim[i] = std::sin(0.2*static_cast<double>(i));
  }
  ncrystal_probe_cmul( re, im, cre, cim, n );

  //Independent reference computation (no target_clones/fma involved), to
  //check the dispatched function actually ran and gave a plausible (if not
  //necessarily bit-identical, since it is not compared to a forced
  //std::fma-only reference here) result -- catching e.g. an empty/no-op
  //resolver rather than a genuine dispatcher:
  double maxdiff = 0.0;
  for ( std::size_t i = 0; i < n; ++i ) {
    double a = re_ref[i], b = im_ref[i], c = cre[i], d = cim[i];
    double diff = std::fabs(re[i]-(a*c-b*d)) + std::fabs(im[i]-(a*d+b*c));
    if ( diff > maxdiff ) maxdiff = diff;
  }

  if ( maxdiff > 1e-9 ) {
    std::printf("FAIL: maxdiff=%.6g\n", maxdiff);
    return 1;
  }
  std::printf("OK\n");
  return 0;
}
]]
      )
      message( STATUS "FMA dispatch probe: compiling and running a self-check..." )
      #Using only the long-supported try_run keywords here (COMPILE_OUTPUT_VARIABLE
      #/ RUN_OUTPUT_VARIABLE), not the newer (CMake>=3.25) CXX_STANDARD/SOURCES
      #form, to keep working down to this project's CMake 3.16 minimum:
      try_run(
        run_result compile_ok
        "${CMAKE_CURRENT_BINARY_DIR}/ncrystal_fmadispatch_probe_bld"
        "${testsrc}"
        COMPILE_OUTPUT_VARIABLE compile_output
        RUN_OUTPUT_VARIABLE run_output
      )
      if ( NOT compile_ok )
        message( STATUS "FMA dispatch probe: FAILED to compile (fine, just means the technique is unavailable here):\n${compile_output}" )
      elseif ( NOT "${run_result}" STREQUAL "0" )
        message( STATUS "FMA dispatch probe: compiled but did not run correctly (result=${run_result}, output=${run_output}). Disabling." )
      else()
        message( STATUS "FMA dispatch probe: OK (compiled and ran correctly)" )
        set( supported ON )
      endif()
    endif()
    set( NCRYSTAL_FMADISPATCH_SUPPORTED "${supported}" CACHE INTERNAL
      "Whether an actual compile+link+run test confirmed the target_clones-based FMA dispatch technique works on this machine/toolchain" )
  endif()

  if ( NCRYSTAL_FMADISPATCH_SUPPORTED )
    #Due to a confirmed Clang codegen bug where an ODR-used inline function is
    #only reachable via a target_clones function's throw path (e.g. nc_assert's
    #LogicError construction) can get an undefined reference. Since we want to
    #allow nc_assert's inside functions with target_clones attribute, we disable
    #target_clones in clang debug builds:
    set(
      "${resvar_defs}"
      "#if defined(__clang__) && !defined(NDEBUG)\n#define NCRYSTAL_FMADISPATCH_ATTR\n#else\n#define NCRYSTAL_FMADISPATCH_ATTR __attribute__((target_clones(\"default,fma\")))\n#endif\n#define NCRYSTAL_FMADISPATCH_ENABLED 1\n"
      PARENT_SCOPE
    )
  endif()
endfunction()

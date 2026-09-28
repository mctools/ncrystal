#ifndef NCrystal_FMADispatch_hh
#define NCrystal_FMADispatch_hh

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

//Helper macros for writing a NCRYSTAL_FMADISPATCH_ATTR function that also
//gets an MSVC/x86-only "Windows-fast" twin, dispatched to at runtime. Usable
//unconditionally on every platform: on GCC/Clang, or on MSVC/ARM64, all the
//twin-related macros below expand to nothing, so NCRYSTAL_FMADISPATCH_DECLARATOR
//alone reduces to the same plain "NCRYSTAL_FMADISPATCH_ATTR void name" you
//would write by hand for a function with no Windows twin at all.
//
//MSVC has no equivalent of GCC/Clang's NCRYSTAL_FMADISPATCH_ATTR
//(target_clones): a decorated function's std::fma() calls fall back there to
//a slow software emulation (confirmed via a real Windows CPU trace of
//app_perfvdos: 65% of all runtime). MSVC *can* lower std::fma to hardware
//vfmadd, but only under /arch:AVX2, which is not a safe *global* flag (raises
//the whole binary's baseline, risking SIGILL on older/EVC-capped CPUs --
//checked against both our PyPI cibuildwheel and conda-forge Windows builds,
//neither of which sets any /arch: override today). So, mirroring the
//GCC/Clang mechanism, we instead compile a second, /arch:AVX2 clone of just
//the hot function(s) in a dedicated translation unit, and dispatch to it at
//runtime only if the CPU actually supports FMA3, falling back to the
//ordinary NCRYSTAL_FMADISPATCH_ATTR body otherwise. x86-64 only (ARM64 has
//hardware FMA unconditionally -- see docs/devel_fma_attribute.md's "ARM /
//non-x86" section); CMake only ever defines NCRYSTAL_DISPATCH_TO_WIN_FMA for
//MSVC targeting x86_64.
//
//Mechanism (see NCFastConvolve.cc/_FMA.hh/_WINFMA.cc for a worked example):
// - A function's ordinary body and its Windows-fast twin are kept side by
//   side in one NC<Name>_FMA.hh: the macros below make each function's
//   parameter list and body appear exactly once as plain source, with only
//   the leading declarator (linkage/attribute/name) picked by the ambient
//   NCRYSTAL_WIN_FMA (see NCRYSTAL_FMADISPATCH_DECLARATOR).
// - NC<Name>.cc includes that header directly (NCRYSTAL_WIN_FMA is never
//   defined there), so NCRYSTAL_FMADISPATCH_DECLARATOR always expands to the
//   ordinary, NCRYSTAL_FMADISPATCH_ATTR-decorated name, and (only when
//   NCRYSTAL_DISPATCH_TO_WIN_FMA is defined, i.e. MSVC/x86)
//   NCRYSTAL_FMADISPATCH_WINFMA_FORWARD checks ncrystalHasWinFma() and
//   forwards to the twin.
// - NC<Name>_WINFMA.cc is a no-op unless NCRYSTAL_WIN_FMA is defined for it
//   specifically (true only for this one filename, only on MSVC/x86 -- see
//   ncrystal_core/CMakeLists.txt), in which case it #includes the same
//   header, where NCRYSTAL_FMADISPATCH_DECLARATOR now expands to the extern
//   "C" + NCRYSTAL_APPLY_C_NAMESPACE-wrapped twin name instead, compiled
//   under CMake's per-file /arch:AVX2. Being a whole, separately compiled
//   function (not per-element intrinsics), MSVC remains free to
//   auto-vectorise its loop into packed hardware FMA.
// - The twin needs external, cross-TU linkage to be callable from the
//   ordinary function's body, so (mirroring rule 5 of
//   docs/devel_fma_attribute.md) it is a private extern "C" +
//   NCRYSTAL_APPLY_C_NAMESPACE-wrapped free function, never a bare
//   extern "C" one.

#ifdef NCRYSTAL_DISPATCH_TO_WIN_FMA

//No <windows.h> needed: IsProcessorFeaturePresent has no PF_* flag for
//AVX/FMA (confirmed the hard way -- an earlier version of this function
//tried PF_FMA3_INSTRUCTIONS_AVAILABLE, which does not exist and failed a
//real MSVC CI build with "undeclared identifier"). Detecting FMA3 support
//correctly needs both a CPUID check (ECX bits 12/28/27 for FMA/AVX/OSXSAVE)
//and an XGETBV check of XCR0 (bits 1+2, whether the OS actually
//saves/restores the YMM register state) -- the standard, widely documented
//Intel/Microsoft pattern for detecting AVX-family support at runtime:
#include <intrin.h>
#include <immintrin.h>

namespace NCRYSTAL_NAMESPACE {
  namespace detail {
    //Cached (checked once): whether the running CPU and OS actually support
    //hardware FMA3, i.e. whether it is safe to call into a /arch:AVX2-built
    //NC*_WINFMA.cc twin:
    inline bool ncrystalHasWinFma()
    {
      static const bool avail = []() -> bool {
        int cpuInfo[4] = { 0, 0, 0, 0 };
        __cpuid( cpuInfo, 1 );
        const bool osxsave = ( cpuInfo[2] & (1<<27) ) != 0;
        const bool avx     = ( cpuInfo[2] & (1<<28) ) != 0;
        const bool fma     = ( cpuInfo[2] & (1<<12) ) != 0;
        if ( !osxsave || !avx || !fma )
          return false;
        const unsigned long long xcr0 = _xgetbv( 0 );
        return ( xcr0 & 0x6 ) == 0x6;
      }();
      return avail;
    }
  }
}

#endif

//NCRYSTAL_FMADISPATCH_DECLARATOR(rettype,name): the leading declarator for
//a dispatched function's definition (everything up to, but not including,
//the parameter list). Conditionally (re)defined per translation unit,
//since NCRYSTAL_WIN_FMA is a per-file CMake property, not a global one:
#ifdef NCRYSTAL_WIN_FMA
#  define NCRYSTAL_FMADISPATCH_DECLARATOR(rettype,name) \
     extern "C" rettype NCRYSTAL_APPLY_C_NAMESPACE(detail_##name##_winfma)
#else
#  define NCRYSTAL_FMADISPATCH_DECLARATOR(rettype,name) \
     NCRYSTAL_FMADISPATCH_ATTR rettype name
#endif

//NCRYSTAL_FMADISPATCH_DECLARATOR_C(rettype,cname,name): variant of
//NCRYSTAL_FMADISPATCH_DECLARATOR for a function that must be extern "C" +
//NCRYSTAL_APPLY_C_NAMESPACE-wrapped even in its ordinary (non-winfma)
//form, e.g. one already using that wrapping for docs/devel_fma_attribute.md
//rule 5 (Apple/Mach-O) reasons independent of Windows -- cname is that
//function's own (already-established) unmangled name, while name is only
//used to derive the twin's name (detail_<name>_winfma), same as elsewhere:
#ifdef NCRYSTAL_WIN_FMA
#  define NCRYSTAL_FMADISPATCH_DECLARATOR_C(rettype,cname,name) \
     extern "C" rettype NCRYSTAL_APPLY_C_NAMESPACE(detail_##name##_winfma)
#else
#  define NCRYSTAL_FMADISPATCH_DECLARATOR_C(rettype,cname,name) \
     extern "C" NCRYSTAL_FMADISPATCH_ATTR rettype NCRYSTAL_APPLY_C_NAMESPACE(cname)
#endif

//NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(rettype,name,params): forward-declares
//the /arch:AVX2 twin (params must include its own parentheses, e.g.
//"(double* p, std::size_t n)"). No-op (and must NOT be followed by a
//semicolon) except in the one translation unit that both can and needs to
//call the twin:
#if defined(NCRYSTAL_DISPATCH_TO_WIN_FMA) && !defined(NCRYSTAL_WIN_FMA)
#  define NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(rettype,name,params) \
     extern "C" rettype NCRYSTAL_APPLY_C_NAMESPACE(detail_##name##_winfma) params;
#else
#  define NCRYSTAL_FMADISPATCH_WINFMA_DECLARE(rettype,name,params)
#endif

//NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(name,args): at the top of the ordinary
//function's body, forwards to the twin at runtime if supported (args must
//include its own parentheses, e.g. "(data1,data2,n)"), via a bare "return
//<call>;" -- legal (and a no-op on the returned value) even when the
//function returns void. Always safe to call (a no-op where dispatch is not
//applicable), and always needs its own trailing semicolon, like an
//ordinary statement:
#if defined(NCRYSTAL_DISPATCH_TO_WIN_FMA) && !defined(NCRYSTAL_WIN_FMA)
#  define NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(name,args) \
     do { \
       if ( ::NCRYSTAL_NAMESPACE::detail::ncrystalHasWinFma() ) \
         return NCRYSTAL_APPLY_C_NAMESPACE(detail_##name##_winfma) args; \
     } while(0)
#else
#  define NCRYSTAL_FMADISPATCH_WINFMA_FORWARD(name,args) do {} while(0)
#endif

#endif

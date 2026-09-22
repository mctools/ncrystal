# NCRYSTAL_FMADISPATCH_ATTR

## Overview

`NCRYSTAL_FMADISPATCH_ATTR` (see `ncapi.h.in`) is a GCC/Clang attribute
letting a portably built binary (as for Python wheels/conda packages)
still run a hot loop's `std::fma()` calls at hardware FMA speed via runtime
CPU dispatch. Internal NCrystal use only. Quick rules:

1. Use on small hotspot helper functions in anon namespace in .cc.
2. Only decorate a function that is **entirely** audited: every
   floating-point expression is either an explicit `std::fma(...)` call, or
   provably safe if silently contracted.
3. Never remove existing `std::fma()` calls; this only speeds them up.
4. Never use on a virtual member function or in headers.
5. Never use on `extern "C"` functions.
6. Never use `ncrange(...)`/`ValRange` (the `for (auto i : ncrange(n))`
   idiom) inside the decorated function's body; use a plain hand-written
   `for (std::size_t i = 0; i < n; ++i)` loop instead.
7. `nc_assert` (direct or indirect) are OK, since we do not apply the fma
   attribute in debug builds.

Detailed reasoning given below for reference.

## Details

**Why this exists.** Reproducibility across platforms requires writing
`a*c-b*d` as explicit `std::fma(a,c,-(b*d))`: a plain expression may or may
not be silently fused into a hardware FMA instruction depending on compiler,
flags and target, giving a different (if still correctly rounded) answer on
different machines. But `std::fma()` compiled for a generic/baseline x86-64
target is a genuine, fairly expensive library call, even on hardware that
does support FMA, since the compiler cannot know that at build time.
`__attribute__((target_clones("default,fma")))` solves this: it compiles
both a baseline and an FMA-using clone of the function into the same binary,
and an indirect function (ifunc), resolved once by the loader, picks the
right one for the actual runtime CPU. Measured ~5-15x speedup for a
`NCFastConvolve.cc`-shaped hot loop; no cost where unsupported.

**The contraction hazard (rule 1).** A `target_clones("default,fma")` clone
is compiled as if with `-mfma`, which also changes that clone's default
floating-point contraction behaviour. This means it can silently fuse a
*plain* (non-explicit-`std::fma`) expression into hardware FMA too, purely
because that specific clone's target ISA happens to have one -- confirmed by
measurement: an unrelated plain expression in the same function differed
from the same code built without the attribute, on both GCC and Clang.
Decorating a function that mixes
audited `std::fma` calls with other, unaudited plain arithmetic is therefore
unsafe: it would reintroduce, self-inflicted within a single binary, exactly
the platform-dependent-rounding problem this attribute exists to avoid.

**Class methods.** Confirmed working, at full speed, for free functions and
for `static` and non-static (non-virtual) member functions. Virtual member
functions are rejected outright by GCC ("sorry, unimplemented: virtual
function multiversioning not supported") -- do not use it there.

**Never `extern "C"` (rule 4).** Not just unneeded -- actively harmful.
`extern "C"` on a function inside an anonymous namespace defeats that
namespace's local linkage in the object NCrystal actually ships (confirmed
with `nm -D`/`readelf --dyn-syms` on a real build, with this project's
ordinary, non-`-fvisibility=hidden` flags): the symbol leaks into the
shared library's *dynamic* symbol table as a `GLOBAL`/`WEAK` export, even
though the equivalent plain (mangled) C++ symbol correctly stays local.
This is true independent of `target_clones` -- a plain `extern "C"`
function with no attribute at all leaks the same way. An unnamespaced,
`WEAK` symbol is exactly what
`NCRYSTAL_NAMESPACE_PROTECTION`/`NCRYSTAL_C_NAMESPACE` (see `ncapi.h.in`)
exist to prevent: it lets two differently namespaced NCrystal builds
loaded into the same process (e.g. simplebuild's `NCrystalDev` and a
pip-installed `NCrystal`) silently bind to *each other's* definition of
the same-named function -- fine while the two builds happen to agree, a
silent correctness bug the moment they do not. Confirmed working, at full
speed, without `extern "C"`, on both a plain free function (ordinary
mangled C++ linkage) and on member functions (which cannot be given C
linkage at all -- see above). If a stable, unmangled name is ever
genuinely needed (e.g. calling the function from ctypes, or across a
shared-library boundary the way the stable C API in `ncrystal.cc` does),
that function is not a private, file-local kernel any more, and needs the
same namespacing guard `ncrystal.h`/`ncrystal.cc` already use for that
purpose -- a different, bigger undertaking than this document covers.

**Declaration vs. definition, and translation-unit scope (rule 5).** The
attribute does not technically need to appear identically on both a
declaration and its definition (confirmed for free functions and static
member functions, across TUs and even across a shared-library boundary: the
dispatch is resolved at the symbol level by the dynamic linker, independent
of what any given caller's compiler saw). Keeping the definition and its
(only) use together in one file is therefore not required for correctness --
it is recommended because it is the simplest way to guarantee the audit in
rule 1 cannot drift out of sync across a header/`.cc` split, and it avoids
quietly growing the library's exported symbol table.

**No `ncrange`/`ValRange` inside the decorated function (rule 6).**
Confirmed Clang bug (21.1.8, `-O3 -DNDEBUG`, minimal repro): `target_clones`
+ `for (auto i : ncrange(n))` can produce an undefined reference to
`ValRange`'s `inline constexpr` ctor instead of the usual weak definition.
Heuristic-dependent (doesn't always trigger, GCC unaffected) -- avoid
unconditionally rather than rely on non-reproduction.

**`nc_assert` inside the decorated function, and why it needed a global
fix rather than a per-function one (rule 7).** Confirmed Clang bug (21.1.8,
Debug/`-UNDEBUG`): the same "target_clones + ODR-used inline function ->
undefined reference instead of a weak definition" pattern as rule 6, but
triggered via `nc_assert`'s throw path (`LogicError`'s in-class-defined,
implicitly-inline constructor) rather than `ncrange`. First seen via
`ncclamp`/`nclerp`, both of which carry `nc_assert` internally, but the
underlying bug is general: *any* implicitly-inline call reachable only from
a throw path can trigger it, not just these two. Auditing every current and
future dispatched function's transitive call graph for this is exactly the
kind of thing that's easy to miss (as happened here) and easy to regress.
Fixed at the root instead: `NCRYSTAL_FMADISPATCH_ATTR` itself expands to
nothing when `defined(__clang__) && !defined(NDEBUG)` (in both
`ncrystal_fmadispatch.cmake` and simplebuild's `sbgen/main.py`, which have
independent copies of this logic). `NCRYSTAL_FMADISPATCH_ENABLED` is left
at `1` regardless -- it already only gates a *speedup* assertion, and that
assertion is itself only enforced in optimised (`NDEBUG`) builds, which
this condition never touches. Debug builds decorated with the (now empty)
attribute simply don't multiversion, so `nc_assert` and everything else
work exactly as in any ordinary function; Release/RelWithDebInfo/MinSizeRel
are unaffected for `nc_assert` specifically, since it compiles to nothing
there regardless of the bug, and GCC is unaffected on any config (bug is
Clang-specific). This fix is deliberately scoped to `!defined(NDEBUG)`
because that is exactly where `nc_assert` (unlike `nc_assert_always`,
which is not `NDEBUG`-gated) can be reached at all -- `nc_assert_always`
is confirmed (same minimal-repro method) to hit the identical undefined
reference in an optimised Clang build, completely unprotected by this
fix, which is why rule 7 keeps it as a real per-function audit item
rather than folding it into "handled automatically."

**ARM / non-x86.** Not a rule to apply -- this needs no attention during
everyday work, since the macro is already a safe no-op there and code must
always be written to give correct results on x86 too. Noted here only for
background: AArch64 has fused-multiply-add as a mandatory part of its
baseline ISA (unlike x86, where FMA3 is an optional extension needing
runtime detection), so `std::fma()` there already compiles to a hardware
instruction without needing this technique -- confirmed by measurement:
dispatched and never-dispatched reference timings were identical on
aarch64.

**How support is detected.** Not by guessing from compiler-version
preprocessor macros -- those are unreliable across the platforms NCrystal
targets (e.g. Apple's Clang has its own version numbering; musl-libc systems
lack the ifunc mechanism this relies on regardless of compiler version). See
`ncrystal_fmadispatch.cmake` for the actual compile+link+run probe used
instead.

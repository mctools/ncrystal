# Session summary: VDOS/SAB floating-point reproducibility, and a git history cleanup

This documents a Claude Code session on the `tk_volatile` branch, spent
chasing cross-platform floating-point reproducibility failures in the
VDOS/SAB pipeline, then rewriting the branch's commit history into a
coherent form. Written as a record of what was found and why, for anyone
picking this back up later.

## Background: the standing policy

NCrystal's CI matrix spans many compiler/platform/build-type
combinations (GCC/Clang, x86/ARM, Linux/macOS, Debug/Release). A
recurring class of CI failure was two platforms computing slightly
different results for the same VDOS/SAB calculation -- not a logic bug,
but **silent floating-point contraction**: a compiler fusing a plain
`a*b+c` into a single hardware fused-multiply-add instruction, which is
correctly rounded but not bit-identical to separate multiply-then-add.
Baseline x86-64 lacks hardware FMA and can never contract; AArch64 has it
as mandatory baseline ISA and contracts by default; x86 needs `-mfma` (or
runtime CPU dispatch) to do the same.

The project's standing rule (see `doc/devel_fma_attribute.md`): **never**
paper over this with `-ffp-contract=off` / `-DNCRYSTAL_DISABLE_FPCONTRACTION=ON`,
not even scoped to CI. That treats the symptom, not the cause, and leaves
the offending expression un-audited. The fix is always to find the
specific expression and rewrite it as an explicit `std::fma(...)` call,
which gives the same, single-rounded result regardless of whether the
compiler would have contracted it or not. This was tried once as a CI
workaround earlier in the branch's history and explicitly reverted; that
revert (and the workaround it undid) were dropped entirely in the later
history rewrite (see below) since together they carry no lasting
information.

A useful trick used throughout: **`-mfma`-forced local builds simulate ARM**.
Compiling normally on x86_64 never contracts (no hardware FMA in the
baseline target), so bugs that only manifest on AArch64 CI legs are
invisible locally -- unless you build with `CXXFLAGS=-mfma CFLAGS=-mfma`,
which makes GCC/Clang contract on x86 too, reproducing the same class of
bug without needing ARM hardware.

## What was fixed this session

Continuing earlier work in the same vein (`8f920e5a` FastConvolve,
`a407b5c3`/`24fe4cc4` VDOS density interpolation via `NC::nclerp`,
`6ebdc1ea` the SAB alpha/beta accumulate loop, `4f5d6fa3` the
`safe_x{cothx,3cothx,divsinhx}` Taylor series, `fc00c38c` the
`NCSABUtils.hh` Horner polynomials and a `NCVDOSEval.cc` bin-edge grid
computation), this session found and fixed two more:

1. **`NC::regulariseVDOSGrid`'s energy-grid interpolation**
   (`NCVDOSEval.cc`): the `new_emax` and per-point `eval` computations
   used plain `emin + binwidth * i` instead of `std::fma`.
2. **`stirlingsSeriesSum9thOrder`** (`NCVDOSToScatKnl.cc`): a 9-term
   Horner polynomial, used to estimate `n!` for high phonon orders, had
   no `std::fma` at all.

The second one mattered more than a typical last-digit fix: it can flip a
*discrete* downstream decision (how many phonon orders/how much of the
beta grid gets included) once the input curve or the material has enough
orders in play (heavy elements, high `vdoslux`/`knllux`). This was the
root cause of `tests/scripts/n2endf_bad.py`'s ThO2/Th beta-range
divergence, which had looked "gcc-10-only" but turned out to reproduce on
any platform/compiler with active FMA contraction (confirmed later to
also affect some macOS/Clang and ARM legs in a fresh CI run).

### Diagnostic tooling added: `tests/src/app_vdos2knldiag`

Built specifically to chase the `n2endf_bad` divergence without needing
gcc-10 (unavailable in the sandbox) or real ARM/macOS hardware: a verbose
C++ test app calling the exact same two functions the C API's
`ncrystal_raw_vdos2kernel` (and hence `NCrystal.vdos.extractKnl`, and
hence `ncmat2endf`) uses --
`NC::VDOS::createScatteringKernel()` + `NC::SABUtils::transformKernelToStdFormat()`
-- for ThO2's Th and O elements, dumping full alpha/beta/sab grids at
%.17g precision to a golden log.

This paid off immediately: forcing `-mfma` locally reproduced the exact
CI divergence in the tool's own output (beta range landing at ~11.99 vs
the reference's ~12.02), *before* either of the two fixes above existed.
After both fixes, the tool's output is bit-identical across
GCC/Clang x default-vs-`-mfma` builds (verified by building the same
source in the Release cache with plain flags and the Debug cache with
`-mfma` forced, and diffing the two ~1.9MB logs -- zero differences).

The golden log is deliberately huge (~1.9MB) for now, per explicit
instruction: added to the size-override table in
`devel/pypath/ncrystal_repo_tools/_check_misc.py` with a `#fixme` note;
trimming it down is left for later, once the tool has done its job on a
real CI run if further divergences ever show up.

### Verification methodology used throughout

For every fix: (1) confirm the *un-fixed* code diverges between a plain
build and an `-mfma`-forced build (via disassembly showing `vfmadd*`
instructions in supposedly-unaudited functions, and/or diffing computed
values); (2) apply the explicit-`std::fma` fix; (3) re-verify plain vs
`-mfma` (GCC and Clang, and Release-vs-Debug-cache to get genuinely
separate object files) are now bit-identical; (4) only then regenerate
any stale reference data via each test's own update mechanism, having
independently confirmed the new values aren't a regression.

Reference data regenerated this session: `tests/data/sabxs_A` (several
`vdoslux=2003/2004` cases -- the ones exercising enough phonon orders to
hit the Stirling-series bug), `tests/scripts/n2endf_cliex.log`.
`tests/scripts/n2endf_bad.log` needed no change: it already encoded the
"correct" (non-contracted) value, and the fixes simply stopped other
platforms from disagreeing with it.


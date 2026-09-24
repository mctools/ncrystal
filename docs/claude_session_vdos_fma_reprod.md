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

## Follow-up session: a residual sabxs failure not caused by FP contraction

A later session kept chasing a *different* residual: after the
`determineEMinDivKTminRange` EMax-anchoring fix (`69c55b0e`) and the
`sIntegralAtE` threshold taper (`2bb1714c`, itself real but confirmed by a
real CI re-run to not fix this residual), CI still showed exactly two
failing legs -- `ubuntu-22.04.gcc-10.python-3.9-Release` and
`macos-latest.clang.python-3.13-Release` -- both on the same
`py_rl_sabxs` case (`Li_from_Li2O.ncmat;vdoslux=2004;knllux=4;temp=10`),
with reldiffs (1.29e-05, 1.16e-04) far too large to be last-bit noise and
completely unmoved by any local `-mfma`/`-ffp-contract=off` toggling
(expected: this residual has nothing to do with FMA contraction).

Two more local repro attempts, both genuine negatives:
- **A real gcc-10 toolchain**, extracted from Ubuntu 22.04 `.deb`
  packages into an isolated prefix (no system install) and verified via
  `readelf -p .comment` on the actual `.o` files to rule out a stale-cache
  false positive. Byte-identical results vs. the local gcc-15 build.
- **glibc's own internal math dispatch**: since glibc 2.34, `log`/`exp`/`pow`
  live in `libc.so.6` and use runtime IFUNC dispatch to pick an
  FMA-optimised implementation based on detected CPU features --
  completely independent of how NCrystal itself is compiled or flagged.
  Toggled locally via `GLIBC_TUNABLES=glibc.cpu.hwcaps=-FMA,-FMA4` (further
  hwcaps such as AVX2/AVX512 made no additional difference on this CPU).
  This *did* produce a real, reproducible difference in the final sabxs
  value for the failing case, but only at the ~7e-8 relative level --
  2-3 orders of magnitude below the actual CI-observed failure, and
  unaffected by the eventual code fix below (see next section). Also,
  since macOS uses Apple's own libm rather than glibc, this diagnostic can
  never reach the macOS leg's own root cause even in principle -- it is a
  real, independent source of last-few-ULP noise, not *the* mechanism.

### Root cause traced to grid-point selection, not `e_touch` arithmetic

Using the glibc toggle plus targeted (fully-reverted-before-commit) debug
instrumentation, the divergence was traced through
`determineEMinDivKT`'s `f_of_e` evaluation, through every point of the
`m_eGrid`/`m_sIntegral` table (all ~449 points differing, 1.5e-5 to 4e-5
relative -- uniform, not tie-like), into `SABCellSurvey`'s constructor,
where a **different alpha/beta cell** (`b2` differing by ~0.4%, not
last-ULP) was found at the same iteration position between builds. This
first looked like it implicated `NCSABSurveyor.cc`'s `CellInfo::e_touch`
computation (`(alpha-beta)^2/(4*alpha)` at grid corners -- naive and,
unlike `getAlphaMinus`/`getBetaMinus`, with no available Taylor-expansion
hardening, since alpha and beta come from two unrelated grids rather than
one analytic formula with a removable singularity). An `e_touch`-cutoff
safety-margin fix was drafted for `NCSABProcessor.cc`'s `sIntegralAtE`
along these lines, but on reflection (and after the user pointed out the
flaw directly) this was the wrong layer: if two builds end up with
genuinely *different* alpha/beta grids, no downstream margin on a fixed,
already-different `survCells` array can reconcile them -- each build only
ever sees its own grid's cells. The fix was reverted before being tested
against CI.

The actual grid-construction call chain for next-gen SAB data is
`NCVDOSToScatKnl.cc` -> `VDOS::determineAlphaBetaGridFromGn`
(`NCVDOSKnlGrid.cc:155`) -> `NC::reducePtsByEquidistribution`
(`NCMath.cc:816`). The Gn expansion function values feeding this are
computed via `std::exp(-x + n*std::log(x) - minus_log_nfactorial)`
(`NCVDOSKnlGrid.cc:451`); `reducePtsByEquidistribution` then takes
`std::log` of those again to build a curvature-based point density,
integrates it into a cumulative distribution `cum[]` (Kahan-summed), and
picks which of the candidate grid points survive via, for each target
quantile `q`, a `std::lower_bound(cum.begin(), cum.end(), q)` search
followed by a "closest neighbour" tie-break (`NCMath.cc:956`, pre-fix:
`if (hi>0 && (q-cum[hi-1]) <= (cum[hi]-q)) idx = hi-1;`). A last-ULP
difference in `cum[]` (ultimately from glibc's exp/log dispatch, or any
other libm implementation difference) can flip which side of a near-tie
this raw comparison lands on, causing a **genuinely different discrete
grid point** to be selected between builds/platforms -- not noise on a
shared point. This matches the observed "different cell, `b2` off by
0.4%" finding exactly, and is the same structural bug pattern as two
already-fixed cases in this history: the VDOSGn truncation-edge
`std::round` tie, and `determineEMinDivKTminRange`'s `i1`/`i2` threshold
search.

### Fix applied: bias the tie-break away from the noise floor

`NCMath.cc`'s `reducePtsByEquidistribution` tie-break was changed to:
```cpp
constexpr double tieBreakRelTol = 1e-9;
if (hi > 0 && (q - vectAt(cum, hi - 1))
    <= (vectAt(cum, hi) - q) + tieBreakRelTol*mtot)
  idx = hi - 1;
```
i.e. only pick the upper neighbour when it is *clearly* (by more than a
fixed fraction of the total cumulative measure, chosen orders of
magnitude above the ~1e-15-level noise a Kahan-summed `cum[]` can carry)
closer, rather than on a raw floating-point comparison. This doesn't
remove the underlying discreteness (some exact tie boundary still
exists), but moves it far away from where realistic cross-platform libm
noise can ever reach it -- the same "push the flip point away from the
noise floor" idea as the earlier VDOSGn taper and `sIntegralAtE` taper,
just implemented as a biased comparison rather than a continuous blend,
since here the output truly is a discrete choice of which input grid
point to keep.

Locally verified: full `sb --long -t` suite green with zero reference-log
changes (this fix doesn't perturb any currently-tested grid), `ncdevtool
check` clean, and bit-for-bit identical vs. an `-mfma`-forced full
rebuild+test (single unrelated failure: `sb_nctestapps_testfmadispatch`,
a benchmark whose ">1.5x speedup" assertion is expected to fail under a
blanket `-mfma` CXXFLAGS override since it also FMA-enables the
`target_clones` "default" variant the benchmark compares against --
confirmed via the test's own "bit for bit" correctness check, which still
passes). The local glibc-hwcap-toggle diagnostic remains unchanged by
this fix (still ~7e-8, same as without it) for the specific failing case
tried, but as established above that diagnostic never reached CI-failure
scale in the first place, so this is not evidence against the fix -- it
is a ceiling of the diagnostic, not a property of the fix. As with the
`sIntegralAtE` taper, this is a real, independently-justified hardening
of a genuinely fragile discrete decision, but whether it resolves the
still-open two-platform CI failure can only be confirmed by an actual CI
re-run.


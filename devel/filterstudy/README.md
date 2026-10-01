# Cross section tables for an NCrystal filter component: findings

Study for mctools/ncrystal#338. The code studied is the new component
`ncrystal_core/src/filter`, with the JSON query `filtertable`, the C API
function `ncrystal_filtertable` (example: `examples/ncrystal_example_filter.c`),
and the Python class `NCrystal.filter.NCrystalFilter` (tests:
`tests/scripts/filtertable.py`, `tests/scripts/filtertable_api.py`,
`tests/src/app_filtertable`, `tests/src/app_capifilter`,
`tests/src/app_filterexample`). A draft of a page for the NCrystal wiki is in
`wiki_filter_tables.md`.

Sections 1-5 describe the study of the first version, with tables for
0.01-100 Aa. Section 6 describes the final design (tables for 0-500 Aa, and
safeguards against wrong tables), and its validation. Some of the scripts of
sections 2-5 use query options of the first version (`wlmin`, `edges`), which
no longer exist.

The scripts in this directory are run against the simplebuild development
build (build it first with `ncdevtool sb`), with `./run.sh <script>.py`. They
write their results (JSON, logs, plots) to this directory.

To visualise tables vs. exact cross sections, use:

```
./run.sh plot_filtertable.py [-w|-e] [-x MIN MAX] [--tol TOL] CFGSTR...
```

(see `--help`).

Files:

| File | Purpose |
|---|---|
| `plot_filtertable.py` | Plot tables vs. exact cross sections, and their relative differences |
| `xsmat.py`, `proto.py` | Helpers: exact cross sections and Bragg edges; Python prototypes of the table algorithms |
| `characterise.py`, `plotcurves.py` | Section 1: properties of the cross section curves |
| `cmp1.py`, `study2.py`, `ycase.py`, `rangeknots.py` | Sections 2-3: table variables, algorithms, tolerances, ranges (Python prototypes) |
| `extrap.py` | Section 5: extrapolation |
| `validate_cpp.py`, `summ_cpp.py`, `ndscale.py`, `coarse.py`, `final_val.py`, `edgecases.py` | Section 4: validation of the C++ implementation |
| `limits.py`, `limit0.py`, `validate_full.py` | Section 6: table end points, and validation of the final design |
| `wiki_filter_tables.md` | Draft of a page for the NCrystal wiki (for users) |

All numbers are for NCrystal 4.4.6 (dev build), with 147 material
configurations: the 134 files of the standard library, plus gas mixtures (air,
He3, BF3), multiphase materials, and variations in temperature (20 K, 80 K,
600 K), density, dcutoff, vdoslux, and with some processes disabled. The
default range is 0.01-100 Aa. The cross section is the total (scattering +
absorption) in barn per atom.

## 1. What the curves look like

* **Absorption** is exactly proportional to lambda (1/v) for the current
  absorption model of NCrystal. The table code does not assume this (it
  evaluates absorption like scattering, via the process objects).
* **Bragg scattering** (powders) is lambda^2 times a step function, with edges at
  lambda = 2d. Materials have up to 29 000 edges in range (CaSiO3), but at most
  about 450 with a jump larger than 0.1% of the total.
* **Inelastic and incoherent elastic scattering** are smooth. For
  hydrogen-rich materials they dominate everywhere.
* **Asymptotics:** at short lambda the total tends to a constant (the free-atom
  cross section). At long lambda, beyond the last Bragg edge, it tends to a + b*lambda
  (incoherent elastic plus 1/v inelastic and absorption). Hydrogen-rich
  materials are not yet fully asymptotic at 100 Aa (log-log slope 0.6-0.9).
* **Special cases:** sigma = 0 (`void.ncmat`); sigma down to ~0.001 barn at short lambda
  when inelastic scattering is disabled; the number of edges depends on
  temperature.

## 2. Table variable: wavelength

Greedy reduction at a tolerance of 1e-3 (number of points):

| Material | (lambda, sigma) | (E, sigma) | (ln lambda, sigma) | (ln lambda, ln sigma) |
|---|---|---|---|---|
| Al | 104 | 232 | 169 | 115 |
| Be 80 K | 280 | 412 | 358 | 297 |
| Polyethylene | 42 | 114 | 70 | 61 |
| B4C | 5 | 244 | 144 | 31 |
| air | 16 | 121 | 76 | 43 |

(Douglas-Peucker numbers for the other columns; Greedy is ~10% below those.)

Linear in lambda is the best or close to the best everywhere. Energy is always the
worst, by a factor of 1.1-50. The reasons are that absorption and the 1/v
parts are exactly linear in lambda, and that the Bragg part (~ lambda^2 between edges) is
closer to linear in lambda than in E (~ 1/E). Log-log is similar to linear-lambda for
crystals, but needs more points for 1/v-dominated materials, and needs
exp/log in the lookup.

## 3. Point selection

* **Existing utilities:** NCrystal has no general curve reduction utility (the
  VDOS thinning is specific to phonon expansions). New code was needed.
* **Adaptive bisection** from a seed grid (as done by e.g. NJOY for ENDF
  linearisation) is not suited: seeding with all Bragg edges keeps up to 58 000
  points for complex materials, and a midpoint-only check misses the
  tolerance in smooth regions.
* **Dense sample + reduction** works well: evaluate the cross section at a
  dense log-spaced grid plus just below and above each Bragg edge, and select
  the subset of points needed for the tolerance.
  * **Greedy forward** (longest valid segment from each point; exponential +
    binary search) needs ~10% fewer points than Douglas-Peucker.
  * **Bragg edges must be included explicitly.** Without them, errors reach a
    factor of 85 near edges (diamond, Be).
  * **Verification at midpoints** between neighbouring dense points (adding
    failing midpoints to the sample and repeating) is needed for a real
    guarantee: without it the tolerance is exceeded by up to 3% at tol=1e-3
    and 20% at tol=1e-4 between the dense points.
  * **The dense sample must scale with the tolerance:** ndense =
    20000*sqrt(1e-3/tol) (the default). Coarse starting samples (2-1000
    points) plus refinement are much cheaper, but miss the tolerance by up to
    16x.
* **Edges dominate:** in crystalline materials, 82-95% of the table points
  are at Bragg edges. A free-knot (non-interpolating) fit of the smooth
  regions, which could save ~30% there, is therefore not worth the complexity.
* **Edges are represented by two table points with the same wavelength.** For
  a lookup "x >= edge uses the value above the edge" (e.g. binary search for
  the last point with x <= lambda), the pair behaves correctly.

## 4. Results of the C++ implementation

With defaults (greedy, auto ndense, edges included), all 147 configurations:

| tol | Worst error | Points: median / 90% / max | Build time: median / max |
|---|---|---|---|
| 1e-2 | 1.001 x tol | 38 / 89 / 134 | 1.2 / 11 ms |
| 1e-3 | 1.006 x tol | 114 / 388 / 936 | 4.2 / 20 ms |
| 1e-4 | 1.009 x tol | 283 / 1139 / 4932 | 15 / 47 ms |

(Build time excludes loading the material, which is cached by NCrystal.)
The worst error is measured at 300 000 random wavelengths plus points 1e-8 and
1e-5 (relative) from each Bragg edge. Tighter tolerances also work (Al at 1e-6:
2312 points in 0.18 s).

The table range has little influence: most points are in the Bragg region.
For example Al has 101 points in 0.01-100 Aa, and 57 in 1-6 Aa.

## 5. Extrapolation beyond the table

Relative error of sigma extrapolated from a table of 0.1-20 Aa:

| Method | Long lambda: median at 40/100 Aa | Long lambda: worst at 100 Aa | Short lambda: median at 0.05/0.01 Aa |
|---|---|---|---|
| Linear continuation in lambda | 0.2% / 0.4% | 6% (Bi) | 0.03% / 0.05% |
| Clamping | 50% / 80% | 80% | 0.4% / 0.7% |
| Power law | 1% / 3% | 36% | 1% / 4% |

Linear continuation in lambda is by far the best, as expected from the
asymptotics. It must be clamped at >= 0 (with inelastic scattering disabled,
sigma vanishes at short lambda and linear extrapolation overshoots a factor of 7).
Since wide tables are cheap, the component can use a wide default range, and
extrapolate only as a fallback.

## 6. Final design: tables for 0-500 Aa, safeguards, and API

**Range.** The tables cover 0-500 Aa, with no wavelength range to choose:

* **Long wavelengths:** extending the table from 100 to 300 Aa costs a median
  of 0 extra points (max 11), and reduces the median error of linear
  extrapolation to 1000 Aa from 0.05% to 0.003% (worst case, liquid water,
  0.8% for both; `limits.py`). Hydrogen-rich materials are much closer to
  asymptotic at 300 Aa (log-log slope 0.83-0.89) than at 100 Aa (0.62-0.74).
  To be conservative, the default is 500 Aa.
* **Wavelength 0:** the first point is the limit for wavelength -> 0, so no
  extrapolation is needed at short wavelengths. The limit is estimated by
  linear extrapolation to 0 of the total cross section (both processes
  evaluated as usual) from h and 2h, which is exact for any a + b*lambda
  behaviour there. This is done at three scales, h = 1e-7, 1e-8 and 1e-9 Aa
  (the estimates agree to 2e-14 or better for all configurations;
  `limit0.py`), and the last estimate is used. Higher order terms only give
  tiny, rapidly decreasing differences: e.g. for a cross section c * lambda^2
  vanishing at 0 (as with UCN scattering only, `ucnmode=only`), the estimates
  are -2c * h^2. The estimates are therefore accepted if they agree within
  1e-3 * tol * max(|limit|, 1e-12 barn), or if they converge (the difference
  between the last two estimates is at most 5% of the difference between the
  first two, and the remaining error is within tol). A cross section without
  a limit (e.g. with 1/lambda, log(lambda) or sqrt(lambda) behaviour) results
  in an exception. (A first version, using two scales and requiring
  agreement within 1e-3 * tol * max(|limit|, 1e-9 * sigma_max), wrongly
  rejected `ucnmode=only`.)
* **Short wavelengths:** a straight line from the limit to the point at
  0.01 Aa is within 2e-5 of the exact cross section for all configurations
  except Al with `inelas=0` (7%), where the scattering vanishes at short
  wavelengths. The dense sample therefore also covers 1e-9-0.01 Aa, so the
  reduction and refinement handle such cases automatically. (It started at
  1e-5 Aa in a first version, but then a straight line from 0 to the region
  around 1e-5 Aa could exceed the tolerance slightly for cross sections
  vanishing at 0 with a large curvature, e.g. 10 * lambda^2 / (1 + lambda^2)
  in the synthetic tests. With the absolute floor of 1e-12 barn, see below,
  the region where the relative tolerance applies extends to even shorter
  wavelengths, hence 1e-9 Aa.)
* **Dense sample:** 10 000 log-spaced points per decade (twice as dense as in
  sections 1-5) from 1e-9 to 500 Aa (117 000 points for tol=1e-3, scaling as
  sqrt(1e-3/tol)).

**Safeguards.** Wrong tables must not be produced silently:

* **Known sharp features:** Bragg edges and the boundaries of the energy
  domains of all physics processes are included as discontinuities (two
  points each).
* **Unknown sharp features:** a feature narrower than the dense sampling
  (e.g. an absorption resonance) can not be found reliably by sampling. The
  physics processes known to have no such features are `NullScatter`,
  `NullAbsorption`, `FreeGas`, `ElIncScatter`, `SABScatter`, `PowderBragg`
  and `AbsOOV` (all 148 configurations use only these). Any other process
  (e.g. the UCN processes, or a future process with resonances) gives a
  warning (`NCRYSTAL_WARN`), until it has been studied and added.
* **Exceptions** are also thrown for cross sections which are negative or not
  finite, for a limit at wavelength 0 which does not converge, for a
  discontinuity found during the refinement (which is not one of the known
  ones), and if the refinement does not converge.
* **Tolerance:** the table is built with 0.98 x tol, and then verified at 4
  random points in each segment (exception if not within tol).
* **Absolute floor:** the tolerance is relative to max(sigma, 1e-12 barn per
  atom), since cross sections vanishing at wavelength 0 (e.g. with neither
  inelastic scattering nor absorption) could otherwise not be tabulated with
  a relative tolerance. (A first version used 1e-9 * sigma_max, which depends
  on the table range. In practice, sigma was below that floor only for
  wavelengths below 5e-7 Aa, for Al with `inelas=0`: for all other
  configurations, sigma / sigma_max is at least 2.6e-5 everywhere.)

**Validation** (`validate_full.py`, 148 configurations, worst error relative
to max(sigma, 1e-12 barn) at 300 000 random wavelengths, log-uniform in
1e-10-500 Aa and uniform in 0-1e-5, 0-1e-7 and 0-1e-9 Aa, and next to all
Bragg edges):

| tol | All OK | Worst error | Points: median / 90% / max | Time: median / max |
|---|---|---|---|---|
| 1e-2 | yes | 0.982 x tol | 40 / 89 / 135 | 7 / 22 ms |
| 1e-3 | yes | 0.987 x tol | 116 / 389 / 947 | 24 / 60 ms |
| 1e-4 | yes | 0.985 x tol | 290 / 1148 / 4974 | 109 / 273 ms |

(The times are larger than for the first version of the final design, with the
dense sample starting at 1e-5 Aa: 4, 14 and 56 ms median.) Materials with UCN
scattering (`ucnmode=remove` and `ucnmode=only`), which give a warning, also
have valid tables (worst error 0.98 x tol), and oriented materials are
rejected.

The synthetic tests (`tests/src/app_filtertable`) check that oscillating
curves, resonance-like peaks (1% wide in energy), declared steps and curves
vanishing at wavelength 0 (also with a large curvature) are tabulated within
the tolerance, and that undeclared steps, curves without a limit at wavelength
0, and negative, NaN or infinite values result in exceptions.
`tests/src/app_filterexample` checks that the evaluation function of
`examples/ncrystal_example_filter.c` (a verbatim copy, checked by
`ncdevtool check copies`) agrees with the internal evaluation of the tables.

**API:**

* C: `ncrystal_filtertable(cfgstr, &n, &wl, &macroxs, options)` returns the
  table of the macroscopic cross section in 1/cm, with a tolerance of 1e-3
  (the options are reserved for future extensions, e.g. a tolerance, and must
  be NULL or empty). `examples/ncrystal_example_filter.c` shows how to
  evaluate it, with a function suitable for GPUs (linear extrapolation of the
  last segment beyond 500 Aa, clamped at 0).
* Python: `NCrystal.filter.NCrystalFilter(cfgstr, options=None)` creates the
  table (available via `.table`), and `.xsect(ekin=..., wl=...)` evaluates it
  (also for arrays).
* The JSON query `["filtertable", CFGSTR, "tol=...", "wlmax=...", ...]` gives
  the table in barn per atom, with diagnostics (e.g. the leaf processes).

## 7. Open questions and next steps

1. **"NPTS" mode:** the algorithm is tolerance driven. A maximum number of
   points could be supported by searching for the smallest tolerance that
   satisfies it (each table only takes milliseconds).
2. **Where the code lives:** the component is small (~600 lines); it could be
   moved to `extd_utils` later, as you suggested.
3. **Performance:** the refinement re-evaluates all midpoints in each
   iteration; caching the already verified ones would save some time, but the
   current cost (tens of ms) is already small.
4. **More processes:** the UCN processes and `SANSSphereScatter` could be
   studied and added to the supported processes.
5. **Local build notes:** gcc 15 from conda-forge reports a
   `-Werror=maybe-uninitialized` false positive in the existing
   `NCSANSUtils.cc` (worked around with `CXXFLAGS=-Wno-error=maybe-uninitialized`),
   and the latest ruff (0.16) reports 1305 issues in existing Python code, so
   `ncdevtool check` fails at the ruff step. Also, with Python 3.14, argparse
   colours its help output, which makes the tests of the command line tools
   fail (set `NO_COLOR=1`). None of this is related to this work.

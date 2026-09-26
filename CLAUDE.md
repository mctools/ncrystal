# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

NCrystal is a C++ (C++11+, no third-party deps) library for thermal neutron scattering in crystals and other materials, with C and Python interfaces, CLI tools (`nctool`, `ncrystal-config`, ...) and `.ncmat` data files (format: `docs/ncmat_doc.md`). Materials are selected with cfg-strings such as `Al_sg225.ncmat;temp=10C`.

## Branch

Work happens on the `tk_volatile` branch (not `main`): a long-running development branch, currently version 4.4.7 (work in progress) and far ahead of `main` (new `sab`/`sabscatter`/`vdos` code, new tests such as `app_fft`, `app_tinyvect`, `app_vdosutils`). Only `basictest.yml` remains under `.github/workflows` here.

FIXME: This branch instruction should obviously not be in the version merged onto the main branch.

## Everyday workflow

Edit code (mostly `ncrystal_core/`, `ncrystal_python/`, `tests/`), then build and test with simplebuild via `ncdevtool`, inspect results, repeat:

```
./devel/bin/ncdevtool sb -t                                   # incremental build + all tests
./devel/bin/ncdevtool sb -t --testfilter='*cfgparse*'         # only matching tests
./devel/bin/ncdevtool sb -t --testfilter='*a*,*b*,!*c*' --testexcerpts 30   # comma-separated, "!" negates, show failing logs
./devel/bin/ncdevtool sbrun <cmd> [args]                      # build, then run cmd in the sb environment (do not pass -t here)
```

- Needs a venv with `pip install -r devel/reqs/requirements_all_and_devel.txt` (provides `simple-build-system`, numpy, ASE, ruff, ...). `source devel/setup.sh` puts `ncdevtool` in PATH. If the venv lives inside the repo, put it under `venv/` or `nocheck/venv/` (both gitignored and excluded from `ncdevtool check`'s repo-file scans); elsewhere gets flagged (e.g. bogus "bad license blurb" hits on the venv's own site-packages).
- To build with clang instead of the system default, `export CC=clang CXX=clang++` before `ncdevtool sb -e` (forces re-examination, rebuilds the same cache with the new toolchain -- no separate cache needed, and later `sb`/`sb -t` runs stay incremental).
- First build takes ~1 min, later ones only rebuild what changed; a filtered test run takes seconds.
- `--testfilter` matches **simplebuild test names**, not ctest names: `tests/src/app_<x>/` is `sb_nctestapps_test<x>`, `tests/scripts/<x>.py` is `sb_nctestpynp_test<x>`. A `.log` next to a test is its reference stdout; any diff fails the test. Test binaries can be run directly with `ncdevtool sbrun sb_nctestapps_test<x>`.
- `sb --long -t` (`--long` must be first) also enables `long_*` tests. `sbdbg`/`sbrundbg` use a Debug build in a separate cache, so `ncdevtool sbdbg -t [--testfilter=...]` is the Debug equivalent of `sb -t`.
- In the sb environment the Python module is `NCrystalDev`, e.g. `ncdevtool sbrun python3 -m NCrystalDev nctool --help`; `ncrystal-config` is `sb_nccmd_config`.
- Other useful: `sb -s` (summary), `sb --grep/--find PATTERN`, `sb --pkginfo PKG`, `sb -c` (wipe cache), `sb -i` (rebuild from scratch). Docs: https://mctools.github.io/simplebuild/
- How it works: `ncdevtool sb` runs `devel/simplebuild/sbgenerate.py`, generating git-ignored `devel/simplebuild/autogen/` full of simplebuild packages that are **symlinks to the real sources**. Always edit the real files. New/moved files are picked up on the next run.
- `ncdevtool check` is the quick (~2 s) linter step: run it on its own before the costlier `sb -t`, and always before committing. It runs static checks (`ncdevtool check -l` lists them: `comps`, `copies`, `deps`, `fixme`, `incguards`, `license`, `misc`, `ruff` (ruff linter on all Python code), `toml`, `versions`, `yaml`). Run a subset by name (`ncdevtool check ruff`) or all-but with `-n`. On this branch the marker-word check fails by design, so use `ncdevtool check -n "fix""me"` (spelled so; this is what CI runs) to see the other checks.
- Other modes (`ncdevtool` lists all): `grep|find|replace`, `cppcreate|cppmove|fixdeps|cppana|graph`, and `cmake` (pure CMake build+ctest in a temp dir, slow; options `-b`, `--dbg`, `--long`, `-c @-DVAR=VAL`, `-t @-R@<regex>`).
- FIXME (temporary, not yet decided whether to keep; drop this bullet along with the commit that added the `nofpc` modes if it goes): `sb[env|run][dbg]nofpc` variants of the `sb`/`sbenv`/`sbrun` modes (each with its own cache) add `-ffp-contract=off`, to check whether a numerical difference between platforms is caused by the compiler fusing `a*b+c` into FMA instructions (default off; a `-mfma`-built binary otherwise fuses, a baseline-x86-64 one doesn't). The equivalent CMake option is `NCRYSTAL_DISABLE_FPCONTRACTION` (default OFF).
- Non-test, developer-only benchmark/utility apps (`app_benchsab`, `app_benchsabsample`, `app_minimc`, `app_minimcbench`, `app_query`) live under `devel/simplebuild/static/NCDev/` rather than `tests/src/`, and are run with `ncdevtool sbrun sb_ncdev_<name>` (not built or picked up by `-t`).

## Repository layout

Since 4.0 the distribution is several packages, always released with identical versions (in `VERSION` and stamped into several pyproject/CMake files, hence the separate "Update version numbers" commits):

- `ncrystal_core/`: C++ shared library, headers, `ncrystal-config`. CMake project (options in `cmake/modules/ncrystal_options.cmake`; tests need CMake >= 3.28), also buildable via pip/scikit-build-core.
- `ncrystal_python/`: `NCrystal` Python module and CLI tools (`_cli_*.py`). Talks to core via the C API (`ncrystal.h`, ctypes in `_chooks.py`); no build-time dependency on core.
- `ncrystal_metapkg/`, `ncrystal_pypluginmgr/`, `ncrystal_verify/`: meta-package, plugin manager, and the `ncrystal-verify` package. The latter is generated at build time from the pure-Python scripts in `tests/scripts` (rewriting `NCrystalDev` to `NCrystal`, skipping ones using `NCTestUtils.loadlib`/`hists`), so end-users can run them on an installed NCrystal. Scripts declare optional deps with a `# NEEDS: numpy ...` header line, which must match the `[all]` extras in its pyproject. Rarely relevant to everyday work.
- `data/` (`.ncmat` files, embedded into the library by default), `tests/`, `examples/`, `devel/`.

## Python package (`ncrystal_python/src/NCrystal/`)

- Public modules (`core`, `datasrc`, `cfgstr`, `ncmat`, `vdos`, `plugins`, `cli`, ...) are re-exported by `api.py`, which `__init__.py` star-imports (unless `NCRYSTAL_SLIMPYINIT` is set). Private modules start with `_` (`_chooks.py` = ctypes hooks into the C API, `_locatelib.py` finds the shared library, `*impl.py` = implementation of public modules). Numpy is optional at import time and accessed via `_numpy.py`. Python code uses `print`/`warn` from `_common` (redirectable) rather than the builtins, and `ncgetenv*` for env vars.
- **Command-line tools are ordinary modules, not separate scripts.** Each tool is `_cli_<name>.py` (e.g. `nctool`, `ncmat2cpp`, `cif2ncmat`) with `climod_metadata()` (display group/order/description for the `ncrystal` overview), `create_argparser_for_sphinx(progname)`, and `main` decorated with `@cli_entry_point` (from `_cliimpl.py`, which also provides `create_ArgumentParser`, to be used instead of `argparse.ArgumentParser` directly). The decorator turns exceptions into clean `ERROR:` exits when run from a shell, and when called from Python it lets `NCrystal.cli.run('ncmat2cpp', ...)` call the same `main` in-process (mapping `SystemExit` to `RuntimeError`).
- Tools are **discovered by globbing `_cli_*.py`** (plus the special `config` tool, `_cliwrap_config.py`, which just runs the compiled `ncrystal-config` from core). Canonical names are `ncrystal_<name>` (except `nctool` and `ncrystal-config`). They are reachable as `python3 -m NCrystal <tool>`, as the umbrella `ncrystal <tool>` command (`_clientry.py`), and as standalone commands.
- Standalone commands are wired up in three places that must agree: the `[project.scripts]` table in `ncrystal_python/pyproject.toml` (`ncrystal_<name> = "NCrystal._cli_<name>:main"`), the same table in the root monolith `pyproject.toml`, and, for simplebuild, `sbgenerate.py`, which auto-generates a `sb_nccmd_<name>` wrapper script per `_cli_*.py`. `ncdevtool check` verifies that `[project.scripts]` has exactly one entry per `_cli_*.py`, so **adding a tool means: create `_cli_<name>.py` and add the pyproject entries** (nothing else needs registering).

## Core C++ architecture

Code is split into small components, one directory each in `ncrystal_core/src/<component>/`, with public headers in `ncrystal_core/include/NCrystal/`.

- **Each component's allowed dependencies are listed in its `dep.txt`** and enforced by the build. After adding cross-component `#include`s, update it (`ncdevtool fixdeps` does this).
- Data flow: cfg-string + data source -> immutable `Info` (built by factories, e.g. `.ncmat` parsing in `ncmat/`) -> `Scatter`/`Absorption` objects from scatter/absorption factories (`stdscatfactory`, `absfact`, ...). Factories are registered by name and extendable by plugins.
- The C API (`cinterface/`, `ncrystal.h`) is the stable ABI used by the Python module and by Geant4/OpenMC/McStas integrations; changes there usually need matching edits in `_chooks.py`. `NCrystal.hh` is the C++ entry header.
- **All public headers are stable ABI, not just the C API**: anything under `ncrystal_core/include/NCrystal/` that is *not* under `ncrystal_core/include/NCrystal/internal/`. Changes there (including to widely-used utility templates like `Pimpl` in `NCDefs.hh`) must preserve existing layout/symbols; only headers under `internal/` are free of that constraint.

## C++ style for new code

Target C++11, two-space indent, max 72 columns, `//` comments only (never `/*..*/`), `m_` prefix for data members, short variable names. Code should be concise and efficient, and give identical results across platforms (no `long double`). Project idioms:

- Types: `VectD` (`std::vector<double>`), `PairDD`, `kInfinity`, `Span` (C++11 `std::span`, accepts a `VectD`), `Optional<T>` (like `std::optional`, use `NullOpt`, `.emplace()`, `.has_value()`, `.value()`).
- `NCDefs.hh` and `NCMem.hh` (pulled in essentially everywhere indirectly) already include a large set of standard headers (`<limits>`, `<cmath>`, `<vector>`, `<algorithm>`, `<memory>`, `<type_traits>`, ... -- check them before assuming a `#include` is needed). Don't add a redundant explicit include for something they already provide.
- Use `ncmax`, `ncmin`, `ncabs` instead of `std::fmax/fmin/fabs`.
- Use `nc_assert(..)` for internal invariants rather than throwing explicit exceptions (it is compiled out with `NDEBUG`; `nc_assert_always` is not). Errors caused by user input or data are the exception: they use `NCRYSTAL_THROW2(BadInput, "..." << x)` (also `CalcError`, `DataLoadError`, `MissingInfo`, `FileNotFound`), and `LogicError` is for internal errors.
- Index vectors (not spans) with `vectAt(x,i)` (plain `x[i]` in non-debug builds); `x.front()`/`x.back()` are fine.
- Prefer range-for (`for ( auto& f : foo )`); when an index is really needed use `for ( auto i : ncrange(n) )` (also `ncrange(lo,n)`), not a hand-written `size_t` loop.
- Start functions with an `#ifndef NDEBUG` block validating all inputs (using `nc_assert_always` inside), so validation does not clutter the rest.
- New non-trivial functions get a ~6 line top comment describing inputs and outputs.
- Comments (and commit messages) should be terse: enough for a clever developer or AI to reconstruct the reasoning, not a full explanation for a junior dev. Prefer one dense line over a paragraph.
- Use `fmt(x)` for number formatting in messages/output, `StableSum`/`StableSumKahan`/`StableDbl` where summation precision matters, and `floateq(a,b,rtol=1e-6,atol=1e-6)` to compare doubles.
- `NCRYSTAL_FMADISPATCH_ATTR` (`target_clones`-based FMA dispatch): see `doc/devel_fma_attribute.md` for the rules before decorating a function with it. Notably: never `ncrange(...)` inside such a function (confirmed Clang codegen bug, rule 6).
- Clang `-O3` can optimise a `cond ? f(x) : v` ternary into an unconditional call to `f(x)` (incl. for pole/domain-error inputs) and will re-derive that even past a naive "route the bad input through a safe value first" rewrite; only a `volatile` barrier reliably blocks it. Relevant for anything like `s>0.0 ? std::log(s) : 0.0` -- see `NCrystal::SABUtils::safeLogOrElse` for the pattern.
- Ownership: prefer `shared_obj<T>` (a `shared_ptr` wrapper that guards against nullptr, except after moves, so a `shared_obj<const A>` converts directly to a `const A&` argument) over `std::shared_ptr`, `std::unique_ptr` for unique ownership; use the `MoveOnly` / `NoCopyMove` / `NonInstantiable` helper bases, `final` on classes not meant for inheritance, and `nc_as_const`. Use `NCRYSTAL_MSG`/`NCRYSTAL_WARN` for messages, never iostream directly in library code.
- Unit tests (`tests/src/app_<x>/main.cc`, expected stdout in `test.log`) may define `#define REQUIRE(x) nc_assert_always(x)` and `#define REQUIREFLTEQ(x,y) nc_assert_always(floateq(x,y))`.
- All C++ (including test/scratch apps under `tests/src/`) builds with `-Wall -Wextra -pedantic -Werror` (`/W4 /WX` on MSVC), so even throwaway diagnostic code must be warning-clean: no unused variables, no misleading indentation (e.g. two statements after one `if`/`for` on the same line), and no brace-init list (`Type{a,b}`) as a bare macro argument (the comma inside it splits into two macro arguments; assign to a local first).
- `ncdevtool check` (`misc`) rejects any repository file above ~60KB (300KB for `.log`), unless it is added to the override table in `devel/pypath/ncrystal_repo_tools/_check_misc.py`. The point isn't only to keep individual files small: it's to keep an eye on overall repo size. A large addition of code or data is not automatically wrong, but it should be a deliberate choice that's judged worth it, never something that happens by accident (splitting a test into a new `app_<name>` is only the mechanical fix for one file; it doesn't itself address that judgement).

File layout conventions:

- Header guards are `NCrystal_<Name>_hh` for `NC<Name>.hh` (checked by `ncdevtool check`), file names are `NC<Name>.hh/.cc`, and each file starts with the licence banner then the guard.
- Code lives in `namespace NCRYSTAL_NAMESPACE { ... }`. In `.cc` files, `namespace NC = NCrystal;` is declared at the top, and out-of-namespace definitions are written as `NC::Foo::bar(...)`. File-local helpers go in an unnamed namespace.
- Small classes are typically header-only in a trailing "Inline implementations" section of the header (with `inline`), and only larger or non-hot code goes in the `.cc`. Public headers of a component go in `include/NCrystal/internal/<component>/` (or `core`, `interfaces`, ... for public API), private ones next to the `.cc` in `src/<component>/`.
- Section banners in headers use `////` comment blocks.

## Conventions

- Every source file carries the licence banner (checked by `ncdevtool check`).
- The marker word spelled "fix"+"me" (any capitalisation) is used for work-in-progress notes on this branch, but the marker check in `ncdevtool check` flags every occurrence, so it must all be resolved before merging to `main`. Don't add new ones unless asked.
- User-visible changes get an entry at the top of `CHANGELOG`.
- CI covers many compilers/platforms and Python 3.8-3.14, so stay within C++11 in core and Python 3.8 in the Python package.

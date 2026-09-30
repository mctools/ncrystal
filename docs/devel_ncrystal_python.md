# ncrystal_python: developer notes

Working notes on the `NCrystal` Python package (`ncrystal_python/`, ~24k
lines in `src/NCrystal/`). Keep below 200 lines. Commits touching this file
must contain only this file, with the message
`updating docs/devel_ncrystal_python.md`.

## Packaging and install modes

- `ncrystal_python/pyproject.toml`: setuptools, `ncrystal-python`, version
  read from `NCrystal.__version__`, depends on `numpy>=1.22` (so numpy is a
  hard *install* dependency, although code treats it as optional at import).
  `[project.scripts]` maps `ncrystal` -> `_clientry:main` and each
  `ncrystal_<x>`/`nctool` -> `_cli_<x>:main`.
- Deliberately has no dependency on `ncrystal-core`; the shared library is
  found at runtime (`_locatelib.py`).
- Install mode is signalled by empty marker files next to the modules:
  - `_is_std.py`: normal install (checked in, empty).
  - `_is_monolithic.py`: root `pip install .` (root `CMakeLists.txt` renames
    `_is_std.py`). Both present = broken env -> `SystemExit`.
  - `_is_sblddevel.py`: simplebuild dev mode (added by
    `devel/simplebuild/pypath/sbgen/main.py`, which symlinks all sources into
    package `NCrystalDev`). Changes lib lookup (`sb_nccmd_config`), shell
    command names (`sb_nccmd_<x>`, `sb_nccmd_tool`) and lets
    `_common._lookup_existing_file` resolve `pkg/file` via `$SBLD_DATA_DIR`.
  - `_do_namespace_envvars.py`: env vars become `NCRYSTAL<NS>_X` (see
    `_common.expand_envname`); also set by sbgen.
- Tests import `NCrystalDev`; `ncrystal_verify/generate.py` rewrites that to
  `NCrystal` for the `ncrystal-verify` package.

## Import chain and library loading

- `__init__.py`: version metadata, py>=3.9 guard, then
  `from .api import *` unless `NCRYSTAL_SLIMPYINIT` is set (never namespaced).
  `version_tuple`/`version_num` defined after.
- `api.py` star-imports `exceptions`, `core`, `datasrc`, `_testimpl`
  and selected names from `constants`, `atomdata`, `cfgstr`,
  `ncmat`, `plugins`, `vdos`. `cifutils`, `misc`, `mcstasutils`, `plot`,
  `minimc`, `hist`, `ncmat2endf`, ... must be imported explicitly.
- Importing `core` or `datasrc` loads the C library immediately
  (`_rawfct = _get_raw_cfcts()` at module level), so plain `import NCrystal`
  needs a working lib. `core` also installs the default C++ message handler
  (`_msg._setDefaultPyMsgHandlerIfNotSet`).
- `_locatelib._search()` order: `NCRYSTAL_LIB` (+ optional
  `NCRYSTAL_LIB_NAMESPACE_PROTECTION`, both never namespaced) ->
  `_ncrystal_core[_monolithic].info` module (skipped in sb mode) ->
  `ncrystal-config --show shlibpath namespace version`. Version must match
  exactly. `NCRYSTAL_DEBUG_LIBSEARCH` prints the search.
- `_chooks._load()` wraps the C API with ctypes into a dict of callables
  (`_get_raw_cfcts()`), keyed by C name or by a python-side helper name
  (`iter_hkllist`, `raw_vdos2knl`, `jsonquery`, `flexmmcrun`, ...).
  `_wrap()` applies the namespace prefix (`ncrystal<ns>_`), and after each
  call checks `ncrystal_error()` and raises the matching `NC*` exception
  (`exceptions.py`; all derive from `NCException(RuntimeError)`).
  Library runs with halt-on-error off and quiet-on-error on.
- Adding a C API function: declare in `ncrystal.h`, then `_wrap(...)` in
  `_chooks._load` (plain) or add a `hide=True` raw wrap + python helper
  registered in `functions`. Arrays go through `as_contiguous_double_array` /
  `ndarray_to_dblp`; library-allocated results are copied with
  `_cptr_to_nparray` and freed via `ncrystal_dealloc_*`. Python callbacks
  passed to C must be kept alive (`_keepalive`) if C++ stores them.
- `_numpy.py`: `_np` (or None), `_ensure_numpy()`, endpoint-exact
  `_np_linspace/_np_geomspace/_np_logspace`, `_np_trapezoid` (numpy 1/2).

## Main public modules

- `core.py` (1.9k): `RCBase` (holds raw handle, unrefs in `__del__`),
  `AtomData`, `Info` (phases, composition, structure, `AtomInfo`, HKL lists,
  `DynamicInfo` subclasses `DI_Sterile/FreeGas/ScatKnlDirect/VDOS/VDOSDebye`,
  custom sections), `Process` -> `Scatter`/`Absorption`, `LoadedMaterial`
  (`load()`), `directLoad`, `TextData`, plus `createInfo/...`,
  `clearCaches`, `enableFactoryThreads`, `setDefaultRandomGenerator`.
  AtomInfo/DynamicInfo hold weakrefs to their `Info` and raise if it died.
  Scalar calls go to the C scalar function; array/`repeat` calls to the
  `_many` variants (numpy). Oriented `crossSection` on arrays is vectorised
  in Python (slow). Deprecated methods warn (`genscat`, `getCalcName`, ...).
- `datasrc.py`: search dirs, in-memory files (`registerInMemoryFileData`),
  enabling/disabling data sources, `browseFiles`.
- `cfgstr.py` (`normaliseCfg`, `decodeCfg`, `generateCfgStrDoc`),
  `atomdata.py`, `constants.py` (unit conversions), `plugins.py`.
- `ncmat.py` + `_ncmatimpl.py` (2.2k): `NCMATComposer`, the builder behind
  `cif2ncmat`, `hfg2ncmat`, `vdos2ncmat` etc. (`from_cif/from_info/...`,
  `set_*`, `create_ncmat`, `write`, `register_as`, `load`).
- `vdos.py` (1.1k): Debye VDOS/MSD helpers, `analyseVDOS`, `extractGn`,
  `extractKnl` (-> C `ncrystal_raw_vdos2kernel`), `PhononDOSAnalyser`.
- `cifutils.py` (1.8k): `CIFSource`, `CIFLoader` (gemmi, spglib, ase;
  online DBs: COD, Materials Project via `MATERIALSPROJECT_USER_API_KEY`).
- `ncmat2endf.py` + `_ncmat2endf_impl.py` (1.6k): ENDF export via
  `endf_parserpy`.
- `misc.py` + `_miscimpl.py`: type-erasure wrappers `MaterialSource`,
  `AnyTextData`, `AnyVDOS`; `evaluate_query` (JSON queries to C++,
  `huge_arrays` fast path through `ncrystal_fill_jsonarray`).
- `minimc.py`, `minimc_objects.py`, `_mmc_impl.py`, `_mmc_doc.py`: MiniMC
  runs via C `ncrystal_flexmmcrun`, with results/tallies as `Hist1D`
  (`hist.py`). `mmc.py`/`_mmc.py` are obsolete shims.
- `plot.py` (matplotlib), `mcstasutils.py`, `hfg2ncmat.py` + `_hfgdata.py`,
  `ncmat2cpp.py` + `_ncmat2cpp_impl.py`, `_hklobjects.py`, `_sabutils.py`.
- `_testimpl.py`: `NCrystal.test()` / `python -m NCrystal.test` smoke tests.

## Command-line tools

- One module per tool, `_cli_<name>.py`, with `climod_metadata()`
  (displaygroup `main|conv|misc`, displayorder, descr ending in "."),
  `create_argparser_for_sphinx(progname)`, and
  `@cli_entry_point def main(progname, arglist)`. Tools: nctool, browse,
  minimc, cif2ncmat, hfg2ncmat, ncmat2hkl, ncmat2endf, endf2ncmat, vdos2ncmat,
  verifyatompos, query, ncmat2cpp, mcstasunion; plus `config`
  (`_cliwrap_config.py`, runs the compiled `ncrystal-config` in a subprocess).
- Discovery = glob of `_cli_*.py` (`_cliimpl.cli_tool_list_impl`). Names:
  short (`ncmat2cpp`), canonical (`ncrystal_ncmat2cpp`, but `nctool`,
  `ncrystal-config`), shell command (differs in sb mode).
- Entry routes: standalone scripts, `ncrystal <tool>` (`_clientry.py`, which
  also prints the grouped overview), `python -m NCrystal <tool>`
  (`__main__.py`; `--unblock` runs as `cli.run`), and `NCrystal.cli.run()`.
- `cli_entry_point`: from a shell, turns `NCException` into
  `"<Type> ERROR: msg"` `SystemExit`, other exceptions into `"ERROR: msg"`,
  and shows `NCrystalUserWarning` as `WARNING:` lines. The hidden
  `--show-exceptions` flag keeps the traceback. In `cli.run` mode, argparse
  errors raise instead of exiting, argparse output goes through NCrystal
  print, and `SystemExit` becomes a clean return or a `RuntimeError`.
- Parsers must come from `_cliimpl.create_ArgumentParser`, which also
  patches argparse so `--help` / "invalid choice" output is identical on
  py3.9-3.14 (reference logs depend on this).
- `browse` (`_cli_browse.py`, new in 4.4.7) is meant to take over nctool's
  `--browse/--extract/--plugins` (kept for now). It is a thin CLI over the
  public `browse.py` API (`DataBrowser`: chainable immutable selections,
  lazy per-factory physics loading, all formatting; `query_data`; element
  types `DataEntry`, `PhysicsProps`; `AtomDBBrowser`/`query_atomdb` for
  `["util","atomdb"]`), built on the C++ JSON queries
  `["util","browsedb",FACT,(I,N,)("cheap")]`, `["util","browsefactories"]`
  and `["util","factorythreads"]` (`src/query/NCBrowseQuery.cc`). Loading
  runs in parallel via `FactoryJobs`, with temporary threads ("nthreads=N"
  query arg) only if the user did not configure factory threads (explicit
  call or `NCRYSTAL_FACTORY_THREADS`, which the thread pool reads on first
  use). CLI adds grep-like `--color` (`_use_color`), tables
  (`--columns/--sort`), `--json`, `--info`, `--count`, `--path`, and name
  suggestions. `nctool --help` points to it (keep nctool options for long).
  Tests: `browsedbquery.py`, `browseapi.py`, `clibrowse.py`.
- New tool: add `_cli_<x>.py` + `[project.scripts]` entry in
  `ncrystal_python/pyproject.toml` and the root `pyproject.toml`
  (`ncdevtool check` enforces the match). sbgen creates `sb_nccmd_<x>`.

## Conventions (Python side)

- Python 3.9 syntax only; ruff runs via `ncdevtool check ruff` with a long
  ignore list (`devel/pypath/ncrystal_repo_tools/_check_ruff.py`).
- Output through `_common.print` (redirectable with
  `set_ncrystal_print_fct`, `modify_ncrystal_print_fct_ctxmgr`,
  `capture_print_ctxmgr`) and `_common.warn`
  (`NCrystalUserWarning`). `_cli_*` modules import both from `_cliimpl`.
- Env vars: `ncgetenv[_bool|_int|_int_nonneg]('X')`/`ncsetenv` use
  `expand_envname('X')`; never pass the `NCRYSTAL_` prefix. Namespace
  protection (`NCRYSTAL_NAMESPACE`) only renames libs/symbols, so several
  builds can coexist; env vars stay `NCRYSTAL_X` for end users. Only the
  sb dev build also namespaces env vars (`NCRYSTALDEV_X`), via the C++
  define `NCRYSTAL_NAMESPACED_ENVVARS` + Python `_do_namespace_envvars.py`
  (kept in sync by sbgen). Must match C++ `raw_getenv` (`NCString.cc`).
- Heavy or optional imports are done inside functions to keep import fast.
  Module state lives in 1-element lists (`_cache = [None]`) instead of
  `global`.
- Files are written with `_common.write_text` (utf8, `\n`); tests fake the
  current time with `FixedFakeDatetimeNow`.

## Tests

- `tests/scripts/*.py` (~100 scripts, `.log` = expected stdout, run as
  `sb_nctestpynp_test<x>`) use `import NCrystalDev as NC` and helpers from
  `tests/pypath/NCTestUtils` (`env.ncsetenv`, `enable_fpe`, `printnumpy`,
  `loadlib`, ...). A `# NEEDS: numpy ...` header declares optional deps.
- Quick loop: `ncdevtool sb -t --testfilter='sb_nctestpynp_*'`.

## Review findings (2026-09-25)

All fixed on `tk_volatile_ub24` (details in commit messages and CHANGELOG):
`ncpprint(do_sort)`, `_setMsgHandler(None)`, MiniMC callback `'error'`,
`NCRYSTAL_LIB` namespace inference, builtin `print` uses, failed-lib-load
retry, `ncsetenv` namespacing, `str(AtomInfo)` with missing dt/msd,
unknown dyninfo types, docstrings; `obsolete.py` removed. Threaded use
exposed two C++ bugs, also fixed: static file read buffer
(`app_mtfileread`) and global C API error state (`app_capierrmt`).
Lessons: verify by running code, and check the C++ side too.

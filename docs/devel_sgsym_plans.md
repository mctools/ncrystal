# The sgsym component: status and plans

Developer notes on the internal `sgsym` component
(`ncrystal_core/src/sgsym/`, headers in
`ncrystal_core/include/NCrystal/internal/sgsym/`), which provides space group
settings and symmetry operations. Nothing outside `sgsym` uses it yet.

## Status

- **Settings:** `SpaceGroupHallNumber` (a strongly typed 1..530 index, i.e.
  the row of ITVB Table A1.4.2.7) and the immutable `SpaceGroup` class
  (number, ITVB choice code, Hall symbol, strict string parsing like
  `"227:2"`, where a bare number selects the default setting). The hardwired
  table of the 530 settings was extracted from and verified against ITB
  (2010), the SgInfo web tables, spglib, gemmi and ASE (see the provenance
  comment in `NCSGSymHallTable.cc`).
- **Symmetry operations:** `SymOp` (12 bytes, exact). Rotation parts are
  restricted to the 64 rotations of the 530 settings (signed permutations,
  and 6/mmm in hexagonal axes), translations are stored in units of 1/24
  (modulo 1). Composition checks its result, and strings are CIF-style like
  `-x,y+1/2,-z+1/2` (fractions only, no decimals).
- **Operations of each setting:** `SGSymmetry::get(sg)`, derived at run time
  from the Hall symbol and cached lazily. It provides the operations (in a
  canonical order: sorted coset representatives with the identity first,
  then one block per centring vector), representatives, centring vectors,
  order, lattice symbol and centrosymmetry. The Hall parser is only reachable
  (for tests) as `detail::rawSymOpsFromHallSymbol` in the `.cc`.
- **Derived properties:** `SpaceGroup::crystalSystem()` (`SGCrystalSystem`,
  since the name `CrystalSystem` is taken by `NCLatticeUtils.hh`) and
  `settingInfo()` (hexagonal-family axes, monoclinic unique axis and cell
  choice, origin choice, orthorhombic axis permutation), plus
  `SGSymmetry::cellConstraints()`, derived from the operations, with
  `check()` and `complete()` for cell parameters (exact equality, today's
  completion rules but driven by the setting).
- **Site expansion:** `expandSiteToOrbit(sym, site)` returns an `SGSiteOrbit`
  (positions wrapped into [0,1) with the site first, and the site symmetry
  order). Positions are compared per fractional coordinate, periodically.
  Images within 5e-4 are merged and the site is symmetrised (averaged over
  the operations fixing it, so e.g. `0.3333, 0.6667` becomes exactly
  (1/3, 2/3)). Images which are not merged, but are within 1e-2, give a
  BadInput about the site being close to a special position but not on it
  (catching input like `0.333`). The thresholds are based on an analysis of
  rounding noise and genuinely distinct positions in CIF files from the COD.

Everything is verified against spglib, gemmi, ASE, the ITB table and
`EqRefl`, by the tests `app_sgsym` and `sgsym.py` (the latter also
reconstructs all stdlib crystals with space groups from one site per orbit).

## Open items

- **NCMAT integration:** `@SPACEGROUP` with setting suffixes (e.g. `227:2`)
  and a new section with one site per orbit, in a new NCMAT format version.
  Expansion would happen in `NCMATData` before `validate()`, so `@DYNINFO`
  fraction checks see the expanded counts.
- **CIF export** (e.g. as a JSON query): needs Hermann-Mauguin symbols in the
  settings table (step 4, to be verified against several sources like the
  Hall symbols), and probably grouping positions into orbits and detecting
  the setting of existing structures (`Info` only stores the number).
- **Deferred derived properties:** point group classification (incl.
  orientation, e.g. `-3m1` vs. `-31m`), Bravais type, chirality, and
  transformations between settings (e.g. origin shifts 227:1 to 227:2).
- **Reciprocal space utilities** (step 5): reflection conditions (systematic
  absences) and Laue orbits of hkl, possibly eventually replacing `EqRefl`.
- **Migration:** move `crystalSystem()`, `checkAndCompleteLattice[Angles]()`,
  `isRhombohedralSpaceGroup()` and `usesRhombohedralAxes()` from
  `NCLatticeUtils.hh` into `sgsym` (TODO comments in the code).
- **Info object (NCrystal 5.x):** decide what to store about space groups.
  Not all sources know the setting (LAZ/LAU files, NCMAT files with a bare
  number, plugins), so one idea is to always store the number, plus an
  optional `SpaceGroupHallNumber` when the setting is specified or
  unambiguously determined.

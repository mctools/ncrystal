# Roadmap: partial site occupancy

Status: **postponed to NCrystal 5.0**, since it requires C++ ABI changes (see
below) and a new NCMAT format version. This note records the agreed design so
work can resume from it.

## Motivation and current state

Neither the NCMAT format nor the C++ core supports partial site occupancy. The
Python `NCMATComposer` (and hence `ncrystal_cif2ncmat`) emulates an occupancy
*o* with an `@ATOMDB` mixture: *o* of the real element plus (1-*o*) of a
hijacked isotope (`Og299`, `Og298`, ...) acting as a "sterile" atom with the
mass of the element and zero cross sections (see
`allow_siteoccu_ncmatv5_hack` in `_ncmatimpl.py`).

- Scattering is correct: the coherent scattering length becomes o*b, and the
  mixture variance term 4pi*o(1-o)*b^2 is exactly the disorder ("Laue
  diffuse") scattering of randomly placed vacancies.
- Bookkeeping is wrong: vacancies count as atoms with full mass (density too
  high), the composition contains a fake element (breaks Geant4/OpenMC and
  other downstream users), and the output is a hack other tools can not
  interpret.

## Agreed design

- **Occupancy is a property of the atom role (`AtomInfo`)**, not of
  `AtomData`. The `AtomData` (possibly a mixture of elements/isotopes) keeps
  meaning "what sits on the site when it is occupied". A vacancy
  pseudo-component in `AtomData`/`@ATOMDB` was considered and rejected (it
  would e.g. make the averaged mass o*m, which is wrong for the dynamics).
- **Vacancies only.** Substitution (several elements sharing sites) already
  works via `@ATOMDB` mixtures, including the disorder term. A partially
  occupied mixed site combines both mechanisms.
- **Format: new NCMAT version (v8) with a per-role `@ATOMSITEOCCUPANCIES`
  section**, listing each atom label once (e.g. `Al 0.8`, fractions like `2/3`
  allowed, default 1, valid range 0<o<=1). Occupancy is thus shared by all
  positions of a role, like `@DYNINFO`. Sites of the same element with
  different occupancies become different roles via `@ATOMDB` aliases (the
  composer already creates such labels). Older format versions reject the
  section, and existing files are unaffected.

## Physics

For a role with occupancy o, coherent scattering length b and Debye-Waller
factor exp(-W):

| Quantity | Treatment |
|---|---|
| Structure factor | Absorb o into the per-role coherent scattering lengths used by `NCFillHKL.cc` (the `csl` arrays), so no per-position weights or special fast paths are needed. |
| Atoms per unit cell | sites*o, used for composition, density and number density. |
| Incoherent elastic | o*sigma_inc*exp(-2W), plus the disorder term 4pi*o(1-o)*b^2*exp(-2W), with b the (mixture-averaged) coherent scattering length of the `AtomData`. |
| Inelastic, absorption | Unchanged per real atom (they follow the occupancy-weighted composition). |

The disorder term assumes randomly placed vacancies (no short-range order),
the standard approximation and what the current workaround gives.

## Why this needs an ABI break

- `AtomInfo` (public, `NCInfoTypes.hh`) needs a new data member, changing its
  memory layout.
- `AtomInfo::numberPerUnitCell()` should return the average number of atoms
  per cell (sites*o, e.g. 3.6 for 4 sites at 90%), so its return type changes
  from `unsigned` to `double`. A new `sitesPerUnitCell()` returns the number
  of positions (4).
- `StructureInfo::n_atoms` (`unsigned`) needs the same decision (sites or
  average atoms).
- The C API (`ncrystal_info_getatominfo` in `cinterface/ncrystal.h`) returns
  an unsigned atom count, and must expose occupancy (with matching
  `_chooks.py` changes, or preferably JSON queries where feasible).

The existing `NCRYSTAL_ALLOW_ABI_BREAKAGE` define can be used to develop the
public-interface parts before 5.0.

## Plan

Each step its own commit(s), with tests:

1. **Done:** refactoring of `NCFillHKL.cc` (shared F^2 calculation, so the
   occupancy only needs to enter the per-role scattering lengths).
2. **No ABI impact:** NCMAT v8 parsing of `@ATOMSITEOCCUPANCIES` into the
   parser's internal data, the scattering length multiplier and the disorder
   term, written against internal interfaces.
3. **Behind `NCRYSTAL_ALLOW_ABI_BREAKAGE`:** the `AtomInfo` member,
   `sitesPerUnitCell()`, and the `numberPerUnitCell()`/`n_atoms` semantics.
   Without the define, NCMAT data using the new section gives a clear "not
   supported in this build" error.
4. **NCrystal 5.0:** C API and Python exposure (`Info.AtomInfo`), and
   replacing the composer/`cif2ncmat` workaround by the new section (writing
   v8 only when needed, reading it back in `from_ncmat`/`from_info`).
5. Documentation: a v8 section in `docs/ncmat_doc.md`, and CHANGELOG entries.

## Validation

- The same structure written with the new section and with the old
  workaround must give identical coherent and incoherent elastic cross
  sections, while density and composition differ (the new ones correct).
- Per-plane F^2 checked against an independent calculation including
  occupancy (reusing `tests/pypath/NCTestUtils/latticecells.py`).
- Limits: o=1 identical to no section; continuity for small o; combination
  with `@ATOMDB` mixtures.
- `cif2ncmat` round trips of partially occupied CIF data (occupancy support
  to be added to `tests/pypath/NCTestUtils/cifgen.py`): composition and
  density must match the CIF formula, without fake elements.
- Errors: o<=0 or o>1, the section in pre-v8 files, `@DYNINFO` fractions
  inconsistent with the weighted composition.

## Related items (also postponed)

- Anisotropic (per-site) displacements also need new NCMAT syntax and could
  share the same format version bump.
- The single crystal d-spacing spread (`delta_d`) is only a stub which throws
  if used; implementing it needs a new broadening model.

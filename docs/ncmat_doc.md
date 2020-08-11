
Description of the NCMAT format
===============================

This document describes the various versions of the NCMAT file format associated
with NCrystal and the extension `.ncmat`. Note documented here is merely the
format itself, and that different versions of NCrystal might choose to do
slightly different things with the same input file - ideally because newer
versions of the code will introduce improved modelling capabilities, but also
because bugs might get fixed (usually) or introduced (hopefully
rarely). Currently the two versions of the format are *NCMAT v1* which is
supported by all NCrystal releases, and *NCMAT v2* which can be read with
NCrystal releases 2.0.0 and beyond.

# The NCMAT v1 format #

The *NCMAT v1* format is the original supported format, in which crystal unit
cells are defined in @CELL, @SPACEGROUP and @ATOMPOSITIONS sections, while a
@DEBYETEMPERATURE section must provide global or per-element Debye
temperatues.

The publication associated with NCrystal release v1.0.0 describes this format in
more detail, including how both elastic and inelastic physics is modelled based
on what is essentially just the layout of the unit cell and a single dynamic
parameter, the Debye temperature:

> X.-X. Cai and T. Kittelmann, NCrystal: A library for thermal neutron
> transport, Computer Physics Communications 246 (2020) 106851,
> https://doi.org/10.1016/j.cpc.2019.07.015

Example file in *NCMAT v1* format, units in the @CELL section are in angstrom
and degrees respectively, while the values in the @DEBYETEMPERATURE section are
in kelvin:

```
NCMAT v1
# Some comments about the file can be
# optionally added here.
@CELL
    lengths 4.913437 4.913437 5.405118
    angles 90. 90. 120.
@SPACEGROUP
    154
@ATOMPOSITIONS
    Si 0.47 0. 0.666666666667
    Si 0. 0.47 0.333333333333
    Si 0.53 0.53 0.
    O 0.4146 0.2678 0.78543
    O 0.2678 0.4146 0.21457
    O 0.7322 0.1468 0.452096666667
    O 0.5854 0.8532 0.881236666667
    O 0.8532 0.5854 0.118763333333
    O 0.1468 0.7322 0.547903333333
@DEBYETEMPERATURE
    Si   515.5240
    O   515.1032
```

The file must always begin with the format designation, `NCMAT v1`, on a single
line, with no space added before the initial `N`. All NCMAT files thus begin
with 5 bytes spelling out `NCMAT` in ASCII. The entire file must be encoded in
ASCII or UTF-8 encodings, with non-ASCII parts of UTF-8 only allowed in
comments. Line endings can be Unix (LF) or DOS (CFLF), but not the ancient
Macintosh endings (CF).

The four sections appearing in the example file above, @CELL, @SPACEGROUP,
@ATOMPOSITIONS, and @DEBYETEMPERATURE must all be present exactly once in all v1
files - with the exception of @SPACEGROUP which is optional. There is no strict
requirement on the order in which they appear. Comments must be on lines
starting with `#` marks and are allowed only after the initial line (with the
format version, "NCMAT v1"), and before the first section. Empty lines are
ignored, as is any extraneous white-space within sections. The section markers
(@CELL, etc.) must be on their own lines. When reading *NCMAT v1* files with
NCrystal v2.0.0 and later, inelastic scattering will be dealt with as for any
*NCMAT v2* file where @DYNINFO sections are absent (see below for what this
means).

## The @CELL section ##

The lattice parameters of the crystal. The format is as shown in the example,
but the order of `length` and `angles` can be reversed if desired.

## The @ATOMPOSITIONS section ##

Atomic unit cell positions in relative lattice coordinates. Each element
appearing must be on separate line, as shown in the example above, and the name
must be Capitalized and use the standard chemical element notation.

## The @SPACEGROUP section ##

This optional section contains the space group number of the crystal symmetry,
which must be mathematically consistent with the unit cell layout provided in
the @CELL and @ATOMPOSITIONS sections. The space group number is given as a
single integer with a value from 1 to 230.

## The @DEBYETEMPERATURE section ##

Either a single number, the global Debye temperature, must be provided, or one
for each of the elements in the material. Each element must be on separate line.

## Historical note about the v1 format parsers ##

The original parsing code associated with the NCMAT v1 format was somewhat
forgiving for certain syntax errors. Thus some clearly broken NCMAT v1 files
could be read with the old code, while the new code will raise an error. It is,
however, the intention that any v1 file which can be read with newer code will
also be readable with older NCrystal code all the way back to NCrystal
v1.0.0. All official v1 files in the data libraries published with NCrystal can
be read with both new and old parsing code.

# The NCMAT v2 format #

The *NCMAT v2* format is similar to the *NCMAT v1* format, but introduces the
possibility of using other models for the modelling of, in particular, inelastic
scattering. Support for it is first introduced in NCrystal release 2.0.0. It
adds two new optional sections, @DYNINFO and @DENSITY, and makes the @CELL,
@ATOMPOSITIONS, and @DEBYETEMPERATURE sections optional in some
circumstances. It also introduces support for the element name "D" (deuterium).

## The @DYNINFO section ##

The main new feature is the addition of optional @DYNINFO sections. There should
exist exactly one for each type of element in the material. Each @DYNINFO
section can have one of five different types, influencing what modelling should
be used for that element: *scatknl*, *vdos*, *vdosdebye*, *freegas* or
*sterile*, all described in the following sections. For the special case of
crystalline materials with both @ATOMPOSITIONS and @DEBYETEMPERATURE sections,
it is allowed to omit all @DYNINFO sections from the file, which then implies
*vdosdebye* type for all elements (which is how *NCMAT v1* files will be read in
NCrystal v2.0.0 and later). Fields always present in any @DYNINFO sections are:

```
@DYNINFO
  element  O
  fraction 2/3
  type     <name-of-type>
```

Here the *element* field indicates the name of the element, the *type* field the
aforementioned type, while the *fraction* field indicates the relative count of
this type of element in the material as a value between 0 and 1. The fractions
in all DYNINFO sections should add up to to 1, and it is possible to express
fractions as rational numbers for higher precision and readability, e.g. `2/3`
(with no spacing allowed, thus no `2 /3`). The order of fields in @DYNINFO
sections does not matter, and unless otherwise noted, the values should appear
on the same line as the field name.

### Dynamic model: Sterile ###

The @DYNINFO type *sterile*, declares a modelling where there is no
scattering by thermal neutrons on that element at all (depending on the context,
there can still be absorption or scattering by particles other than thermal
neutrons).  This is of course does not represent a realistic model of
scattering, but it can nonetheless occasionally be useful if other more detailed
information is not available. For instance it might be used when an element of a
polyatomic material is believed to only contribute significantly to absorption
cross-sections or density calculations. Example section:

```
@DYNINFO
  element  O
  fraction 2/3
  type     sterile
```

### Dynamic model: Free gas ###

The simplest dynamic model with actual scattering of thermal neutrons, indicates
that a basic free-gas model should be chosen for that element. For instance, to
enable free-gas modelling of the oxygen in water molecules (which is a typical
strategy if one has detailed scatter kernels for the hydrogen but no specific
knowledge about the contribution of the oxygen), one would add a section with:

```
@DYNINFO
  element  O
  fraction 2/3
  type     freegas
```

It will depend on the use case or material in question whether a *freegas* or
*sterile* type represents a more suitable fall-back strategy for a given
material.

### Dynamic model: Scattering kernel ###

The most realistic modelling is in the form of 2D scattering kernels, in the
form of either S(alpha,beta) or S(Q,omega) tables. For instance an S(alpha,beta)
table is specified with the *alphagrid*, *betagrid*, and *sab* keywords like:

```
@DYNINFO
  element     Be
  fraction    1
  type        scatknl
  temperature 296.3
  alphagrid   8.9330622e-02 1.3293247e-01 2.2332656e-01 2.9245144e-01
              <...>
              5.9819613e+02 6.2035154e+02 6.6466237e+02 7.0897319e+02
  betagrid    -7.9349821e+01 -7.4390457e+01 -6.9431093e+01 -6.4471730e+01
              <...>
              6.4471730e+01 6.9431093e+01 7.4390457e+01 7.9349821e+01
  sab         3.3877900e-37 1.3833500e-36 4.5036200e-36 1.4220200e-35
              <...>
              1.2642200e-19 1.5147000e-19 1.5901900e-19
```


Of course, the number of entries in the *sab* field must be equal to the product
of the number of entries in the *alphagrid* and *betagrid* fields. To allow for
more efficient processiong, the *alphagrid* and *betagrid* fields can have at
most 65534 entries, and must have at least 5.  In addition to requiring large
amounts of space, scatter kernels are complicated by the fact that they describe
the relevant dynamics at a certain material temperature only. Thus, a
*temperature* field exists in which this value must be specified. Furthermore,
all @DYNINFO sections of type *scatknl* in a given file must contain the same
value of the temperature field, and NCrystal will default the material
temperature to that value as well when loading the input file (and might decide
to raise an error on attempts to override it when loading the file). It is
therefore recommended that *NCMAT* files with @DYNINFO sections of type
*scatknl* should indicate the temperature in the file-name, to avoid user
confusion.

At initialisation time, NCrystal will typically integrate the scattering kernel
for neutrons at a range of predefined energies, in order to be able to later
provide cross-sections and sampling capabilities. By default NCrystal will
analyse the scattering kernel at initialisation time in order to select a
suitable energy grid. This default behaviour can be modified through
specification of the optional *egrid* field inside the @DYNINFO section. If
desired, the exact energy grid can be directly specified by providing at least
10 positive sorted energy values (unit eV):


```
  egrid 1e-5 2e-5 4e-5
        <...>
        4.0 4.5 5.0
```

Alternatively, one can specify just the endpoints and number of points of the
energy grid by providing just 3 values, corresponding respectively to Elow,
Eupper and N of an energy grid running from Elow to Eupper with a total of N
grid points spaced evenly in logarithmic space. An example of this would be
(first two values are in eV, the third must be integral):

```
  egrid 1e-5 5.0 1000
```

Replacing any of the three values with 0 indicates that NCrystal should
determine it automatically. For instance, the following entries would keep the
automatically determined range of the energy grid, but increase the number of
grid points:

```
  egrid 0 0 1000
```

As a final alternative, a single value provided to *egrid* will be interpreted
as Eupper, implying with automatic determination of Elow and number of
points. Thus, Eupper=20eV can be requested with:

```
  egrid 20.0
```

Finally, it is possible to specify `S'(alpha,beta)=S(alpha,beta)*exp(-beta/2)`
values instead of `S(alpha,beta)` values, by using the `sab_scaled` keyword
instead of `sab`. When using `sab_scaled` it is optionally possible to specify
just the table at non-negative values of `beta`, and assume
`S'(alpha,-beta)=S'(alpha,beta)`. In this case, the `betagrid` values must start
at zero and contain only non-negative entries. Irrespective of the input format,
at initialisation time, NCrystal will transform the scattering kernel into a
standard S(alpha,beta) table, for consistent internal treatment.

For convenience and readability, it is allowed to specify parameters of any of
the 1D or 2D "array data" fields on multiple lines. This includes all the
various *XXXgrid* fields as well as *sab* and *sab_scaled*. It is also allowed
to save space by specifying two or more consecutive identical values using the
syntax `<value>r<count>` (i.e.  `0r2000` means the value 0 repeated 2000 times).

Although not supported in the *NCMAT v2*, it might be interesting to note
already that future versions of the NCMAT format are planned to support the
optional specification of scattering kernels as S(Q,omega) tables. This is
planned to simply use `qgrid`, `omegagrid` and `sqw` in place of `alphagrid`,
`betagrid` and `sab`.

### Dynamic model: Vibrational Density of States ###

As a more convenient alternative to providing a scattering kernel directly,
modelling of solid materials might instead be based on a phonon spectrum, or
more specifically a vibrational density of state (vdos) spectrum. This is done
by specifying the density of phonon states as a function of energy. The density
values themselves must be provided in the *vdos_density* field, while the range
of the corresponding energy grid will be given as values in the *vdos_egrid*
field (unit eV). As a convenient and efficient shorthand, when only two points
are provided in the *vdos_egrid* field, these will be assumed to be the first
and last points of the actual energy grid with the remaining points spaced
evenly between those (the total number of points will then be inferred by the
number of points in the *vdos_density* field (which must be at least 5).

The overall normalisation or unit of the values in the *vdos_densiy* field does
not matter, as NCrystal will renormalise them automatically. The density
function will be formed by linear interpolation within the grid. For energies
below the first grid point, the density values will be assumed to scale as the
square of the energy - an approximation which should be suitable for solid
materials. For energies above the last grid point, densities will be assumed to
be 0. To ensure convergence of subsequent calculations in NCrystal, the first
point in the *vdos_egrid* is required to be at least 1e-5 (i.e. 0.01meV).

An example vdos specification might look like:

```
@DYNINFO
  element  V
  fraction 1
  type     vdos
  vdos_egrid   0.013300161256 0.38127128934268
  vdos_density 0.2594369576 0.3556097936 0.4375132376 0.4857022268 0.4915645468
               0.444171744 0.340902404 0.2143757752 0.1140493096 0.0607265688
               <...>
               0.0593700356 0.0693096064 0.0812168444 0.0885391804 0.082003864
               0.0920026688 0.0224714976 0.0040273484
```

For readability, it is allowed to use multiple lines to specify the parameters
of the *vdos_egrid* and *vdos_density* fields. At initialisation time, NCrystal
will use model-dependent assumptions (to be described elsewhere) to expand the
provided density of states into a complete S(alpha,beta) scattering kernel, with
subsequent treatment identical to that of any scattering kernel. For that
reason, an optional *egrid* field can be specified as well, with the same
meaning as when specifying scatter kernels directly.

As a final note, it should be mentioned that NCrystal calculations internally
require the VDOS to be specified on an energy grid with equidistantly spaced
points and which if extended downwards would eventually coincide exactly
with 0. If the VDOS specified in the input does not exactly conform to this
requirement, it will be re-parameterised upon loading (using a procedure which
is unlikely to introduce significant numerical errors).

### Dynamic model: Debye Model for Vibrational Density of States ###

For crystalline materials with a @DEBYETEMPERATURE section, it is allowed to
substitute a realistic VDOS spectrum with an idealised one, in which the density
of states increases quadratically up to a cutoff energy, after which it
becomes 0. In this model due to Debye, phonons states are considered uniformly
distributed in momentum space, up to the cutoff - which is the Debye energy, or
Boltzmann's constant times the Debye temperature. This is represented with a
@DYNINFO type of *vdosdebye*, e.g.:

```
@DYNINFO
  element  V
  fraction 1
  type     vdosdebye
```

### The @DENSITY section ###

In order to support non-crystalline materials for which a vibrational density of
state or scatter kernel is present (like water), it is now optionally allowed to
not specify a unit cell at all via @CELL, @ATOMPOSITIONS and @SPACEGROUP
sections, and instead configure non-crystalline materials exclusively via the
new @DYNINFO sections. However, such materials lack the material densities
otherwise calculated from the unit cell information, which is necessary to
translate per-atomic cross-sections to macroscopic interaction lengths. Thus,
materials without a unit cell, must instead provide their density in the new
@DENSITY section:

```
@DENSITY
  0.4 atoms_per_aa3
```

Other valid units are `kg_per_m3` or `g_per_cm3` (units must be specified just as
written here, no spaces allowed).

## Changes in existing sections ##

### The @SPACEGROUP, @CELL and @ATOMPOSITIONS sections are now optional ###
As mentioned, it is since *NCMAT v2* possible to omit all of the @SPACEGROUP,
@CELL and @ATOMPOSITIONS sections, in order to describe non-crystalline
materials. If at least one of the @CELL and @ATOMPOSITIONS is present, they must
both be present. The @DEBYETEMPERATURE section should never be specified if
there are no @CELL/@ATOMPOSITION sections, and is for now actually also
*required* when there are @CELL/@ATOMPOSITIONS in the file (in later versions
this might become optional for elements for which a @DYNINFO section of vdos
type is present).

### The @ATOMPOSITIONS section supports fractions ###

Positions like 0.666667 can now be more precisely described via fractions like
"2/3" (no spaces are allowed on either side of the fraction).

## Other changes for NCMAT v2 ##

Comments can now appear anywhere in the file, even at the end of other lines
(except the first). Any parts of the line after the `#` should simply be ignored
by the parser.

# Effects of @DYNINFO sections on inelastic NCrystal physics modelling #

Modelling of inelastic scattering with `.ncmat` files will depend on the
NCrystal material configuration parameter, `inelas`. By default, modelling will
be based on the content of the @DYNINFO sections (including implicitly defined
@DYNINFO sections when all are omitted for a crystalline
material). Additionally, it is possible to override the @DYNINFO sections
specified in the files, by setting `inelas=freegas`, `inelas=none` or
`inelas=vdosdebye`. In this case, the indicated type will be applied to all
elements. Of course, `inelas=vdosdebye` can only be used for crystalline
materials with Debye temperatures. Setting `inelas=freegas;elas=0` can be used
to enable a pure freegas model for any input file.

# EMACS Syntax highlighting #

If using EMACS, add the following to your ~/.emacs configuration file to enable
a basic syntax highlighting of .ncmat files:

```
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Major mode for .ncmat files:
(setq ncmat-sections '("CELL" "ATOMPOSITIONS" "SPACEGROUP" "DEBYETEMPERATURE" "DYNINFO" "DENSITY"))
(setq ncmat-fields '("lengths" "angles" "element" "fraction" "type" "temperature"
                     "sab_scaled" "sab" "alphagrid" "betagrid" "egrid" "vdos_egrid" "vdos_density"))
(setq ncmat-highlights
      `(("^NCMAT\s*v[0-9]+[a-z0-9]*" . font-lock-function-name-face)
        (,(concat "^@\\(" (mapconcat 'identity ncmat-sections "\\|") "\\)") . font-lock-keyword-face)
        (,(concat "^\s*\\(" (mapconcat 'identity ncmat-fields "\\|") "\\)") . font-lock-constant-face)
        ))
(setq ncmat-mode-syntax-table
      (let ( (synTable (make-syntax-table)))
        (modify-syntax-entry ?# "<" synTable)
        (modify-syntax-entry ?\n ">" synTable)
        synTable)
      )
(define-derived-mode ncmat-mode fundamental-mode "NCMAT"
  "major mode for editing .ncmat files."
  (setq font-lock-defaults '(ncmat-highlights));keywords
  (set-syntax-table ncmat-mode-syntax-table);comments
  )
(add-to-list 'auto-mode-alist '("\\.ncmat$" . ncmat-mode))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
```

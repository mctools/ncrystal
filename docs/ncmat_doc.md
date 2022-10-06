
Description of the NCMAT format
===============================

This document describes the various versions of the NCMAT data format associated
with NCrystal and the file extension `.ncmat`. Note that documented here is
merely the format itself, and that different versions of NCrystal might choose
to do slightly different things with the same input file - ideally because newer
versions of the code will introduce improved modelling capabilities, but also
because bugs might get fixed (usually) or introduced (hopefully rarely).
Currently the versions of the format are *NCMAT v1* which is supported by all
NCrystal releases, *NCMAT v2* which can be read with NCrystal releases 2.0.0 and
beyond, *NCMAT v3* which can be read with NCrystal releases 2.1.0 and beyond,
*NCMAT v4* which can be read with NCrystal releases 2.6.0 and beyond, *NCMAT v5*
which can be read with NCrystal releases v2.7.0 and beyond, *NCMAT v6* which can
be read with NCrystal releases v3.0.0 and beyond, and *NCMAT v7* which can be
read with NCrystal releases v3.2.0 and beyond.

# The NCMAT v1 format #

The *NCMAT v1* format is the original supported format, in which crystal unit
cells are defined in @CELL, @SPACEGROUP, and @ATOMPOSITIONS sections, while a
@DEBYETEMPERATURE section must provide global or per-element Debye
temperatures.

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
comments. Line endings can be Unix (LF) or DOS (CRLF), but not the ancient
Macintosh endings (CR).

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
appearing must be on a separate line, as shown in the example above, and the name
must be Capitalised and use the standard chemical element notation.

## The @SPACEGROUP section ##

This optional (but highly recommended) section contains the space group number
of the crystal symmetry, which must be mathematically consistent with the unit
cell layout provided in the @CELL and @ATOMPOSITIONS sections. The space group
number is given as a single integer with a value from 1 to 230.

## The @DEBYETEMPERATURE section ##

Either a single number, the global Debye temperature, must be provided, or one
for each of the elements in the material. Each element must be on a separate
line. Note that the concept "global Debye temperature" in this sense, just
implies that the same value will be applied to all elements.

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

The @DYNINFO type *sterile*, declares a modelling where there is no scattering
by thermal neutrons on that element at all (depending on the context, there can
still be absorption or scattering by particles other than thermal neutrons).
This of course does not represent a realistic model of scattering, but it can
nonetheless occasionally be useful if other more detailed information is not
available. For instance it might be used when an element of a poly-atomic
material is believed to only contribute significantly to absorption
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
strategy if one has detailed scattering kernels for the hydrogen but no specific
knowledge about the contribution of the oxygen), one would add a section with:

```
@DYNINFO
  element  O
  fraction 2/3
  type     freegas
```

It will depend on the use case or material in question whether a *freegas* or
*sterile* type represents a more suitable fall-back strategy for a given
element.

### Dynamic model: Scattering kernel ###

The most realistic modelling is in the form of 2D scattering kernels, in the
form of either S(alpha,beta) or S(Q,omega) tables. For instance, an S(alpha,beta)
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
more efficient processing, the *alphagrid* and *betagrid* fields can have at
most 65534 entries, and must have at least 5.  In addition to requiring large
amounts of space, scattering kernels are complicated by the fact that they describe
the relevant dynamics at a certain material temperature only. Thus, a
*temperature* field exists in which this value must be specified. Furthermore,
all @DYNINFO sections of type *scatknl* in a given file must contain the same
value of the temperature field, and NCrystal will default the material
temperature to that value as well when loading the input file (and might decide
to raise an error on attempts to override it when loading the file). It is
therefore recommended that *NCMAT* files with @DYNINFO sections of type
*scatknl* should indicate the temperature in the file name, to avoid user
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
as Eupper, implying automatic determination of Elow and number of points. Thus,
Eupper=20eV can be requested with:

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

Although not supported in the *NCMAT v2* format, it might be interesting to note
already that future versions of the NCMAT format are planned to support the
optional specification of scattering kernels as S(Q,omega) tables. This is
planned to simply use `qgrid`, `omegagrid` and `sqw` in place of `alphagrid`,
`betagrid` and `sab`.

### Dynamic model: Vibrational Density of States ###

As a more convenient alternative to providing a scattering kernel directly,
modelling of solid materials might instead be based on a phonon spectrum, or
more specifically a vibrational density of state (vdos) spectrum. This is done
by specifying the density of phonon states as a function of energy (which is
hbar times the phonon frequency). The density values themselves must be provided
in the *vdos_density* field, while the corresponding energy grid will be given
as values in the *vdos_egrid* field (unit eV). As a convenient and efficient
shorthand, when only two points are provided in the *vdos_egrid* field, these
will be assumed to be the first and last points of the actual energy grid with
the remaining points spaced evenly between those (the total number of points
will then be inferred by the number of points in the *vdos_density* field -
which must always be at least 5).

The overall normalisation or unit of the values in the *vdos_densiy* field does
not matter, as NCrystal will renormalise them automatically. The actual density
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
meaning as when specifying scattering kernels directly.

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
becomes 0. In this model due to Debye, phonon states are considered uniformly
distributed in momentum space, up to the cutoff - which is the Debye energy, or
Boltzmann's constant times the Debye temperature. This is represented with a
@DYNINFO type of *vdosdebye*, e.g.:

```
@DYNINFO
  element  V
  fraction 1
  type     vdosdebye
```

Note that in *NCMAT v5* (described below), a new optional keyword *debye_temp*
is added for this section, allowing it to be used in files without a separate
@DEBYETEMPERATURE section.

### The @DENSITY section ###

In order to support non-crystalline materials for which a vibrational density of
state or scattering kernel is present (like water), it is now optionally allowed to
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
both be present. The @DEBYETEMPERATURE section should likewise never be specified if
there are no @CELL/@ATOMPOSITION sections in the file, and is also required if there
are @CELL/@ATOMPOSITIONS in the file (as will be mentioned below this changes in
*NCMAT v4* where the @DEBYETEMPERATURE section becomes optional in some cases).

### The @ATOMPOSITIONS section supports fractions ###

Positions like 0.666667 can now be more precisely described via fractions like
"2/3" (no spaces are allowed on either side of the division symbol).

## Other changes for NCMAT v2 ##

Comments can now appear anywhere in the file, even at the end of other lines
(except the first). Any parts of the line after the `#` should simply be ignored
by the parser.

## Effects of @DYNINFO sections on inelastic NCrystal physics modelling #

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

# The NCMAT v3 format #

The *NCMAT v3* format is similar to the *NCMAT v2* format, but introduces two
new features. The first new feature is the introduction of more flexibility in
how atoms can be defined, which both concerns how atoms can be specified in
existing sections and also adds a new optional @ATOMDB section.  Although not
pertaining to the NCMAT format as such, it should nonetheless be noted here that
starting with NCrystal v2.1.0 it is also possible to modify atomic definitions
for any ncmat file (even those in NCMAT v1 or NCMAT v2 formats), through the
usage of the atomdb configuration parameter, using a syntax similar to the one
used in @ATOMDB sections (which will be described below). For details about the
atomdb configuration parameter see elsewhere (for instance at
https://github.com/mctools/ncrystal/wiki/CfgRefDoc).

The other new feature in *NCMAT v3* is that it is now possible to add custom
sections to NCMAT files. The content of these custom sections will be loaded
into the internal material Info objects, but the core NCrystal code will
otherwise ignore them. This is primarily intended as a way to facilitate the
development of new experimental physics models which might need additional
data. Of course, if a such an experimental model matures (and is deemed to be of
interest to a wider community), and makes its way into the core NCrystal code,
the NCMAT format should evolve further in order to include proper well defined
sections for the new data.

## Changes in existing sections ##

In addition to specifying atoms via their element names (H, He, Li, ...) in
@DYNINFO, @DEBYETEMPERATURE, and @ATOMPOSITIONS sections, it is now also
possible to specify atoms as specific isotopes by postfixing the nucleon number
to the element name (H1, H2, He3, Li6, Gd157, ...). Two special aliases are
provided: The name D is shorthand for H2 and T is shorthand for H3. This is not
only handy, it also provides backwards compatibility with NCMAT content
originally written for the NCMAT v2 format, in which D was made available to
describe deuterium.

For more complicated atomic setups, one must use an @ATOMDB section (or `atomdb`
configuration parameter as noted above), which might redefine the meaning of
element names throughout the files (e.g. define particular isotopic abundances
in "B", or define "Al" as actually being a mixture of 99% Al and 1% Cr), or
provide/override physics constants associated with a given isotope or
element. It can also define any of the generic markers: X, X1, X2, ..., X99. If
such are defined, they may be used in place of element names throughout the
existing sections.

## The @ATOMDB section ##

The simplest usage of the @ATOMDB section is to provide physics constants
associated with certain natural elements or isotopes. The four constants which
must be provided after the atom symbol are respectively: mass, coherent
scattering length, incoherent cross section, and absorption cross
section. Masses and cross sections must be non-negative. For readability, units
must be indicated by postfixes: Daltons (postfix `u`), femtometres (postfix
`fm`), and barns (postfix `b`). There should be no space between the number and
the unit. An example of this would be:

```
@ATOMDB
  Si    28.09u 4.1491fm 0.004b 0.171b
  Si29  28.97649466525u 4.7fm 0.001b  0.101b
  Rn222 222.018u -1.23fm 4.56b 7.89b
```

The first line provides data for silicon in natural abundance, the second line
provides data for the silicon isotope with 29 nucleons, and the last line
provides (nonsense) data for the radon isotope with 222 nucleons. More
accurately, the first two lines *override* the values for Si and Si29 which
would otherwise be taken from the internal database shipped with NCrystal. On
the contrary, Rn222 is an example of an isotope which is not available in the
internal database (at least not in NCrystal release v2.1.0), so NCMAT files
referring to that isotope *must* provide values for it in the @ATOMDB
section. It is possible (perhaps for validation purposes) to disable the inbuilt
database, which is done by placing the keyword `nodefaults` on the first line of
the @ATOMDB section:

```
@ATOMDB
  nodefaults
  Si29 28.97649466525u 4.7fm 0.001b  0.101b
```

In this example, only `Si29` will be available for usage in the file. Concerning
isotope labels (as was already mentioned) the labels, D and T are simply taken
to be aliases for the isotope labels H2 and H3 respectively, wherever they are
encountered.

In addition to specifying physics constants associated with elements and
isotopes, the @ATOMDB section also allows the definition of mixtures of elements
and isotopes. This allows both flexible isotopic compositions, but can also be
used to model impurities where for instance a contaminant replaces a primary
element on certain positions in a crystal lattice. Here is an example which
specifies a particular enrichment of boron (95% B10):

```
@ATOMDB
  B is 0.95 B10 0.05 B11
```

Obviously, the specified fractions must add up to 1. Next is an example which
introduces trace amounts of chromium into the positions otherwise occupied by
aluminium atoms (assuming the other sections of the NCMAT file define a crystal
involving aluminium atoms):

```
@ATOMDB
  Al is 0.99 Al 0.01 Cr
```

Lines in the @ATOMDB section are evaluated in order, and in essence can be
thought of as defining or redefining the variable which is named at the
beginning of the line. This can be used to build up more complicated
scenarios. If for instance one would like to model boron carbide (B4C) in which
the boron is enriched to 90% boron-10, and where 1 permille of boron atoms are
replaced by carbon atoms, one can write (assuming the other sections of the
NCMAT file defines a B4C lattice in the obvious manner, putting a B or a C atom
at each lattice position):

```
@ATOMDB
  B is 0.9 B10 0.1 B11
  B is 0.999 B 0.001 C
```

As the lines are evaluated in order, a mathematically equivalent but less handy
version would be:

```
@ATOMDB
  B is 0.8991 B10 0.0999 B11 0.001 C
```

To avoid confusion, it is not allowed to use isotope labels to define mixtures:

```
@ATOMDB
  B10 is 0.999 B10 0.001 B11 # This is invalid syntax!
```

On the other hand, a number of "generic labels" are available, in case one would
like to avoid using a standard chemical symbol to denote a mixture. Exactly 100
of such generic markers are available and they are: X, X1, X2, X3, ..., X98, and
X99. Apart from Xe (xenon), no standard element starts with X, which should make
it easy to identify such custom mixtures in a given NCMAT file. It makes no
difference which of the markers X, X1, ..., X99 are used, and the data loaded
from a given NCMAT file will contain no indication of the name used in the input.

With generic labels, one can thus define:

```
@ATOMDB
  X is 0.666666666666666666667 H 0.333333333333333333333 O
```

and then use X in place of an atom name in the other sections of the NCMAT
file. Again, to avoid confusion, generic labels (X, X1, ..., X99) can only be
used to denote mixtures, and it is not allowed to assign physical properties
directly to such a label:

```
@ATOMDB
  X is 12.5u 0.5fm 3b 0.6b # This is invalid syntax!
```

A mixture consisting of only a single component can be better thought of as an
alias:

```
@ATOMDB
  X5 is 1.0 Al
```

And in fact, for readability, it is allowed to leave out the redundant 1.0 from
such a single-component definition, making it more clear that it is essentially
just an alias definition:

```
@ATOMDB
  X5 is Al
```

Finally, it should be noted that some lines in the @ATOMDB section might have no
effect. For instance:

```
#(assume that the rest of the file refers to Al and O only)
@ATOMDB
  B is 0.9 B10 0.1 B11  #<---- has no effect, B is not used in file
  Al is 0.99 Al 0.01 Cr #<---- has no effect, next line redefines Al.
  Al 26.98u 3.449fm 0.0082b 0.231b
```

It is syntactically OK to include lines which have no effect, but it is possible
that a future version of NCrystal could emit warnings when encountering such
lines.

## Custom sections (@CUSTOM_*) ##

Any number of custom sections can be added. They must have a name of the form
@CUSTOM_<name> where `<name>` must consist of 1 or more upper case characters,
A-Z. Multiple custom sections with the same name are allowed. The contents will
be parsed like in all other sections, meaning that comments, superfluous
white-space and empty lines will all be stripped away. The remaining lines and
their contents will be available after loading as lists of lists of strings.
Example of how a custom section might look:

```
@CUSTOM_MYNEWPHYSICSDATA
  #This section is needed for my great new physics model
  6.789 Xe 0.95 0x0003450 true
```

Note that NCrystal might produce warning messages when loading NCMAT files with
custom sections, serving as a reminder that such files are likely only intended
for a particular development setup and might not be sensible to share with a
wider community.

# The NCMAT v4 format #

The *NCMAT v4* format is similar to the *NCMAT v3* format, but with two changes.
Firstly, Debye temperatures are no longer required to be specified in crystals
which have actual phonon VDOS spectra available, and global Debye temperatures
are no longer allowed. Secondly, unit cell definitions in the @CELL section
can now optionally use a syntax which is slightly more terse and avoids needless
repetition.

## Changes for the @DEBYETEMPERATURE section ##

Based on the fact that it is possible to estimate mean-squared-displacements
(and corresponding Debye temperatures) from a phonon VDOS curve, elements which
have such curves available in @DYNINFO sections of type vdos are no longer required
to have entries in the @DEBYETEMPERATURE section. And if this is the case for all
elements, the @DEBYETEMPERATURE section can be left out entirely.

Note that it is still *allowed* to specify @DEBYETEMPERATURE entries for such
elements. If this is done, the values specified in the @DEBYETEMPERATURE section
will take precedence, but NCrystal might emit a warning message to make sure this
was intended. Note that the NCrystal developers would in general advice against
specifying such Debye temperature values for elements with VDOS curves, as the
values estimated from VDOS curves seems to provide more robust and trustworthy
physics performance.

Finally, global Debye temperatures are no longer allowed in *NCMAT v4*.

## Changes for the @CELL section ##

To reduce repetion, *NCMAT v4* introduces options for specifying data in the
@CELL section more tersely. Firstly, cubic crystals can now be specified
by just providing one lattice length. Thus, writing:

```
@CELL
  cubic 4.04958
```

has the same effect as writing:

```
@CELL
 lengths 4.04958 4.04958 4.04958
 angles 90 90 90
```

When using the `cubic` keyword in the @CELL section, the space group number
in the @SPACEGROUP section (if present) must indicate a cubic material as
well (i.e. be in the range 195..230).

The other new option that is introduced is that the second or third parameters
after the `length` keyword, can be replaced by two exclamation marks, `!!`,
which has the effect of repeating the previous value on the line (for those
curious, `!!` was chosen for this since in many Unix shells `!!` is a command
which means "repeat previous command"). Thus,

```
@CELL
  lengths 2.2866 !! 3.5833
  angles 90 90 120
```

Has the same meaning as:

```
@CELL
  lengths 2.2866 2.2866 3.5833
  angles 90 90 120
```

Although this syntax is not quite as terse as is what possible for cubic
materials with the `cubic` keyword, it still means that e.g. hexagonal
crystals can now be defined without repetition of the same value.

# The NCMAT v5 format #

The *NCMAT v5* format is similar to the *NCMAT v4* format, but with two changes
meant to facilitate modelling of amorphous solids (i.e. non-crystalline solids).
Firstly, Debye temperatures can now be specified directly in @DYNINFO sections
of *vdosdebye* type, meaning that such sections can now also be used in
non-crystalline materials where @DEBYETEMPERATURE sections are disallowed.
Secondly, a new optional @STATEOFMATTER section can be used to explicitly
indicate the state of matter if needed.

## Changes for the @DYNINFO section ##

In files without a @DEBYETEMPERATURE section, it is now allowed to specify the
Debye temperatures directly inside @DYNINFO sections of type *vdosdebye*,
providing them in units of kelvin via the new keyword `debye_temp`, as for
instance:

```
@DYNINFO
  element  V
  fraction 1
  type     vdosdebye
  debye_temp 350
```

Note that files for crystalline materials (those having @CELL and @ATOMPOSITIONS
sections) which do not have a @DEBYETEMPERATURE section, are instead required to
have all @DYNINFO sections being of type *vdos* or *vdosdebye*. This restriction
is needed to ensure adequate information required for Debye Waller factor
estimation in all crystalline materials.

## The @STATEOFMATTER section ##

This new section is optional and is intended to provide a way to indicate the
state of matter in cases where this can not be automatically inferred. Possible
states of matter are at present just *solid*, *liquid*, and *gas*, but that list
might be extended or detailed in the future. The syntax is very simple as
illustrated by the following example:

```
@STATEOFMATTER
  liquid
```

Crystalline materials (those having @CELL and @ATOMPOSITIONS sections) or
materials where one or more @DYNINFO sections are of *vdos* or *vdosdebye* type,
are automatically classified as *solid*, so for such materials the
@STATEOFMATTER section is optional. It can nonetheless be provided, but is then
only allowed to specify the state as *solid*. Future NCMAT versions might change
this requirement (for instance to allow VDOS usage for liquids).

## Effect of state of matter on NCrystal physics modelling #

At present (NCrystal v2.7.0), the effect of the state of matter type is mostly
cosmetic, with one important exception: Non-crystalline solids will get elastic
scattering added based on the assumption that atoms vibrate around fixed
positions in the material. Technically, this will use atomic mean squared
displacement information inferred from the @DYNINFO sections (for now just those
of *vdos* or *vdosdebye* type, but it could possibly at some point also be
inferred from other types such as *scatknl*). Furthermore, coherent-elastic
scattering will for now simply be accounted for via the incoherent
approximation.

# The NCMAT v6 format #

The *NCMAT v6* format is similar to the *NCMAT v5* format, but adds a new
optional @OTHERPHASES section which can be used to add one or more secondary
material phases along with the primary one defined in the file itself. This
support for definition of multiphase materials in an NCMAT file should be
considered a cautious first step, and it is expected that future versions of the
NCMAT format will expand upon this with a more evolved syntax concerning phase
composition and phase-contrast physics (SANS).

## The @OTHERPHASES section ##

The @OTHERPHASES section allows the definition of one or more secondary phases,
by specifying phase volume fractions and cfg strings on one line per secondary
phase. For instance the following section can be used to add 5% (by volume) of
aluminium and 20% of copper with an increased dcutoff. The remaining 75% of the
volume will be composed of the primary phase, defined in the rest of the file in
the usual manner:

```
@OTHERPHASES
  0.05 Al_sg225.ncmat
  0.2  Cu_sg225.ncmat;dcutoff=0.5
```

The cfg-strings provided after the fractions in the @OTHERPHASES section will be
used to set up the phases by directly passing them onto the the plugin and
factory infrastructure available in the NCrystal installation, relying on the
global functions such as createInfo(..) and createScatter(..). This means that
the cfg-strings can even refer to other file-types than `.ncmat` files or employ
a multi-phase syntax like (e.g. `phases<...>`).

Once loaded, the primary phase defined in the NCMAT file itself will become the
first phase in a multiphase Info object, and the secondary phases from the
@OTHERPHASES section will follow it in the order specified.

It should be noted that, due to the constraints of the NCMAT format, there are a
few restraints concerning whitespace and non-ASCII filenames in the cfg-strings
specified in the @OTHERPHASES section. In practice these restraints are likely
to be only rarely relevant, but we list them here in any case for
completeness. Firstly, newline or `#` characters can not be used at all in the
cfg-strings for obvious reasons. Secondly, whitespace is "normalised"
during parsing. In essence this means that all groups of internal whitespace in
the cfg-strings in the @OTHERPHASES section will be parsed as a single space
character. Finally, as non-ASCII characters are not allowed in NCMAT data (the
UTF-8 encoding is only allowed in comments), it is not possible to use non-ASCII
characters in data names.

# The NCMAT v7 format #

The *NCMAT v7* format is similar to the *NCMAT v6* format, but adds a new
optional @TEMPERATURE section.

## The @TEMPERATURE section ##

The @TEMPERATURE section can be used to modify the default temperature value, or lock
the temperature to a specific value. For instance, changing the default
temperature to 400K is done by adding:

```
@TEMPERATURE
  default 400.0
```

If the default keyword is omitted, the value is instead locked in the sense that
any attempts at using the cfg-level variable `temp`  to modify the temperature
will result in an error:

```
@TEMPERATURE
  400.0
```

In either case, if a @DYNINFO section in the same file also specifies a
temperature value, the value must be the same as specified in the @TEMPERATURE
section (in this case the @TEMPERATURE section has no effect and is at most
syntactic sugar). Valid temperature values must be above 0K and at most 1e6K.

# EMACS Syntax highlighting #

If using EMACS, add the following to your ~/.emacs configuration file to enable
a basic syntax highlighting of .ncmat files:

```
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Major mode for .ncmat files:
(setq ncmat-sections '("CELL" "ATOMPOSITIONS" "SPACEGROUP" "DEBYETEMPERATURE" "DYNINFO" "TEMPERATURE"
                       "DENSITY" "ATOMDB" "STATEOFMATTER" "OTHERPHASES" "CUSTOM_[A-Z]+"))
(setq ncmat-fields '("lengths" "angles" "cubic" "element" "fraction" "type" "debye_temp" "temperature"
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

# Browsing and searching NCrystal data

NCrystal ships with a large library of materials (the "standard library" of
`.ncmat` files), and can also use your own files, in-memory data, and data
created on-demand from cfg-strings like `solid::B4C/2.52gcm3` or
`gasmix::air`. This document describes two tools for finding your way in all
of that data:

* The `ncrystal browse` command (also available as `ncrystal_browse`).
* The `NCrystal.browse` Python module, which provides the same features in
  Python scripts and Jupyter notebooks.

They can also show the data in NCrystal's database of isotopes and natural
elements (see the end of this document).

Both let you list the available data, search in names and descriptions, and
even select materials based on their physics, for instance "all crystals
containing boron with an absorption cross section above 100 barn".

In the following, the term *entries* is used for all kinds of data NCrystal
can find by name, since not all of them are actual files.

## Listing available data

Running the command without arguments lists everything, grouped by the source
providing it (standard library, current directory, in-memory data, ...), with
a short description for each NCMAT entry, taken from the comments at the top
of the data:

```
$ ncrystal browse
==> 134 entries from "stdlib" (/path/to/data, priority=120):
    AcrylicGlass_C5O2H8.ncmat               Polymethyl-methacrylate (PMMA).
    AgBr_sg225_SilverBromide.ncmat          Silver Bromide (AgBr, cubic, SG...
    Ag_sg225.ncmat                          Silver (Ag)
    ...
==> 9 entries from "gasmix" (examples of on-demand gas mixtures, ...):
    gasmix::0.72xCO2+0.28xAr/massfractions/1.5atm/250K
    ...
```

The location shown for the standard library depends on your installation
(the data is often embedded directly in the NCrystal library). Long
descriptions are shortened to fit the terminal; use `--no-truncate` to see
them in full.

Use `-f` (`--factory`) to only see data from one source, e.g. `-f stdlib` for
the standard library, or `-f relpath` for files in the current directory.

## Finding data by name

Give one or more patterns to select entries by name. By default a pattern is a
case-insensitive part of the name:

```
$ ncrystal browse -f stdlib Al
==> 21 entries from "stdlib" (/path/to/data, priority=120):
    Al2O3_sg167_Corundum.ncmat              Corundum (alpha alumina, alpha-Al...
    Al4C3_sg166_AluminiumCarbide.ncmat      Aluminium carbide (Al4C3, trigona...
    AlN_sg186_AluminumNitride.ncmat         Aluminum nitride (AlN, hexagonal,...
    Al_sg225.ncmat                          Aluminium (Al, fcc, cubic, SG-225...
    CaF2_sg225_CalciumFlouride.ncmat        Calcium Flouride (CaF2, cubic, SG...
    ...
```

Note how `Al` also matched `CaF2_..._CalciumFlouride.ncmat`. To be more
precise, use glob patterns with `*` and `?` (remember the quotes, so your
shell does not expand them):

```
$ ncrystal browse "Al*"          # names starting with Al
$ ncrystal browse "*_sg225*"     # all materials with space group 225
$ ncrystal browse "stdlib::Be*"  # "::" matches against the full name
```

If nothing matches, you might get a suggestion:

```
$ ncrystal browse Diamnd
No matching data found.
Did you mean: C_sg227_Diamond.ncmat?
```

For more complicated cases, `-E` (`--regex`) turns the patterns into
(case-insensitive) Python regular expressions, e.g. `-E "^b.*_sg1[0-9]{2}"`.

## Searching in the descriptions

Most NCMAT files start with comments describing the material, where it comes
from, and how it was made. Use `-s` (`--search`) to find entries containing
certain words in their name or these comments. When several words are given,
all of them must be present, and the matching lines are shown:

```
$ ncrystal browse -f stdlib -s togo -s lithium
==> 4 entries from "stdlib" (/path/to/data, priority=120):
    Li2O_sg225_LithiumOxide.ncmat    Lithium oxide (Li2O, cubic, SG-225 / Fm-3m)
        | Lithium oxide (Li2O, cubic, SG-225 / Fm-3m)
        | Atsushi Togo and Isao Tanaka, Scr. Mater., 108, 1-5 (2015)
    Li3N_sg191_LithiumNitride.ncmat  Lithium nitride (Li3N, hexagonal, SG-191...
    ...
```

In a terminal, the hits are highlighted with colors, just like with `grep`
(control this with `--color=always|never|auto`; the usual `NO_COLOR` and
`GREP_COLORS` environment variables are respected).

Words are matched literally, so to search for "boron" *or* "b4c", use a
regular expression with `-E`:

```
$ ncrystal browse -f stdlib -E -s "boron|b4c"
```

To read the full comments, add `-c` (`--comments`), and to see the complete
content of an entry, use `-x` (`--extract`):

```
$ ncrystal browse -c Al_sg225.ncmat
$ ncrystal browse -x Al_sg225.ncmat
```

## Selecting materials by their physics

With `-w` (`--where`) you can select entries based on the physics of the
materials, using Python expressions. This requires loading the materials,
which is done in parallel and typically takes less than a second for the
whole standard library.

```
$ ncrystal browse -w "'B' in elements"
$ ncrystal browse -w "state == 'liquid'"            # all liquids
$ ncrystal browse -w "absxs > 100 and crystal"
$ ncrystal browse -w "'scatknl' in dyninfo"        # has a full scattering kernel
$ ncrystal browse -w "crystalsystem == 'hexagonal'" -w "braggthreshold > 8"
$ ncrystal browse -w "max(debyetemps.values()) > 1000"
```

Several `-w` options must all be true. The most useful properties are:

| Property                 | Meaning                                          |
|--------------------------|--------------------------------------------------|
| `elements`, `atoms`      | sets of element names and atom labels, e.g. `{'H','O'}` and `{'D','O'}` |
| `formula`                | chemical formula, e.g. `'Al2O3'`                 |
| `absxs`, `scatxs`        | absorption (at 2200m/s) and free scattering cross sections per atom [barn] |
| `cohxs`, `incohxs`       | bound coherent and incoherent scattering cross sections per atom [barn] |
| `density`, `numdens`     | density [g/cm3] and number density [atoms/Aa3]   |
| `temp`, `state`          | temperature [K], and `'solid'`, `'liquid'`, or `'gas'` |
| `crystal`, `crystalsystem`, `sg` | crystallinity, crystal system (e.g. `'cubic'`) and space group |
| `a`, `b`, `c`, `alpha`, `beta`, `gamma`, `volume` | unit cell parameters |
| `braggthreshold`         | longest wavelength with Bragg diffraction [Aa]   |
| `debyetemps`, `msds`     | per-atom Debye temperatures and mean-squared-displacements (dicts) |
| `dyninfo`                | the kinds of dynamics present: `'vdos'`, `'vdosdebye'`, `'scatknl'`, `'freegas'`, `'sterile'` |

Run `ncrystal browse --help` for the complete list. Expressions which fail
because a value is not available (for instance `sg > 200` for a liquid) are
simply considered false. To see all the property values of some entries, use
`--props`, or `--info` (see below).

## Tables, sorting, and exporting

Instead of the usual listing, you can get a table with the properties of your
choice with `--columns`, and sort with `--sort` (and `--reverse`):

```
$ ncrystal browse -f stdlib -w "absxs > 50" --sort absxs --reverse --columns formula,density
NAME                               FORMULA  DENSITY    ABSXS
B4C_sg166_BoronCarbide.ncmat           CB4  2.48754  613.601
Dy2O3_sg206_DysprosiumOxide.ncmat    Dy2O3  8.25017    397.6
BO3H3_sg2_BoricAcid.ncmat            BH3O3  2.28803  109.714
Au_sg225.ncmat                          Au  19.2877    98.65
Ag_sg225.ncmat                          Ag  10.5013     63.3
Li3N_sg191_LithiumNitride.ncmat       Li3N  1.29522    53.35

$ ncrystal browse -f stdlib -w "crystalsystem=='hexagonal' and braggthreshold > 8" \
                  --columns formula,sg,a,c,braggthreshold
NAME                                FORMULA   SG       A       C  BRAGGTHRESHOLD
GaSe_sg194_GalliumSelenide.ncmat       GaSe  194    3.75   15.92           15.92
LaBr3_sg176_LanthanumBromide.ncmat    Br3La  176  7.9648  4.5119         13.7954
SiO2-beta_sg180_BetaQuartz.ncmat       O2Si  180   5.013    5.47         8.68277
```

Tables are never cut to fit the screen (except for a `description` column).
For use in other programs, the same tables can be written as CSV (with full
numerical precision), JSON, or HTML, and all available information about the
selected entries can be exported as JSON:

```
$ ncrystal browse "Be_*" --columns formula,sg,density --csv > beryllium.csv
$ ncrystal browse "Be_*" --columns formula,sg,density --json
$ ncrystal browse "Be_*" --columns formula,sg,density --html > beryllium.html
$ ncrystal browse "Be_*" --json
```

Finally, a few options are handy in scripts:

```
$ ncrystal browse -w "'Gd' in elements" --count   # just the number of entries
$ ncrystal browse "*_sg225*" --names              # one name per line
$ ncrystal browse -f stdlib Al_sg225 --path       # location of on-disk files
```

## All details about an entry

The `--info` option shows everything known about the selected entries,
including the physics properties, the full header comments, and examples of
how to use them:

```
$ ncrystal browse --info Al_sg225.ncmat
==> Al_sg225.ncmat
    Description   : Aluminium (Al, fcc, cubic, SG-225 / Fm-3m)
    Full key      : stdlib::Al_sg225.ncmat
    Source        : /path/to/data (factory "stdlib", priority 120)
    Data type     : ncmat
    On-disk path  : /path/to/data/Al_sg225.ncmat
  Physics properties:
    elements={'Al'}  atoms={'Al'}  nelements=1  formula='Al'  mass=26.9815
    absxs=0.231  scatxs=1.39667  cohxs=1.49485  incohxs=0.0082  density=2.69865
    ...
  Header comments:
    # Aluminium (Al, fcc, cubic, SG-225 / Fm-3m)
    ...
  Usage examples:
    # Plot cross sections:
    nctool "Al_sg225.ncmat"
    ...
```

## Hidden entries

The same name can be provided by several sources. For instance, a file
`Al_sg225.ncmat` in your current directory is used instead of the one in the
standard library, since files in the current directory have a higher
priority. The listing then marks the other entry with `(hidden)` and shows its
full name, which can always be used to select it explicitly, e.g. in a
cfg-string like `stdlib::Al_sg225.ncmat;temp=20K`.

## Using the Python API

The `NCrystal.browse` module provides all of the above in Python. The main
tool is the `DataBrowser` class:

```python
import NCrystal.browse as nb

b = nb.DataBrowser( factory = 'stdlib' )
b.match('Al*').dump()                      # same listing as the command
sel = b.where('absxs > 50').sorted('absxs', reverse = True)
print( sel.table('formula,absxs,density') )
print( len(sel), sel.names() )
```

Selections (`match`, `search`, `where`, `filter`, `from_factory`, `sorted`,
and slicing like `sel[:5]`) return new `DataBrowser` objects and never modify
the original, so they are easy to combine and reuse:

```python
crystals = b.where('crystal')
hydrides = crystals.search('hydride')
light = crystals.where( lambda p : p.mass < 12 )   # functions work too
by_density = light.sorted('density')
```

Iterating over a browser gives `DataEntry` objects, with the name,
description, comments, and more. Their `props` attribute holds the physics
properties, which are loaded automatically when needed:

```python
for e in b.search('togo').match('Li'):
    print( e.display_name, e.description )
    print( '   density:', e.props.density, 'elements:', sorted(e.props.elements) )

al = b.match('Al_sg225')[0]
mat = al.load( ';temp=20K' )               # an NCrystal.LoadedMaterial
```

The output methods return strings, apart from the `dump*` methods which print:

```python
sel.dump( comments = True )                # also dump( props = True )
sel.dump_info()                            # like --info
csv_text = sel.to_csv('formula,sg,density')
html = sel.to_html('formula,crystalsystem,description')
data = sel.to_dicts()                      # everything, as Python objects
```

In a Jupyter notebook, a table can be shown with:

```python
from IPython.display import HTML
HTML( sel.to_html('formula,sg,density,description') )
```

If you just want the raw data, `nb.query_data()` returns everything as a list
of dictionaries (or as a JSON string with `as_json=True`).

Materials are only loaded when physics properties are needed, and then in
parallel. To prevent loading altogether, create the browser with
`physics=False`. The number of threads used for loading can be set with the
`nthreads` parameter. It is only used temporarily, and not at all if you
have already configured NCrystal's factory threads, for instance with
`NCrystal.enableFactoryThreads(..)` or the `NCRYSTAL_FACTORY_THREADS`
environment variable.

## The database of isotopes and elements

NCrystal contains a database with the masses, scattering lengths and cross
sections of natural elements and many isotopes. Add `--atomdb` to browse it
instead of the materials. Patterns select elements (including all their
isotopes), specific isotopes, or glob patterns:

```
$ ncrystal browse --atomdb He B10
LABEL  Z   A  ELEMENT  NATURAL     MASS    COHSL       COHXS  INCOHXS   SCATXS    ABSXS
He     2   0       He     True   4.0026  3.26548        1.34        0     1.34  0.00747
He3    2   3       He    False  3.01603     5.74     4.14032      1.6  5.74032     5333
He4    2   4       He    False   4.0026     3.26      1.3355        0   1.3355        0
B10    5  10        B    False  10.0129     -0.1  0.00125664        3  3.00126     3835
```

Here `cohsl` is the coherent scattering length [fm], the cross sections are
in barn (absorption at 2200m/s), and isotopes of hydrogen are labelled `H2`
and `H3`. The `-w`, `--sort`, `--columns`, `--csv`, `--json`, `--html`,
`--count`, and `--names` options work as for materials:

```
$ ncrystal browse --atomdb -w "absxs > 1000 and natural" --sort absxs --reverse \
                  --columns mass,absxs
LABEL     MASS  ABSXS
Gd     157.251  49700
Sm      150.35   5922
Eu     151.965   4530
Cd      112.41   2520

$ ncrystal browse --atomdb --csv > atomdb.csv      # the whole database
```

In Python, use the `AtomDBBrowser` class (or `query_atomdb()` for the raw
data):

```python
a = nb.AtomDBBrowser()
strong = a.where('absxs > 1000').sorted('absxs', reverse = True)
print( strong.table('mass,absxs') )
for e in nb.AtomDBBrowser('Li'):
    print( e['label'], e['mass'], e['absxs'] )   # entries are dictionaries
```

## Relation to nctool

The `--browse`, `--extract`, and `--plugins` options of `nctool` still work,
but `ncrystal browse` is the recommended tool for this. It can also list the
loaded plugins with `ncrystal browse --plugins`.

# Filters and windows: tables of macroscopic cross sections

*(Draft of a page for the NCrystal wiki.)*

Some parts of a neutron instrument simply attenuate the beam: filters (e.g. a
cooled beryllium or sapphire filter) and windows (e.g. aluminium). For these,
the scattered neutrons are usually not of interest. What matters is the
probability that a neutron passes through without interacting:

    T = exp(-Sigma * L)

where L is the path length in the material, and Sigma is the macroscopic total
cross section (scattering plus absorption) at the wavelength of the neutron.

NCrystal can provide Sigma as a table vs. the neutron wavelength, which is
suitable for fast lookups during a simulation, also on GPUs.

## The tables

* **Values:** the macroscopic total cross section, in 1/cm, vs. the neutron
  wavelength, in Aa.
* **Interpolation:** the table is piecewise linear. With linear interpolation,
  it reproduces the cross section of the material within a relative tolerance
  of 1e-3. (The tolerance is relative to the largest of the cross section and
  1e-12 barn per atom, so cross sections vanishing at wavelength 0 can be
  tabulated.)
* **Discontinuities:** these are mainly Bragg edges. Each is represented by two
  points with the same wavelength, where the first has the value just below
  the discontinuity and the second the value just above it.
* **Range:** the first point is at wavelength 0, with the limit of the cross
  section for wavelength -> 0. Beyond the last point, the last segment must be
  extrapolated linearly (clamped at 0).

The tables are created by evaluating the cross sections of the material, and
are verified against them. An error is raised if the cross section can not be
tabulated reliably. If a material uses a physics process which has not been
validated for such tables, a warning is emitted.

Only isotropic materials are supported, i.e. no oriented single crystals.
Multiphase materials are supported. A sapphire filter, where usually no
reflections satisfy the Bragg condition, can be approximated by disabling Bragg
diffraction: `"stdlib::Al2O3_sg167_Corundum.ncmat;bragg=0;temp=200K"`.

## Python

```python
import numpy as np
from NCrystal.filter import NCrystalFilter

f = NCrystalFilter('stdlib::Be_sg194.ncmat;temp=80K')
f.xsect(wl=4.0)                       # macroscopic cross section [1/cm] at 4 Aa
f.xsect(ekin=0.025)                   # ... or at 25 meV
wl = np.linspace(1.0, 10.0, 10)
transmission = np.exp(-f.xsect(wl=wl) * 5.0)   # through 5 cm of Be
wl_table, macroxs_table = f.table     # the table itself
```

The `.xsect` method accepts numbers or arrays, like the `.xsect` methods of the
scatter and absorption objects.

## C

```c
unsigned n;
double * wl;
double * macroxs;
ncrystal_filtertable( "stdlib::Be_sg194.ncmat;temp=80K", &n, &wl, &macroxs, NULL );
/* ... use the table ... */
ncrystal_dealloc_doubleptr( wl );
ncrystal_dealloc_doubleptr( macroxs );
```

The last argument (options) is reserved for future use, and must be NULL or
empty.

The example
[ncrystal_example_filter.c](https://github.com/mctools/ncrystal/blob/main/examples/ncrystal_example_filter.c)
shows how to use the table:
* It copies the table to arrays owned by the application (e.g. memory
  accessible on a GPU).
* It evaluates the table with a small self-contained function, which can be
  copied into other code. That function only reads the arrays and calls no
  other functions, so it can also be used on GPUs (e.g. with
  `#pragma acc routine seq` for OpenACC).

## McStas

The `NCrystal_filter` component (in McStas 3.x, from 2026) is a box or cylinder
of a material, which attenuates the beam passing through it, e.g.:

```
COMPONENT befilter = NCrystal_filter(cfg = "stdlib::Be_sg194.ncmat;temp=80K",
                                     xwidth = 0.05, yheight = 0.05, zdepth = 0.1)
```

* **How it works:** it creates the table in INITIALIZE, and only uses the table
  in TRACE, so it also works on GPUs.
* **Older NCrystal versions:** with versions which do not provide
  `ncrystal_filtertable`, the component creates a (larger and less verified)
  table itself, and prints a warning.
* **Example:** the instrument `NCrystal_filter_example` demonstrates it.
* **For other components:** the functions behind it
  (`mccode_init_ncrystal_xstable`, `mccode_eval_ncrystal_xstable` and
  `mccode_free_ncrystal_xstable` in `mccode-ncrystal-lib`) can also be used by
  other McStas components.

For the full NCrystal physics, including the scattered neutrons, use
`NCrystal_sample` or the Union components instead.

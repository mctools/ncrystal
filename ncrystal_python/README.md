NCrystal : A library for thermal neutron transport in crystals and other materials
----------------------------------------------------------------------------------

The `ncrystal-python` package provides the Python API and command-line scripts
for NCrystal. To function it needs to acess NCrystal shared library behind the
scenes, and as such it requires the `ncrystal-core` package to function.

In conda and Python environments, it is recommended that most users simply
install and depend on the package named `ncrystal` which itself depends on both
`ncrystal-core` and `ncrystal-python`.

Note that to support some esoteric installations, the `ncrystal-python` package
does not itself depend on the `ncrystal-core` package, which is another reason
that most users are recommended to simply install or depend only on the package
named `ncrystal`.

Installation of the `ncrystal-python` package more specifically provides (in
addition to what is found in the `ncrystal-core` package):

- The NCrystal python module.
  - This allows you to do `import NCrystal` and use the Python API of NCrystal.
- Command-line utilities.
  - These commands (e.g. `nctool`, `ncrystal_cif2ncmat`, ...) are actually
    themselves written in Python, and can also be accessed from the Python API
    itself via the `NCrystal.cli` module.
  - Note that the special command `ncrystal-config` is provided by the
    `ncrystal-core` package, _not_ the `ncrystal-python` package.
- Note that the `ncrystal-python` package can only be installed via Python
  installation tools like pip.
- For convenience, the `ncrystal-core` package has a dependency on `numpy`.

# Referencing NCrystal in scientific work

A very substantial effort went into developing NCrystal. If you use it for your
work, we would appreciate it if you would use the following primary reference in
your work:

  X.-X. Cai and T. Kittelmann, NCrystal: A library for thermal neutron
  transport, Computer Physics Communications 246 (2020) 106851,
  https://doi.org/10.1016/j.cpc.2019.07.015

For work benefitting from elastic physics (e.g. Bragg diffraction), we
furthermore request that you additionally also use the following reference in
your work:

  T. Kittelmann and X.-X. Cai, Elastic neutron scattering models
  for NCrystal, Computer Physics Communications 267 (2021) 108082,
  https://doi.org/10.1016/j.cpc.2021.108082

For work benefitting from our inelastic physics, we furthermore request that you
additionally also use the following reference in your work:

  X.-X. Cai, T. Kittelmann, et. al., "Rejection-based sampling of inelastic
  neutron scattering", Journal of Computational Physics 380 (2019) 400-407,
  https://doi.org/10.1016/j.jcp.2018.11.043


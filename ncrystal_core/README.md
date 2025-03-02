NCrystal : A library for thermal neutron transport in crystals and other materials
----------------------------------------------------------------------------------

The `ncrystal-core` package currently contains the NCrystal library and
everything else needed to use the C or C++ APIs of NCrystal. It notably does
_not_ include any Python modules and most of the command line scripts
(`ncrystal-config` being the exception), as those are provided in the
`ncrystal-python` package.

In conda and Python environments, it is recommended that most users simply
install and depend on the package named `ncrystal` which itself depends on both
`ncrystal-core` and `ncrystal-python`.

Installation of the `ncrystal-core` package more specifically provides:

- The NCrystal shared library (a binary file)
  - This binary file (.so/.dylib/.dll) is really the "engine" behind NCrystal.
- Header files
  - These files (.hh/.h) are included when using the C/C++ API of NCrystal.
- The `ncrystal-config` command-line utility.
  - This binary executable can be used to get information about NCrystal. As
    long as it is your PATH, you can use it to find all the other components
    of NCrystal.
- CMake configuration files:
  - Enables downstream CMake-based projects to use NCrystal via CMake
    `find_package(NCrystal)` calls (but see the file
    `<reporoot>/downstream_cmake/CMakeLists.txt` for how such `find_package`
    calls should actually be written to work with pip install'ed NCrystal.
- The NCrystal standard data library (NCMAT data).
  - These might be distributed as actual files, or (the default) embedded directly
    into the shared library.

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


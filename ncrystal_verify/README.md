Verification tool for NCrystal installations
--------------------------------------------

The `ncrystal-verify` package is non-binary package which is released along-side
each release of NCrystal, with an identical version number. It provides a single
command, `ncrystal-verify`, which can be used to validate that a particular
installation of NCrystal is functional and produces expected results. In case of
failure, `ncrystal-verify` ends with non-zero exit code.

Note that technically, the `ncrystal-verify` package depends on the
`ncrystal-python` package with the same version number. Therefore, doing
`pip install ncrystal-verify` should normally result in the appropriate version
of `ncrystal-verify` being installed.

Examples of usage:

1. Launch with no arguments to require ALL tests to succeed. Tests that are
   missing optional dependencies will count as failures.
   ```
   $> ncrystal-verify
   ```
2. Only run tests with no missing dependencies. Those with missing optional
   dependencies will be reported as skipped, and will NOT count as failures.
   ```
   $> ncrystal-verify -m all
   ```
3. Get specific usage instructions:
   ```
   $> ncrystal-verify --help
   ```

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


The `requirements_<name>.txt` files here correspond to the optional dependencies
of the ncrystal meta package in `ncrystal_metapkg/pyproject.toml`, with the
addition of any non-ncrystal non-optional dependencies needed by the
ncrystal-core and ncrystal-python packages (meaning in practice just numpy).

The `conda_*.yml` files are similar, but for conda-forge based conda
environments.

The `_all` files contains everything a user might need, while the `_devel` files
contains everything needed for NCrystal development. The `_all_combined` files
are simply a combination of the `_all` and `_devel` files.

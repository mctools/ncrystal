Standard installation
=====================

The recommended manner in which to install NCrystal is via prebuilt Conda or
Python packages which are available for most modern platforms (including Linux,
macOS, and Windows for both Intel, AMD and ARM platforms):

If using conda, NCrystal can be installed from the conda-forge channel, by
running the command:

```
conda install conda-forge::ncrystal
```

![Conda Version](https://img.shields.io/conda/v/conda-forge/ncrystal)
![Conda Platform](https://img.shields.io/conda/pn/conda-forge/ncrystal-core)
![Conda Downloads](https://img.shields.io/conda/dn/conda-forge/ncrystal)

For consistency, it is recommended that your conda environment ONLY uses
packages from the conda-forge channel.

If not using conda, NCrystal can be installed via a Python installation tool
like pip, via the command:

```
pip install ncrystal
```

![PyPI - Version](https://img.shields.io/pypi/v/ncrystal)
![PyPI - Downloads](https://img.shields.io/pypi/dm/ncrystal)

Some additional and optional packages are required for certain optional NCrystal
capabilities, which are not always needed by all users (for instance,
`matplotlib` is required for plotting). For that reason, some dependencies are
marked as optional, but can be added via Python packages optional dependency
support:

```
pip install ncrystal[plot] # for plotting
pip install ncrystal[cif] # for working with CIF files
pip install ncrystal[composer] # for creating new crystalline materials with
                               # NCMATComposer
pip install ncrystal[all] # all of the above
```

In addition to PyPI and conda-forge, other distribution channels for NCrystal
exist. Examples are platform specific packages for Debian/Ubuntu (![Debian
package](https://img.shields.io/debian/v/ncrystal) ![Ubuntu Package
Version](https://img.shields.io/ubuntu/v/ncrystal)) and FreeBSD (packaged as
[ncrystal](https://www.freshports.org/science/ncrystal) and
[py-ncrystal](https://www.freshports.org/science/py-ncrystal/). However, *these
are not tested directly by NCrystal developers*, so we can not make any
guarantee as to their correctness. Furthermore, PyPI and conda-forge packages of
NCrystal are usually ready within a few hours after a new NCrystal release,
ensuring that these packages represent the most updated version of NCrystal.

To test if an installation of NCrystal works, you should be able to run the
command `nctool --test`. You can also run `ncrystal-config -s`,
`nctool --plugins`, and `nctool --browse` to learn more about the version,
locations and capabilities of files in your installation.



Building from source (experts)
==============================

Advanced users might wish to build NCrystal from the sources in the upstream
repository (https://github.com/mctools/ncrystal). This is considered an
expert-only procedure, but it *is* based on standard tools such as CMake and
`pyproject.toml` so should hopefully be somewhat straight forward. A brief
high-level description will be given in the following, but feel free to reach
out and ask questions on https://github.com/mctools/ncrystal/discussions in case
something is not clear, or if you run into issues.

Note that to avoid potential conflicts, it is highly recommended that you
uninstall any existing NCrystal packages from your environment before performing
manual installations! If you are in a conda environment, you should ensure that
only the conda-forge channel is enabled, and that the following packages are
also installed into the conda environment: `cxx-compiler`, `c-compiler`,
`cmake`, and `make` (unix only), `numpy`, `pip`, and `python`. You can find
conda environment files for creating a suitable environment in
`<reporoot>/devel/reqs/conda*.yml`.


Package structure
-----------------

Before release 4.0.0, NCrystal was distributed in a single monolithic package
named `ncrystal`. It is still the case that users should generally just install
a package named `ncrystal`, but people wishing to build from sources need to be
aware that "under the hood", the various components are now actually distributed
in separate packages. Currently, the packages are:

- `ncrystal-core`:
  - The NCrystal shared library (a binary file)
    - This binary file (.so/.dylib/.dll) is really the "engine" behind NCrystal.
  - Header files
    - These files (.hh/.h) are included when using the C/C++ API of NCrystal.
  - The NCrystal standard data library (NCMAT data).
    - These might be distributed as actual files, or (the default) embedded directly
      into the shared library.
  - The `ncrystal-config` command-line utility.
    - This binary executable can be used to get information about NCrystal. As
      long as it is your PATH, you can use it to find all the other components
      of NCrystal.
  - CMake configuration files:
    - Enables downstream CMake-based projects to use NCrystal via CMake
      `find_package(NCrystal)` calls (but see the file
      `<reporoot>/examples/downstream_cmake/CMakeLists.txt` for how such
      `find_package` calls should actually be written to work with pip install'ed
      NCrystal.
- `ncrystal-python`:
  - The NCrystal python module.
    - This allows you to do `import NCrystal` and use the Python API of NCrystal.
  - Command-line utilities.
    - These commands (e.g. `nctool`, `ncrystal_cif2ncmat`, ...) are actually
      themselves written in Python, and can also be accessed from the Python API
      itself via the NCrystal.cli module.
  - Note that the `ncrystal-python` package can only be installed via Python
    installation tools like pip. Note additionally, that it does NOT technically
    depend on the `ncrystal-core` package, but of course does have a runtime
    requirement for it. This setup is intended to allow advanced deployment
    mechanisms where a manual CMake-based build of `ncrystal-core` is paired with
    a pip-based installation of `ncrystal-python`.
  - For convenience, this package has a dependency on `numpy`.
- `ncrystal`:
  - This is a meta-package with no actual contents, but it depends on both the
    `ncrystal-core` and `ncrystal-python` packages, and thus provide a complete
    and self-consistent NCrystal installation when users use a package manager
    to install "ncrystal".
  - Note that the pinning in this package ensures that `ncrystal-core` and
    `ncrystal-python` packages are always installed in the same versions. Thus, if
    you install `ncrystal` version 4.0.1, you will always get `ncrystal-core`
    version 4.0.1 and `ncrystal-python` version 4.0.1 as well.


Building ncrystal-core
----------------------

The `ncrystal-core` package can be built with either CMake or Python build tools
like pip, starting from the `./ncrystal_core` subdirectory which contains both a
CMakeLists.txt and a pyproject.toml file. The pyproject.toml file simply uses
scikit-build-core to wrap the CMake code under the hood, and is ignored if not
using Python build tools.

To build with Python build tools, simply point pip at the appropriate directory,
e.g.:

```
pip install ./ncrystal_core
```

It is beyond the scope of the present instructions to discuss the procedure of
manual building with CMake, although it should be noted that you can also modify
CMake options when installing via pip (e.g. `pip install
--config-settings=cmake.define.NCRYSTAL_ENABLE_DYNLOAD=OFF ./ncrystal-core`. Those
using CMake directly should simply note that it follows the usual procedure
(e.g. supports CMAKE_`INSTALL_PREFIX` etc.). To get started, here is nonetheless
how one might do it in BASH:

```
cd /path/to/ncrystalsource
mkdir build
cmake -S . -B build \
   -DCMAKE_INSTALL_PREFIX=/path/to/ncrystalinstall \
   -DNCRYSTAL_ENABLE_DATA=EMBED
cmake --build build --target install --config Release
```

You should then make sure that the directory where ncrystal-config was installed
is in your path, and you can use `which ncrystal-config` to verify that it is
so. Depending on your platform and options, this might look like:

```
   $> export PATH="/path/to/ncrystalinstall/bin:$PATH"
   $> ncrystal-config --summary
```

To see specific NCrystal CMake options, use CMake tools to do so, or refer to
the file `./ncrystal_core/cmake/modules/ncrystal_options.cmake`

Before reporting any issues with a manual CMake-based build of NCrystal, please
always try to completely clear your temporary CMake build directory and try
again.


Building ncrystal-python
------------------------

The `ncrystal-python` package can be built with Python build tools like pip,
starting from the `./ncrystal_python` subdirectory which contains an appropriate
pyproject.toml file:

```
pip install ./ncrystal_python
```

Note that runtime (i.e. when importing the NCrystal module later), you will get
failures if `ncrystal-core` was not somehow installed. More specifically, if
`ncrystal-core` was not installed with Python tools, then the command
`ncrystal-config` must be available in your PATH, so the Python modules can use
it to locate the binary shared library from `ncrystal-core`.


Quick and dirty monolithic installations
----------------------------------------

It is not actually recommended to do so, but if you are using Python and pip and
simply wish to quick build and install a particular revision of NCrystal into
that environment, and you don't want to have to care about the finer details of
`ncrystal-core` and `ncrystal-python` above, you can do so by simply pointing a
`pip install` command at the *root* (and thus *not* the `./ncrystal_core` or
`./ncrystal_python` folders) of the NCrystal git repository. Here are some
examples, which also show how one might choose a particular version tag og
branch name:

```
pip install git+https://github.com/mctools/ncrystal
pip install git+https://github.com/mctools/ncrystal@v4.0.0
pip install git+https://github.com/mctools/ncrystal@develop
```

This will install a single package called `ncrystal-monolithic-bundle` into your
environment, which contains the contents of both `ncrystal-core` and
`ncrystal-python`. If in a conda environment, be sure you have the appropriate
build tools installed *in that environment* (e.g. the `cxx-compiler`, `cmake`,
etc. packages -- please find the full list earlier in this file).

**WARNING:** Due to the limitations of Python build tools, any existing
`ncrystal-core` or `ncrystal-python` package already in the environment will not
be automatically uninstalled first when installing in this manner. So you MUST
do so manually BEFORE running the above command. If you are unsure what you have
already installed, you might use a command like `pip list` or `conda list` to
investigate.

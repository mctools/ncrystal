NCrystal unit tests
-------------------

This subdirectory adds a suite of unit tests for NCrystal. The primary way to
exercise them is to configure and build (but not install) the CMake project at
the root of the repository (i.e. one level up from here), with the option
`NCRYSTAL_ENABLE_TESTING=ON`, and then subsequently launch the tests via ctest.

Alternatively, one can also use the developer command `ncdevtool` from the
folder at `<reporoot>/devel/bin` (on unix one can inject it into the PATH by
sourcing the script at `<reporoot>/devel/setup.sh`). With that, one can simply
launch the tests via the command `ncdevtool cmake`. Alternatively, one can use
the simplebuild (cf. https://mctools.github.io/simplebuild/) system and use the
command `ncdevtool sb -t` to both build the NCrystal code and launch the tests.

Note that in addition to the tests here, the command `ncdevtool check` provides
several fast checks of the repository based on static code inspection.

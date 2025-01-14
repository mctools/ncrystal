NCrystal developer utilities
----------------------------

This subdirectory adds various utilities for NCrystal developers.

Most importantly, it provides the developer command `ncdevtool` from the folder
at `<reporoot>/devel/bin`. On unix one can inject it into the PATH by sourcing
the script at `<reporoot>/devel/setup.sh`, although one can of course also
simply invoke it with a full path (e.g. running `./devel/bin/ncdevtool
<options>` from the root of the repository).

Running `ncdevtool` with no arguments provides a full list of available modes,
and running `ncdevtool <modename> --help` will usually provide more descriptions
for how to use a particular mode.

For quick reference, a few of the more important modes are mentioned here:

* `ncdevtool check`: This runs several fast checks of the repository based on
  static code inspection and can catch common errors.
* `ncdevtool sb`: This uses the optional simplebuild
  (cf. https://mctools.github.io/simplebuild/) system to build the code. This is
  often more handy during development cycles than a pure CMake based build, and
  offers various features for code development and maintenance.
* `ncdevtool cmake`: Quickly configure and build the code with CMake, and either
  install it or run the tests from the tests/ subdirectory.
* `ncdevtool grep`: Search the code-base for strings.

In addition to the modes above, modes also exist to help with the editing of new
or existing C++ files, analysis of code dependencies, etc.

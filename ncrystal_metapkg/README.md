NCrystal : A library for thermal neutron transport in crystals and other materials
----------------------------------------------------------------------------------

![PyPI - Version](https://img.shields.io/pypi/v/ncrystal)
![Conda Version](https://img.shields.io/conda/vn/conda-forge/ncrystal)
![Conda Platform](https://img.shields.io/conda/pn/conda-forge/ncrystal-lib)

This is a source distribution of NCrystal, a library and associated tools which
enables calculations for Monte Carlo simulations of thermal neutrons in crystals
and other materials. Supported is a range of physics including both coherent,
incoherent, elastic and inelastic scatterings in a wide range of materials,
including crystal powders, mosaic single crystals, layered single crystals,
amorphous solids, and liquids. Multiphase materials or isotopically enriched
material are supported as well, and the framework furthermore supports
phase-contrast (SANS) physics. Written in C++, interfaces and infrastructure
facilitate integration into existing simulation frameworks such as OpenMC
(https://docs.openmc.org/), Geant4 (https://geant4.web.cern.ch/) or McStas
(http://mcstas.org/), as well as allowing direct usage from C++, C or Python
code or via command-line tools. While the C++ library is designed with a high
degree of flexibility in mind for developers, typical end-user configuration is
deliberately kept simple and uniform across various applications and APIs - this
for instance allows tuning and validation of a particular material configuration
to be performed in one tool before it is then deployed in another.

In addition to code and tools, the NCrystal distribution also includes a set of
validated data files, covering many crystals important at neutron scattering
facilities. For more information about the properties and validity of each file,
users are referred to the dedicated page at:

  https://github.com/mctools/ncrystal/wiki/Data-library

Supporting compilation with all modern C++ standards (C++11 and later), the code
has no third-party dependencies and is available under the highly liberal open
source Apache 2.0 license (see NOTICE and LICENSE files for usage conditions and
the INSTALL file for build and installation instructions). NCrystal was
developed in close collaboration by Xiao Xiao Cai (DTU, ESS) and Thomas
Kittelmann (ESS) and was supported in part by the European Union's Horizon 2020
research and innovation programme under grant agreement No 676548 (the
BrightnESS project) and 951782 (the HighNESS project).

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

The rest of this file gives a brief overview of the manners in which NCrystal
capabilities can be utilised. Further instructions and documentation, along with
the latest version of NCrystal, can be found at https://mctools.github.io/ncrystal/



Using the NCrystal installation from the command-line
-----------------------------------------------------

After installing NCrystal and having sourced the setup.sh script mentioned in
the INSTALL file, you can run any of the commands from the $NCRYSTALDIR/bin
directory, which includes example code as well as the "nctool" command. Start by
reading the usage instructions:

$> nctool --help

Assuming you chose to install data files, you can try to let NCrystal load one
of the data files found in $NCRYSTALDIR/data/ (or provide the absolute path to a
data file downloaded from https://github.com/mctools/ncrystal/wiki/Data-library)
and either dump the derived information to the terminal...:

$> nctool --dump 'Al_sg225.ncmat;temp=10C'

Note that this included a choice of temperature. If you leave it out, it will
usually default to room temperature (20C). You can also plot (powder)
cross-sections and sampled scatter angles with:

$> nctool 'Al_sg225.ncmat;temp=10C'



Using the NCrystal installation from C++, C or Python code
-----------------------------------------------------------------------------

If you wish to use NCrystal from Python code, there is no special setup needed,
assuming NCrystal was installed correctly (cf. INSTALL.md). If you on the other
hand wish to use NCrystal from your compiled C++ or C code, you must put
appropriate build flags. The recommended way is using CMake to do this (see next
section), but otherwise you must ensure that the NCrystal header files are in
your compiler's include path, and that the NCrystal library is linked
correctly. Here are some examples of how this could for instance be done, with a
C and a C++ app respectively:

```
export LDFLAGS="${LDFLAGS:-} -Wl,-rpath,$(ncrystal-config --show libdir) $(ncrystal-config --show libpath)"
export CFLAGS="${CFLAGS:-} -I$(ncrystal-config --show includedir)"
export CXXFLAGS="${CXXFLAGS:-} -I$(ncrystal-config --show includedir)"
cc -std=c11 ${LDFLAGS} ${CFLAGS} my_c_code.c -o my_c_app
c++ -std=c++17 ${LDFLAGS} ${CXXFLAGS} my_cpp_code.cpp -o my_cpp_app
```

Then, in your code you can access the relevant APIs with with statements like:

```
#include "NCrystal/NCrystal.hh"     // C++ code, core NCrystal
#include "NCrystal/ncrystal.h"      // C code
import NCrystal                     ## Python code
```

In the ./examples/ directory of your NCrystal distribution that you got after
downloading and unpacking the NCrystal source tar-ball, you will find small
examples of code using NCrystal. For C++/C, there is currently no documentation
beyond the header files and examples. In the case of Python, there is integrated
documentation available via the usual "help" function, accessed with:

```
import NCrystal
help(NCrystal)
```

There are also several jupyter-lab notebooks showcasing the NCrystal python API
at https://github.com/mctools/ncrystal-notebooks



Configuring CMake-based projects to use NCrystal
------------------------------------------------

Assuming NCrystal was built and installed via CMake, it is possible and
recommended for client projects to simply use NCrystal as a CMake package in
order to correctly build their C/C++ code which depends on the NCrystal C++ or C
APIs.

Depending on where NCrystal was installed on the system, it might be necessary
to let CMake know about it via the usual mechanisms (for instance passing
-DNCrystal_DIR=/path/to/ncrystalinstall as an argument to cmake on the command
line).

CMake code for a small project using NCrystal might look like the following
(assume that exampleapp.cc below includes the NCrystal/NCrystal.hh header):

  cmake_minimum_required(VERSION 3.16...3.31)
  project(MyExampleProject LANGUAGES CXX)
  execute_process( COMMAND "ncrystal-config" "--show" "cmakedir"
                   OUTPUT_VARIABLE NCrystal_DIR
                   OUTPUT_STRIP_TRAILING_WHITESPACE )
  find_package(NCrystal 4.0.0 REQUIRED)
  add_executable(exampleapp "${PROJECT_SOURCE_DIR}/exampleapp.cc")
  target_link_libraries( exampleapp NCrystal::NCrystal )
  install( TARGETS exampleapp DESTINATION bin )

Note that the "execute_process( ... )" command above is optional, but is
required before the code can work in an environment where the NCrystal CMake
modules are not automatically injected into the CMake package search path (this
notably includes NCrystal installed via "pip install ncrystal").



Using the NCrystal with Geant4
------------------------------

The NCrystal-Geant4 bindings are developed separately in the ncrystal-geant4
package. Refer to https://github.com/mctools/ncrystal-geant4 for more
information.



Using NCrystal with OpenMC
--------------------------

Using NCrystal materials in OpenMc is supported since OpenMC release 13.3, and
uses a nice simple syntax in the Python API:

```
mat = openmc.Material.from_ncrystal('Polyethylene_CH2.ncmat;temp=50C')
```

which when used in a complete OpenMC project, results in the following material
entry being added to the `materials.xml` produced:

```
<material cfg="Polyethylene_CH2.ncmat;temp=50C" id="1" temperature="323.15">
  <density units="g/cm3" value="0.92" />
  <nuclide ao="0.66656284" name="H1" />
  <nuclide ao="0.00010382666666666666" name="H2" />
  <nuclide ao="0.32964066666666664" name="C12" />
  <nuclide ao="0.003692666666666666" name="C13" />
</material>
```

Temperature, density and material composition were all created automatically
from the cfg-string, and the cfg-string itself was also encoded. Upon launching
the simulation with the OpenMC binary executable `openmc`, it will handle the
material as usual, except that low-energy neutron scattering physics (currently
defined as ($E<5eV$) will be provided by the algorithms in NCrystal.

A few issues might warrent attention:

1. If you try to assemble the above xml manually, it is rather unlikely that you
   will get the base densities and compositions right. It is safest to stick to
   let the Python API compose the xml for you.
2. After creation with `mat=openmc.Material.from_ncrystal(..)`, you can not use
   the usual OpenMC API to modify the material density, temperature, or
   composition. So be sure to reflect the final desired material inside the
   NCrystal cfg-string.
3. The OpenMC binaries must have been built with NCrystal support, or your job
   will fail once you launch the simulation (you can check for this by running
   the command `openmc -v`). Specifically (as documented on
   https://docs.openmc.org/en/stable/usersguide/install.html) you must supply
   the CMake flag `cmake -DOPENMC_USE_NCRYSTAL=on ..` (and make sure NCrystal is
   available already).  Note: we have agreement from OpenMC developers to enable
   NCrystal support by default in the conda-forge version of OpenMC. So in "the
   near future" (summer/fall 2023) conda users will always have NCrystal support
   available in OpenMC.

For more information, please consult the user guide at:

* https://docs.openmc.org/

In particular note the sections concerning installation and usage of NCrystal in
the sections:

* https://docs.openmc.org/en/stable/usersguide/install.html
* https://docs.openmc.org/en/stable/usersguide/materials.html



Using NCrystal with McStas
--------------------------

NOTE: The following discussion concerns the modern McStas 3 branch, and might in
particular not be 100% accurate for releases earlier than McStas 3.3 (probably
OK for v3.2 though).

You can use NCrystal in two ways in McStas. You can either use it for advanced
studies with the McStas Union sub-system through the NCrystal_process component,
or it can be used via the dedicated NCrystal_sample.comp which is simpler but
less feature rich. In any case, the McStas instrument file compilation will need
to build against NCrystal, and it uses the ncrystal-config command to figure out
the correct settings for doing so. Thus, you can always invoke "ncrystal-config
-s" to find out if you have the right NCrystal installation available and
active. Depending on how you installed McStas, NCrystal is most likely already
available. If not, you can try one of the following ways of enabling it:

```
$> conda install conda-forge::ncrystal [if you are in a conda-forge env]
$> python3 -mpip install ncrystal [for non-conda users]
```

It is beyond the scope for this README to provide a full documentation of
McStas, or the Union sub-system, but if you are using McStasScript to compose
your instruments, you can add NCrystal materials into your Union geometry using
code like:

```
from mcstasscript.tools.ncrystal_union import add_ncrystal_union_material
add_ncrystal_union_material(instr, name="myAl", cfgstr="Al_sg225.ncmat;temp=10C")
```

This creates the material and gives it the name "myAl", which you must later
attach to a particular Union volume, like for instance:

```
myvol.set_parameters(radius=0.01, yheight=0.01,
                     material_string='"myAl"', priority=1)
```

If you are instead hand-editing your instrument files, you can generate code
which defines Union materials from an NCrystal cfg-string by invoking:

```
$> python3 -mNCrystal.mcstasutils --union myAl 'Al_sg225.ncmat;temp=250K'
```

It should be noted that McStas 3.3 also provides a new SHELL syntax which can
also be used to faciliate this invocation from with a classic .instr file.

On the other hand, the dedicated NCrystal_sample.comp component, embeds NCrystal
material simulations into simple shapes (currently boxes, cylinders and
spheres), and can be used for components representing samples, filters or
monochromators, entrance windows, etc. The component is since McStas v3.3 part
of the McStas release itself, and can be used in a .instr file - for instance if
you wish to set up an r=1cm sphere with powdered sapphire you would write:

COMPONENT mysample = NCrystal_sample(cfg="Al2O3_sg167_Corundum.ncmat",radius=0.01)
AT (0, 0, 0) RELATIVE PREVIOUS

For more documentation about the NCrystal_sample component, run:

$> mcdoc NCrystal_sample

Or consult the documentation online at
https://www.mcstas.org/download/components/

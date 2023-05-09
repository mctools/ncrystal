################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2023 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

include_guard()

include( ncrystal_utils )

#All supported options are defined below. Additionally note that there is a
#hidden option (for experts) called NCRYSTAL_NOTOUCH_CMAKE_BUILD_TYPE" which can
#be set to ON to stop NCrystal from fiddling with the CMAKE_BUILD_TYPE.

bool_option( NCRYSTAL_ENABLE_EXAMPLES "Whether to build and install various examples." OFF )
bool_option( NCRYSTAL_ENABLE_MCSTAS   "Whether to install the NCrystal mcstas component and related scripts (not needed in McStas 3.2 and later)." ON )
enum_option( NCRYSTAL_ENABLE_GEANT4   "Whether to build the G4 hooks." OFF ON IFAVAILABLE )
bool_option( NCRYSTAL_ENABLE_PYTHON   "Whether to install the NCrystal python module and various python scripts (including ncrystal-config)." ON )
enum_option( NCRYSTAL_ENABLE_DATA     "Whether to include the standard data library files (possibly EMBED'ed into the binary)." ON OFF EMBED )
bool_option( NCRYSTAL_ENABLE_SETUPSH  "Whether to install setup.sh/unsetup.sh which users can source in order to use installation." ON )
bool_option( NCRYSTAL_MODIFY_RPATH    "Whether to try to set RPATH in installed binaries (if disabled all special RPATH handling is skipped)." ON )
enum_option( NCRYSTAL_ENABLE_DYNLOAD  "Enable dynamic library loading capabilities (for plugins)." IFAVAILABLE ON OFF )
bool_option( NCRYSTAL_SKIP_PYMODINST  "Prevents Python module installation and creates setup.py-based skeleton in <blddir>/ncrystal_pypkg (w/o _nclibpath.py file)." OFF )
enum_option( NCRYSTAL_BUILD_STRICT    "Stricter build (primarily for testing). Can optionally select specific C++ standard." OFF ON 11 14 17 20 )
bool_option( NCRYSTAL_ENABLE_CPACK    "Includes CPack and sets relevant meta-data." OFF )
bool_option( NCRYSTAL_QUIET           "Produce less status messages during configuration." OFF )
bool_option( NCRYSTAL_SKIP_INSTALL    "Set to prevent any installation targets from being created." OFF )

string_option(
  NCRYSTAL_BUILTIN_PLUGINS
  "Semicolon separated list of external plugins to statically embed in the NCrystal library (local paths to sources or git <repo_url:tag>)"
  ""
  )


#Older deprecated variables:
if ( NOT NCRYSTAL_NOLEGACYOPTS )
  include( ncrystal_legacyoptions )
endif()

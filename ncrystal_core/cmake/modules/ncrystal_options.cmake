################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2024 NCrystal developers                                   ##
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

bool_option( NCRYSTAL_ENABLE_EXAMPLES  "Whether to build and install various examples." "OFF" )
enum_option( NCRYSTAL_ENABLE_GEANT4    "Whether to build the G4 hooks." "OFF" "ON" "IFAVAILABLE" )
enum_option( NCRYSTAL_ENABLE_DATA      "Whether to include the standard data library files (possibly EMBED'ed into the binary)." "ON" "OFF" "EMBED" )
bool_option( NCRYSTAL_MODIFY_RPATH     "Whether to try to set RPATH in installed binaries (if disabled all special RPATH handling is skipped)." "ON" )
enum_option( NCRYSTAL_ENABLE_DYNLOAD   "Enable dynamic library loading capabilities (for plugins)." "IFAVAILABLE" "ON" "OFF" )
enum_option( NCRYSTAL_BUILD_STRICT     "Stricter build (primarily for testing). Can optionally select specific C++ standard." "OFF" "ON" "11" "14" "17" "20" "23" )
bool_option( NCRYSTAL_ENABLE_CPACK     "Includes CPack and sets relevant meta-data." "OFF" )
bool_option( NCRYSTAL_QUIET            "Produce less status messages during configuration." "OFF" )
bool_option( NCRYSTAL_SKIP_INSTALL     "Set to prevent any installation targets from being created." "OFF" )
enum_option( NCRYSTAL_ENABLE_THREADS   "Enable multithread usage (multithread safety is not affected by this)." "IFAVAILABLE" "ON" "OFF" )
bool_option( NCRYSTAL_ENABLE_TESTING   "Whether to enable the CTest-based tests in the tests subdirectory." "OFF" )
bool_option( NCRYSTAL_ENABLE_CFGAPP    "Whether to build and install the ncrystal-config command" "ON" )

string_option(
  NCRYSTAL_BUILTIN_PLUGINS
  "Semicolon separated list of external plugins to statically embed in the NCrystal library (local paths to sources or git <repo_url:tag>)"
  ""
  )

string_option(
  NCRYSTAL_NAMESPACE
  "Custom string to inject into library symbols (for complex environments with multple NCrystal installations)."
  ""
  )

#Older deprecated variables:
if ( NOT NCRYSTAL_NOLEGACYOPTS )
  include( ncrystal_legacyoptions )
endif()

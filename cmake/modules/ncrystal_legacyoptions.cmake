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

include( ncrystal_options )

#Support for deprecated options (for a little while!!)

add_deprecated_boolvar( BUILD_EXAMPLES  NCRYSTAL_ENABLE_EXAMPLES ON OFF )
add_deprecated_boolvar( INSTALL_MCSTAS  NCRYSTAL_ENABLE_MCSTAS ON OFF )
add_deprecated_boolvar( BUILD_G4HOOKS   NCRYSTAL_ENABLE_GEANT4 ON OFF )
add_deprecated_boolvar( INSTALL_PY      NCRYSTAL_ENABLE_PYTHON ON OFF )
add_deprecated_boolvar( INSTALL_SETUPSH NCRYSTAL_ENABLE_SETUPSH ON OFF )
add_deprecated_boolvar( MODIFY_RPATH    NCRYSTAL_MODIFY_RPATH ON OFF )
add_deprecated_boolvar( INSTALL_DATA    NCRYSTAL_ENABLE_DATA ON OFF )
add_deprecated_boolvar( EMBED_DATA      NCRYSTAL_ENABLE_DATA EMBED ON )
add_deprecated_boolvar( DISABLE_DYNLOAD NCRYSTAL_ENABLE_DYNLOAD OFF IFAVAILABLE )
add_deprecated_boolvar( NO_DIRECT_PYMODINST NCRYSTAL_SKIP_PYMODINST ON OFF )

if ( NOT "x${INSTALL_DATA}" STREQUAL "xUNSET" AND NOT "x${EMBED_DATA}" STREQUAL "xUNSET" )
  #both legacy data vars are set, we can't rely on the mapping by add_deprecated_boolvar:
  if ( EMBED_DATA AND INSTALL_DATA )
    #both on -> embed (old behaviour)
    set( NCRYSTAL_ENABLE_STDLIB EMBED )
  elseif( NOT EMBED_DATA AND NOT INSTALL_DATA )
    #both off -> off (obvious)
    set( NCRYSTAL_ENABLE_STDLIB OFF )
  elseif( EMBED_DATA AND NOT INSTALL_DATA )
    set( NCRYSTAL_ENABLE_STDLIB EMBED )
  elseif( NOT EMBED_DATA AND INSTALL_DATA )
    set( NCRYSTAL_ENABLE_STDLIB ON )
  endif()
endif()

set( BUILD_STRICT UNSET CACHE STRING "Obsolete alias for NCRYSTAL_BUILD_STRICT.")
if ( NOT "x${BUILD_STRICT}" STREQUAL "xUNSET" )
  check_deprecated_var( BUILD_STRICT NCRYSTAL_BUILD_STRICT )
  if ( NOT "x${BUILD_STRICT}" MATCHES "^x(OFF|ON|11|14|17|20)$" )
    message( FATAL_ERROR "The BUILD_STRICT variable is deprecated. Use NCRYSTAL_BUILD_STRICT instead (and provide it a correct value)" )
  endif()
  set( NCRYSTAL_BUILD_STRICT "${BUILD_STRICT}" )
endif()

set( BUILTIN_PLUGIN_LIST "UNSET" CACHE STRING  "Obsolete alias for NCRYSTAL_BUILTIN_PLUGINS.")
if ( NOT "x${BUILTIN_PLUGIN_LIST}" STREQUAL "xUNSET" )
  check_deprecated_var( BUILTIN_PLUGIN_LIST NCRYSTAL_BUILTIN_PLUGINS )
  set( NCRYSTAL_BUILTIN_PLUGINS "${BUILTIN_PLUGIN_LIST}" )
endif()

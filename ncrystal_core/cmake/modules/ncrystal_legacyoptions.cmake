
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
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

add_deprecated_boolvar( NCRYSTAL_ENABLE_MCSTAS "" "OFF" "" )
add_deprecated_boolvar( NCRYSTAL_ENABLE_SETUPSH "" "OFF" "" )
add_deprecated_boolvar( NCRYSTAL_ENABLE_PYTHON "" "OFF" "" )
add_deprecated_boolvar( NCRYSTAL_SKIP_PYMODINST "" "OFF" "" )
add_deprecated_boolvar( NCRYSTAL_ENABLE_SOVERSION "" "OFF" "" )

if ( DEFINED NCRYSTAL_BUILTIN_PLUGINS )
  message( FATAL_ERROR "The NCRYSTAL_BUILTIN_PLUGINS variable is no longer"
    " supported. All custom plugins should now be built separately from"
    " NCrystal ifself and loaded dynamically. See the examples/plugin"
    " directory for an updated plugin template." )
endif()
if ( DEFINED NCRYSTAL_ENABLE_GEANT4 )
  message( FATAL_ERROR "The NCRYSTAL_ENABLE_GEANT4 variable is no longer"
    " supported. The geant4 bindings have moved out of NCrystal and into"
    " https://github.com/mctools/ncrystal-geant4." )
endif()

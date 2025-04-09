
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

function( nccfgapp_create_ncapi_h resvar_includepath )
  #Generate ncconfig.h via both configure_file (for non-config specific
  #variables) and then file(GENERATE ..) (for generator expressions) based on
  #the ncapi.h.in template.
  set( srctemplate "${PROJECT_SOURCE_DIR}/include/NCrystal/ncapi.h.in" )
  set( tgtincdir "${PROJECT_BINARY_DIR}/autogen_include_ncapi" )
  set( tgtdir "${tgtincdir}/NCrystal" )
  set( tgtfile "${tgtdir}/ncapi.h" )
  file( MAKE_DIRECTORY "${tgtdir}" )
  set( t1 "1000000 * ${NCrystal_VERSION_MAJOR}" )
  set( t2 "1000 * ${NCrystal_VERSION_MINOR}" )
  set( t3 "${NCrystal_VERSION_PATCH}" )
  math(EXPR version_int "(${t1})+(${t2})+(${t3})" )
  set( ncapidefs "" )
  string( APPEND ncapidefs "#define NCRYSTAL_VERSION_MAJOR ${NCrystal_VERSION_MAJOR}\n" )
  string( APPEND ncapidefs "#define NCRYSTAL_VERSION_MINOR ${NCrystal_VERSION_MINOR}\n" )
  string( APPEND ncapidefs "#define NCRYSTAL_VERSION_PATCH ${NCrystal_VERSION_PATCH}\n" )
  string( APPEND ncapidefs "#define NCRYSTAL_VERSION ${version_int}\n" )
  string( APPEND ncapidefs "#define NCRYSTAL_VERSION_STR \"${NCrystal_VERSION}\"\n" )
  if ( NCRYSTAL_NAMESPACE )
    string( APPEND ncapidefs "#define NCRYSTAL_NAMESPACE_PROTECTION ${NCRYSTAL_NAMESPACE}\n")
  endif()
  if ( NCRYSTAL_ENABLE_TESTING )
    #Prevent interference in case an ncrystal-pluginmanager command is in the
    #environment:
    string( APPEND ncapidefs "#define NCRYSTAL_DISABLE_CMDLINEPLUGINMGR\n")
  endif()
  if ( WIN32 AND NCRYSTAL_WINEXPORTALL )
    #prevent __declspec(dllexport/import) since compiling with CMake
    #WINDOWS_EXPORT_ALL_SYMBOLS:
    string( APPEND ncapidefs "#define NCRYSTAL_PREVENT_WINDLLEXPORT\n")
  endif()
  set( NCRYSTAL_HOOK_FOR_ADDING_DEFINES
    " -- CMake definitions begin -- */\n\n${ncapidefs}\n/* -- CMake definitions end --" )
  configure_file( "${srctemplate}" "${tgtfile}" @ONLY )
  set( "${resvar_includepath}" "${tgtincdir}" PARENT_SCOPE )
endfunction()

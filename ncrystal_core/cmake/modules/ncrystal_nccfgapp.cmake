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

function( nccfgapp_create_ncconfig_h targetfile )
  #Generate ncconfig.h via file(GENERATE ..) from template. Supports generator
  #expressions in both targetfile and template contents.
  if ( NCRYSTAL_ENABLE_DATA )
    set( nccfgapp_has_data_01 "1" )
  else()
    set( nccfgapp_has_data_01 "0" )
  endif()
  if ( NCRYSTAL_ENABLE_EXAMPLES )
    set( nccfgapp_has_examples_01 "1" )
  else()
    set( nccfgapp_has_examples_01 "0" )
  endif()
  if ( NCRYSTAL_MODIFY_RPATH )
    set( nccfgapp_has_modrpath_01 "1" )
  else()
    set( nccfgapp_has_modrpath_01 "0" )
  endif()
  if ( _ncrystal_actual_enable_threads )
    set( nccfgapp_has_threads_01 "1" )
  else()
    set( nccfgapp_has_threads_01 "0" )
  endif()
  if ( NCRYSTAL_ENABLE_DATA STREQUAL "EMBED" )
    set( nccfgapp_has_dataembed_01 "1" )
  else()
    set( nccfgapp_has_dataembed_01 "0" )
  endif()
  if( _ncrystal_actual_disable_dynload )
    set( nccfgapp_has_dynload_01 "0" )
  else()
    set( nccfgapp_has_dynload_01 "1" )
  endif()
  set( t1 "1000000 * ${NCrystal_VERSION_MAJOR}" )
  set( t2 "1000 * ${NCrystal_VERSION_MINOR}" )
  set( t3 "${NCrystal_VERSION_PATCH}" )
  math(EXPR nccfgapp_intversion "(${t1})+(${t2})+(${t3})" )
  #First generated a configured version of the template:
  configure_file(
    "${NCrystal_SOURCE_DIR}/cmake/template_ncconfig.h"
    "${PROJECT_BINARY_DIR}/template_ncconfig.h.configured"
    @ONLY
  )
  #At build-time, generate the actual files needed for inclusion, expanding any
  #generator-expressions:
  file(
    GENERATE
    OUTPUT "${targetfile}"
    INPUT "${PROJECT_BINARY_DIR}/template_ncconfig.h.configured"
  )
endfunction()

function( create_ncrystal_config_app )
  set( autogenheader_name "ncconfig_autogen_$<CONFIG>.h")
  set ( workdir "${PROJECT_BINARY_DIR}/nccfgapp" )

  nccfgapp_create_ncconfig_h( "${workdir}/include/${autogenheader_name}" )
  configure_file(
    "${PROJECT_SOURCE_DIR}/include/NCrystal/internal/NCCFileUtils.hh"
    "${workdir}/include/NCCFileUtils.h"
    COPYONLY
  )
  configure_file(
    "${PROJECT_SOURCE_DIR}/src/NCCFileUtils.cc"
    "${workdir}/NCCFileUtils.c"
    COPYONLY
  )
  add_executable( ncrystal_cfgapp
    "${PROJECT_SOURCE_DIR}/app_config/main.c"
    "${workdir}/NCCFileUtils.c"
    "${workdir}/include/${autogenheader_name}"
    "${workdir}/include/NCCFileUtils.h"
  )
  set_target_properties(
    ncrystal_cfgapp PROPERTIES
    OUTPUT_NAME "ncrystal-config"
  )

  #C99 and C11 seems to be OK for most systems. Using C99 to be conservative
  #(could be reconsidered if we need it):
  set_target_properties( ncrystal_cfgapp PROPERTIES LANGUAGE C )
  target_compile_features( ncrystal_cfgapp PUBLIC c_std_99 )
  target_compile_definitions(
    ncrystal_cfgapp
    PUBLIC "NCCFGAPPHEADER=\"${autogenheader_name}\""
  )
  target_include_directories( ncrystal_cfgapp PUBLIC "${workdir}/include" )

  if ( SKBUILD_SCRIPTS_DIR )
    #fixme: verify if this actually ends up in PATH
    set( ncrystal_cfgapp_dest "${SKBUILD_SCRIPTS_DIR}" )
  else()
    set( ncrystal_cfgapp_dest "${NCrystal_BINDIR}" )
  endif()
  ncinstall( TARGETS ncrystal_cfgapp DESTINATION ${ncrystal_cfgapp_dest} )

endfunction()

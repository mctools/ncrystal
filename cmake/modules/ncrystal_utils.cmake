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

set( _${PROJECT_NAME}_all_opts "" )
set( _${PROJECT_NAME}_all_opts_vals "" )

macro( check_deprecated_var oldvarname newvarname )
  if ( NOT "x${${oldvarname}}" STREQUAL "xUNSET" )
    message( WARNING "The ${oldvarname} option is deprecated, please stop using it and instead use the new ${newvarname} option instead.")
    if ( NOT "x${${newvarname}}" STREQUAL "x${_${PROJECT_NAME}_vardef_${newvarname}}" )
      message( FATAL_ERROR "Do not set both ${newvarname} and its deprecated older incarnation ${oldvarname}" )
    endif()
  endif()
endmacro()

macro( add_deprecated_boolvar varname newvarname maptrueval mapfalseval )
  set( "${varname}" UNSET CACHE STRING "Deprecated option (use ${newvarname} instead)." )
  if ( NOT "x${${varname}}" STREQUAL "xUNSET" )
    check_deprecated_var( "${varname}" "${newvarname}" )
    if ( "${${varname}}" )
      set( ${newvarname} "${maptrueval}" )
    else()
      set( ${newvarname} "${mapfalseval}" )
    endif()
    set( ${varname} "_do_not_use_" )
  endif()
endmacro()

macro( bool_option varname docstr defaultvalue )
  list( APPEND _${PROJECT_NAME}_all_opts "${varname}" )
  if ( "${defaultvalue}"  )#NB: A simply "if ( defaultvalue )" does NOT work, because we are in a macro?
    set( _${PROJECT_NAME}_vardef_${varname} ON )
  else()
    set( _${PROJECT_NAME}_vardef_${varname} OFF )
  endif()
  option( "${varname}" "${docstr}" "${_${PROJECT_NAME}_vardef_${varname}}" )
  #map all values to ON/OFF for consistency
  if ( "${${varname}}" )
    set( ${varname} ON )
  else()
    set( ${varname} OFF )
  endif()
  list( APPEND _${PROJECT_NAME}_all_opts_vals "${${varname}}" )
endmacro()

macro( enum_option varname docstr defaultvalue )
  list( APPEND _${PROJECT_NAME}_all_opts "${varname}" )
  #To add an option MYOPTION which can only have values, A, B, C, D, or E, write
  #(the first option becomes the default value, here "A"):
  #  enum_option( MYOPTION "Bla bla." A B C D E )
  set( _${PROJECT_NAME}_vardef_${varname} "${defaultvalue}" )
  set(_tmp_allvals "${defaultvalue}" ${ARGN} )
  set( "${varname}" "${defaultvalue}" CACHE STRING "${docstr}" )
  set_property( CACHE ${varname} PROPERTY STRINGS ${_tmp_allvals} )
  if ( NOT "${${varname}}" IN_LIST _tmp_allvals )
    message( FATAL_ERROR "Option ${varname} has forbidden value \"${${varname}}\" (must be one of: ${_tmp_allvals})")
  endif()
  list( APPEND _${PROJECT_NAME}_all_opts_vals "${${varname}}" )
endmacro()

macro( string_option varname docstr defaultvalue )
  list( APPEND _${PROJECT_NAME}_all_opts "${varname}" )
  set( _${PROJECT_NAME}_vardef_${varname} "${defaultvalue}" )
  set( "${varname}" "${defaultvalue}" CACHE STRING "${docstr}" )
  list( APPEND _${PROJECT_NAME}_all_opts_vals "${${varname}}" )
endmacro()

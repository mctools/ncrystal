
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

set( _${PROJECT_NAME}_all_opts "" )
set( _${PROJECT_NAME}_all_opts_vals "" )

macro( check_deprecated_var oldvarname newvarname )
  if ( NOT "x${${oldvarname}}" STREQUAL "xUNSET" )
    if ( NOT "x${newvarname}" STREQUAL "x" )
      message( DEPRECATION "The ${oldvarname} option is deprecated, please stop using it and instead use the new ${newvarname} option instead.")
      if ( NOT "x${${newvarname}}" STREQUAL "x${_${PROJECT_NAME}_vardef_${newvarname}}" )
        message( FATAL_ERROR "Do not set both ${newvarname} and its deprecated older incarnation ${oldvarname}" )
      endif()
    else()
      message( DEPRECATION "The ${oldvarname} option is deprecated, please stop using it.")
    endif()
  endif()
endmacro()

macro( add_deprecated_boolvar varname newvarname maptrueval mapfalseval )
  if ( NOT "x${newvarname}" STREQUAL "x" )
    set( "${varname}" UNSET CACHE STRING "Deprecated option (use ${newvarname} instead)." )
  else()
    set( "${varname}" UNSET CACHE STRING "Deprecated option (simply ignored)." )
  endif()
  if ( NOT "x${${varname}}" STREQUAL "xUNSET" )
    check_deprecated_var( "${varname}" "${newvarname}" )
    if ( NOT "x${newvarname}" STREQUAL "x" )
      if ( "${${varname}}" )
        set( ${newvarname} "${maptrueval}" )
      else()
        set( ${newvarname} "${mapfalseval}" )
      endif()
      set( ${varname} "_do_not_use_" )
    endif()
  endif()
endmacro()

macro( bool_option varname docstr defaultvalue )
  list( APPEND _${PROJECT_NAME}_all_opts "${varname}" )
  if ( "${defaultvalue}"  )#NB: A simple "if ( defaultvalue )" does NOT work, because we are in a macro?
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

function( mctools_determine_strict_comp_flags resvar )
  string(
    SHA256 cacheid
    "${CMAKE_C_COMPILER_ID};${MCTOOLS_EXTRA_STRICT_COMP_FLAGS_MSVC}"
  )
  set( cachevar "MCTOOLS_STRICT_COMPFLAGS_CACHED_${cacheid}" )
  if ( DEFINED "${cachevar}" )
    set( "${resvar}" "${${cachevar}}" PARENT_SCOPE )
    return()
  endif()
  if(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set( flags "/WX" "/W4" )
    if ( MCTOOLS_EXTRA_STRICT_COMP_FLAGS_MSVC )
      #Meant for injecting specific flags to be disabled, like "/WD1234":
      list( APPEND flags ${MCTOOLS_EXTRA_STRICT_COMP_FLAGS_MSVC} )
    endif()
  else()
    set( flags -Wall -Wextra -pedantic -Werror )
  endif()
  include(CheckCCompilerFlag)
  # turn list into space separated string for check_c_compiler_flag:
  string( REPLACE ";" " " tmp "${flags}" )
  check_c_compiler_flag( "${tmp}" flagsok )#assuming same for C and C++
  if ( NOT flagsok )
    set( flags "" )
    message(WARNING "Could not enable strict compilation flags")
  endif()
  set(
    "${cachevar}" "${flags}"
    CACHE INTERNAL "caching strict flags" FORCE
  )
  set( ${resvar} "${flags}" PARENT_SCOPE )
endfunction()

function( mctools_apply_strict_comp_properties targetname )
  if( "${CMAKE_VERSION}" VERSION_GREATER_EQUAL "3.24" )
    set_target_properties(
      ${targetname} PROPERTIES COMPILE_WARNING_AS_ERROR ON
    )
  endif()
  mctools_determine_strict_comp_flags( strictflags )
  if ( strictflags )
    set_property(
      TARGET ${targetname}
      APPEND PROPERTY COMPILE_OPTIONS "${strictflags}"
    )
  endif()
endfunction()

function( mctools_detect_math_libs resvar )
  #Make sure we link in math functions correctly (typically the linker needs
  #libm on unix, but nothing on Windows). Sets the list of libraries to link in
  #resvar.
  string( SHA256 cacheid "${CMAKE_C_COMPILER_ID}" )
  set( cachevar "MCTOOLS_MATH_LIBS_CACHED_${cacheid}" )
  if ( DEFINED "${cachevar}" )
    set( "${resvar}" "${${cachevar}}" PARENT_SCOPE )
    return()
  endif()
  set( result "" )
  set(
    TMP_TESTLIBMSRC
    "#include <math.h>\nint main(int argc,char** argv) { (void)argv;double a=(exp)(argc+1.0); return (int)(a*0.1); }\n"
  )
  set( TMP_TESTDIR "${PROJECT_BINARY_DIR}/test_libm_${cacheid}" )
  file( WRITE ${TMP_TESTDIR}/test.c "${TMP_TESTLIBMSRC}" )
  try_compile(
    ALWAYS_HAS_MATH "${TMP_TESTDIR}" "${TMP_TESTDIR}/test.c"
  )
  if ( NOT ALWAYS_HAS_MATH )
    set( TMP_TESTDIR "${PROJECT_BINARY_DIR}/test_libm2_${cacheid}" )
    file(WRITE ${TMP_TESTDIR}/test.c "${TMP_TESTLIBMSRC}")
    try_compile(
      MATH_NEEDS_LIBM "${TMP_TESTDIR}" "${TMP_TESTDIR}/test.c" LINK_LIBRARIES m
    )
    if (MATH_NEEDS_LIBM)
      set( result "m" )
    else()
      message(
        FATAL_ERROR
        "Could not figure out link flags needed to enable math functions"
      )
    endif()
  endif()
  set(
    "${cachevar}" "${result}"
    CACHE INTERNAL "caching math libs" FORCE
  )
  set( ${resvar} "${result}" PARENT_SCOPE )
endfunction()

function( mctools_detect_extra_cflags resvar )
  #We assume that these flags are identical for C and C++.
  string( SHA256 cacheid "${CMAKE_C_COMPILER_ID};" )
  set( cachevar "MCTOOLS_EXTRA_CFLAGS_CACHED_${cacheid}" )
  if ( DEFINED "${cachevar}" )
    set( "${resvar}" "${${cachevar}}" PARENT_SCOPE )
    return()
  endif()
  include(CheckCCompilerFlag)
  set( flags "" )
  if ( "x${CMAKE_C_COMPILER_ID}" STREQUAL "xIntelLLVM" )
    #Intel defaults to the equivalent of -ffast-math, revert back to proper math:
    list( APPEND flags "-fp-model=precise" )
  elseif( "${CMAKE_VERSION}" VERSION_LESS "3.20" AND "x${CMAKE_C_COMPILER_ID}" STREQUAL "xClang" )
    #Older cmake had the llvm-based intel compiler classified as "Clang". In
    #this case, we check whether or not the -fp-model=precise flag is supported
    #or not:
    check_c_compiler_flag( -fp-model=precise tmp )
    if ( tmp )
      list( APPEND flags "-fp-model=precise" )
    endif()
  endif()
  check_c_compiler_flag( "-fno-math-errno" tmp )
  if ( tmp )
    list( APPEND flags "-fno-math-errno" )
  endif()
  set(
    "${cachevar}" "${flags}"
    CACHE INTERNAL "caching strict flags" FORCE
  )
  set( ${resvar} "${flags}" PARENT_SCOPE )
endfunction()

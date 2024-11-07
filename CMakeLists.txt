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

##################################################################################
#                                                                                #
# CMake file which can be used to compile and link all files distributed with    #
# NCrystal, and which provides CMake configuration files and (with               #
# NCRYSTAL_ENABLE_SETUPSH=ON) ncrystal_setup.sh/ncrystal_unsetup.sh for          #
# subsequent usage (also ncrystal_setup.sh/ncrystal_unsetup.sh).                 #
#                                                                                #
# One way to invoke cmake to build and install would be like this (run this from #
# a temporary build dir)                                                         #
#                                                                                #
#  $> cmake /path/to/sourcedir -DCMAKE_INSTALL_PREFIX=/path/to/installdir        #
#                                                                                #
# Followed by (replace the number 8 by the number of processes you want to       #
# use for the compilation):                                                      #
#                                                                                #
#  $> make install -j8                                                           #
#                                                                                #
# Written 2016-2023 by T. Kittelmann.                                            #
#                                                                                #
##################################################################################

cmake_minimum_required(VERSION 3.16...3.27)

# Respect value of CMAKE_BUILD_TYPE if already defined, otherwise fall back to
# Release. In any case, expose CMAKE_BUILD_TYPE as an explicit cache variable
# (gives drop-down list in gui). This must come before the call to
# project(..). We do not do this in case the generator is multi-cfg, and we also
# provide the hidden NCRYSTAL_NOTOUCH_CMAKE_BUILD_TYPE option to not do it.
#
if( NOT NCRYSTAL_NOTOUCH_CMAKE_BUILD_TYPE )
  get_property( gen_is_multicfg GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG )
  if ( NOT gen_is_multicfg )
    if( DEFINED CMAKE_BUILD_TYPE )
      set( _def_cbt ${CMAKE_BUILD_TYPE} )
    else()
      set( _def_cbt Release )
    endif()
    set( CMAKE_BUILD_TYPE ${_def_cbt}  CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo and MinSizeRel." )
    set_property( CACHE CMAKE_BUILD_TYPE
      PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel None )
  endif()
endif()

#Setup project:

set( _project_metadata LANGUAGES CXX C )

list( APPEND
  _project_metadata DESCRIPTION
  "Library for thermal neutron transport in crystals and other materials" )
list( APPEND
  _project_metadata HOMEPAGE_URL
  "https://github.com/mctools/ncrystal" )

project( NCrystal VERSION 3.9.7 ${_project_metadata} )

unset( _project_metadata )

if( NOT NCRYSTAL_NOTOUCH_CMAKE_BUILD_TYPE )
  if ( NOT gen_is_multicfg )
    if ( NOT CMAKE_BUILD_TYPE )
      #This can happen if parent project called the project(..) function before
      #doing the song and dance we did above.
      set(CMAKE_BUILD_TYPE Release)
    endif()
  endif()
endif()

# Set module path
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

#Define options:
include(ncrystal_options)

function( _ncrystal_option_force varname value )
  message( STATUS "Forcing ${varname}=\"${value}\"" )
  set( "${varname}" "${value}" PARENT_SCOPE )
endfunction()

#Hidden option needed for python setup.py/skbuild:
if ( DEFINED NCRYSTAL_SETUPPY_MODE )
  if ( NCRYSTAL_SETUPPY_MODE )
    set( NCRYSTAL_SETUPPY_MODE "ON" )
    message( STATUS "NCRYSTAL_SETUPPY_MODE detected => force-overriding values of some variables:")
    _ncrystal_option_force( NCRYSTAL_ENABLE_PYTHON "ON" )#plan is to move out of ncrystal repo
    _ncrystal_option_force( NCRYSTAL_ENABLE_EXAMPLES "OFF" )
    _ncrystal_option_force( NCRYSTAL_ENABLE_DYNLOAD "ON" )
    _ncrystal_option_force( NCRYSTAL_NOTOUCH_CMAKE_BUILD_TYPE "ON" )
    _ncrystal_option_force( NCRYSTAL_MODIFY_RPATH "OFF" )
    _ncrystal_option_force( NCRYSTAL_ENABLE_SETUPSH "OFF" )
    _ncrystal_option_force( NCRYSTAL_ENABLE_DATA "EMBED" )#otherwise the c++ library can not locate the files
    _ncrystal_option_force( NCRYSTAL_ENABLE_GEANT4 "OFF" )#plan is to move out of ncrystal repo
    _ncrystal_option_force( NCRYSTAL_ENABLE_THREADS "ON" )
    _ncrystal_option_force( NCRYSTAL_SKIP_PYMODINST "ON" )
    _ncrystal_option_force( CMAKE_INSTALL_LIBDIR "NCrystal/ncrystal_pyinst_data/lib" )
    _ncrystal_option_force( CMAKE_INSTALL_INCLUDEDIR "NCrystal/ncrystal_pyinst_data/include" )
    _ncrystal_option_force( CMAKE_INSTALL_DATADIR "NCrystal/ncrystal_pyinst_data/data" )
    _ncrystal_option_force( NCrystal_DATAFILESDIR "NCrystal/ncrystal_pyinst_data/stdlib_data" )
    #TODO: Postpone adding the following namespace until e.g. mcstas comps stop hardcoding "-lNCrystal":
    #_ncrystal_option_force( NCRYSTAL_NAMESPACE
    # "pip${PROJECT_VERSION_MAJOR}d${PROJECT_VERSION_MINOR}d${PROJECT_VERSION_PATCH}" )
  endif()
else()
  set( NCRYSTAL_SETUPPY_MODE "OFF" )
endif()
set( NCRYSTAL_ENABLE_PYTHON "OFF" )#FIXME
set( NCRYSTAL_SKIP_PYMODINST "OFF" )#FIXME

if ( WIN32 AND NCRYSTAL_ENABLE_TESTING AND NOT NCRYSTAL_SKIP_INSTALL )
  #Testing on windows needs us to export symbols from internal headers.
  message( WARNING "NCRYSTAL_ENABLE_TESTING on Windows implies NCRYSTAL_SKIP_INSTALL=ON for safety" )
  _ncrystal_option_force( NCRYSTAL_SKIP_INSTALL "ON" )
endif()

function( ncinstall )
  if ( NOT NCRYSTAL_SKIP_INSTALL )
    install( ${ARGN} )
  endif()
endfunction()

if ( NCRYSTAL_ENABLE_CPACK )
  set(CPACK_PACKAGE_CONTACT "ncrystal-developers@cern.ch")
  set(CPACK_NSIS_PACKAGE_NAME "${PROJECT_NAME} ${PROJECT_VERSION}")
  set(CPACK_NSIS_DISPLAY_NAME "${PROJECT_NAME} ${PROJECT_VERSION}")
  include(CPack)
endif()

if ( NCRYSTAL_SKIP_PYMODINST )
  if ( NOT NCRYSTAL_ENABLE_PYTHON )
    message( FATAL_ERROR "The NCRYSTAL_SKIP_PYMODINST flag requires -DNCRYSTAL_ENABLE_PYTHON=ON." )
  endif()
  if ( NCRYSTAL_ENABLE_SETUPSH )
    message( FATAL_ERROR "The NCRYSTAL_SKIP_PYMODINST flag requires -DNCRYSTAL_ENABLE_SETUPSH=OFF." )
  endif()
endif()

if ( NCRYSTAL_NAMESPACE )
  string( MAKE_C_IDENTIFIER "${NCRYSTAL_NAMESPACE}" tmp)
  string( TOLOWER "${tmp}" tmp)
  if ( NOT "x${NCRYSTAL_NAMESPACE}" STREQUAL "x${tmp}" )
    message(
      FATAL_ERROR
      "The NCRYSTAL_NAMESPACE must be a lowercased valid C identifier"
      "(meaning: alphanumeric, not leading with digit, lowercased)."
      )
  endif()
else()
  if ( NOT "x${NCRYSTAL_NAMESPACE}" STREQUAL "x" )
    message(
      FATAL_ERROR
      "Do not set NCRYSTAL_NAMESPACE to a value which CMake interprets as \"false\"."
      "To disable, instead set it to an empty string."
      )
  endif()
endif()


#Installation directories (try to follow standard conventions):
include(GNUInstallDirs)
set(NCrystal_BINDIR "${CMAKE_INSTALL_BINDIR}")#e.g. <prefix>/bin>
set(NCrystal_LIBDIR "${CMAKE_INSTALL_LIBDIR}")#e.g. <prefix>/lib>
set(NCrystal_INCDIR "${CMAKE_INSTALL_INCLUDEDIR}")#e.g. <prefix>/include>
set(NCrystal_DATAROOT "${CMAKE_INSTALL_DATADIR}/${PROJECT_NAME}")#e.g. <prefix>/share/NCrystal>
if( NOT NCrystal_DATAFILESDIR )
  set( NCrystal_DATAFILESDIR "${NCrystal_DATAROOT}/data")#e.g. <prefix>/share/NCrystal/data>
endif()
if( NOT NCrystal_PYPATH )
  set(NCrystal_PYPATH "${NCrystal_DATAROOT}/python")#e.g. <prefix>/share/NCrystal/python
endif()
if ( NOT NCrystal_CMAKEDIR )
  set(NCrystal_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")#e.g. <prefix>/lib/cmake/NCrystal>
endif()
if ( NOT NCRYSTAL_SKIP_PYMODINST )
  set(NCrystal_PYMODDIR "${NCrystal_PYPATH}/NCrystal")#e.g. <prefix>/share/NCrystal/python/NCrystal
endif()

#Get a few relative paths, mostly for expansion in various installed files (we
#use PROJECT_BINARY_DIR as prefix here, but it should not matter which as long
#as it is an absolute path):
if ( NOT NCRYSTAL_SKIP_PYMODINST )
  file(RELATIVE_PATH NCrystal_relpath_PYMODDIR2LIBDIR "${PROJECT_BINARY_DIR}/${NCrystal_PYMODDIR}" "${PROJECT_BINARY_DIR}/${NCrystal_LIBDIR}")
  file(RELATIVE_PATH NCrystal_relpath_CMAKEDIR2PYPATH   "${PROJECT_BINARY_DIR}/${NCrystal_CMAKEDIR}" "${PROJECT_BINARY_DIR}/${NCrystal_PYPATH}")
  file(RELATIVE_PATH NCrystal_relpath_BINDIR2PYPATH "${PROJECT_BINARY_DIR}/${NCrystal_BINDIR}"   "${PROJECT_BINARY_DIR}/${NCrystal_PYPATH}")
endif()
file(RELATIVE_PATH NCrystal_relpath_BINDIR2LIBDIR   "${PROJECT_BINARY_DIR}/${NCrystal_BINDIR}"   "${PROJECT_BINARY_DIR}/${NCrystal_LIBDIR}")
file(RELATIVE_PATH NCrystal_relpath_BINDIR2DATAROOT "${PROJECT_BINARY_DIR}/${NCrystal_BINDIR}"   "${PROJECT_BINARY_DIR}/${NCrystal_DATAROOT}")
file(RELATIVE_PATH NCrystal_relpath_BINDIR2CMAKEDIR "${PROJECT_BINARY_DIR}/${NCrystal_BINDIR}"   "${PROJECT_BINARY_DIR}/${NCrystal_CMAKEDIR}")
file(RELATIVE_PATH NCrystal_relpath_BINDIR2INCDIR   "${PROJECT_BINARY_DIR}/${NCrystal_BINDIR}"   "${PROJECT_BINARY_DIR}/${NCrystal_INCDIR}")
file(RELATIVE_PATH NCrystal_relpath_BINDIR2ROOT     "${PROJECT_BINARY_DIR}/${NCrystal_BINDIR}"   "${PROJECT_BINARY_DIR}/")
file(RELATIVE_PATH NCrystal_relpath_CMAKEDIR2ROOT   "${PROJECT_BINARY_DIR}/${NCrystal_CMAKEDIR}" "${PROJECT_BINARY_DIR}/")
file(RELATIVE_PATH NCrystal_relpath_CMAKEDIR2BINDIR   "${PROJECT_BINARY_DIR}/${NCrystal_CMAKEDIR}" "${PROJECT_BINARY_DIR}/${NCrystal_BINDIR}")
file(RELATIVE_PATH NCrystal_relpath_CMAKEDIR2LIBDIR   "${PROJECT_BINARY_DIR}/${NCrystal_CMAKEDIR}" "${PROJECT_BINARY_DIR}/${NCrystal_LIBDIR}")
file(RELATIVE_PATH NCrystal_relpath_CMAKEDIR2INCDIR   "${PROJECT_BINARY_DIR}/${NCrystal_CMAKEDIR}" "${PROJECT_BINARY_DIR}/${NCrystal_INCDIR}")
file(RELATIVE_PATH NCrystal_relpath_CMAKEDIR2DATAFILESDIR   "${PROJECT_BINARY_DIR}/${NCrystal_CMAKEDIR}" "${PROJECT_BINARY_DIR}/${NCrystal_DATAFILESDIR}")

#Dummy interface target for common properties. Note that the interface is always
#built with C++11 compatible compiler (which CMake will choose from a wide
#variety depending on platform and other targets,
#i.e. C++11/C++14/C++17/gnu++11/...). The NCRYSTAL_BUILD_STRICT option only affects
#private non-transitive properties.
add_library( ncrystal_common INTERFACE )
target_compile_features( ncrystal_common INTERFACE cxx_std_11 )

#Properties for executables and G4NCrystal library (can't transfer all
#properties via INTERFACE targets, so we need this variable-based workaround):
set( _nc_executable_rpath "")
set( _nc_libncg4_rpath "")

if ( NCRYSTAL_MODIFY_RPATH )
  #Set RPATH properties. For some annoying reason, this is not possible to do
  #via interface targets, so we have to use a variable-based workaround:

  if ( NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH )
    #TODO: Figure out if we really need this (perhaps only for geant4 targets?)
    set( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )
  endif()
  if( NOT APPLE )
    #Relocatable RPATHS: $ORIGIN in RPATH (including the $-char!!)  means the
    #location of the binary requiring the dependency:
    list( APPEND _nc_executable_rpath "$ORIGIN/${NCrystal_relpath_BINDIR2LIBDIR}" )
    list( APPEND _nc_libncg4_rpath $ORIGIN )
  else()
    #On OSX, rpaths are absolute paths (todo: revisit if this is still the case)
    get_filename_component( tmp "${CMAKE_INSTALL_PREFIX}/${NCrystal_LIBDIR}" ABSOLUTE)
    list( APPEND _nc_executable_rpath  "${tmp}" )
    list( APPEND _nc_libncg4_rpath "${tmp}" )
  endif()

  if ( DEFINED ENV{CONDA_PREFIX} AND EXISTS "$ENV{CONDA_PREFIX}/lib" )
    #In a conda environment we add $CONDA_PREFIX/lib to the rpath:
    list( APPEND _nc_executable_rpath "$ENV{CONDA_PREFIX}/lib" )
    list( APPEND _nc_libncg4_rpath "$ENV{CONDA_PREFIX}/lib" )
  endif()

  #Test if compiler supports -Wl,--disable-new-dtags. If it does, apply it
  #(otherwise RPATH sections in binaries become RUNPATH instead, which can be
  #overridden by users LD_LIBRARY_PATH (CMake>=3.14 is needed for LINK_OPTIONS on
  #try_compile and for the target_link_options function):
  #
  #NB: CMake 3.18 introduces CheckLinkerFlag module which we can eventually use
  #    instead of try_compile!!
  if( NOT MSVC )
    set(TMP_TESTDIR ${PROJECT_BINARY_DIR}/test_dtagflags)
    file(WRITE ${TMP_TESTDIR}/test.c "int main() { return 0; }\n")
    try_compile(LINKER_HAS_DTAGS "${TMP_TESTDIR}" "${TMP_TESTDIR}/test.c" LINK_OPTIONS -Wl,--disable-new-dtags)
    if (LINKER_HAS_DTAGS)
      #target_link_options(NCrystal PUBLIC "-Wl,--disable-new-dtags")
      target_link_options( ncrystal_common INTERFACE -Wl,--disable-new-dtags )
    endif()
  endif()

endif()

# Look for input files. Apparently "file(GLOB ...)" is frowned upon by some
# people. However, the only provided alternative (hardcode all your filenames)
# is rather unappealing. We glob for files, but apply the CONFIGURE_DEPENDS
# keyword when it is supported.

function(file_globsrc output_var pattern)
  file(GLOB tmp LIST_DIRECTORIES false CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/${pattern}" )
  set(${output_var} ${tmp} PARENT_SCOPE)
endfunction()

function(file_glob output_var pattern)
  file(GLOB tmp LIST_DIRECTORIES false CONFIGURE_DEPENDS "${pattern}" )
  set(${output_var} ${tmp} PARENT_SCOPE)
endfunction()

file_globsrc( SRCS_NC "ncrystal_core/src/*.cc")
file_globsrc( EXAMPLES_NC "examples/ncrystal_example_c*.c*")
file_globsrc( DATAFILES "data/*.ncmat")
file_globsrc( SRCS_NCPY "NCrystal/*.py")
file_globsrc( SRCS_PYSCRIPTS "scripts/ncrystal_*")
file_globsrc( tmp "scripts/nctool" )
list(APPEND SRCS_PYSCRIPTS ${tmp})
file_globsrc( HDRS_NCG4 "ncrystal_geant4/include/G4NCrystal/*.*")
file_globsrc( SRCS_NCG4 "ncrystal_geant4/src/*.cc")
file_globsrc( EXAMPLES_NCG4 "examples/ncrystal_example_g4*.cc")

#The ncapi.h header needs expansion of CMake variables, so we copy ALL header
#files to a temporary directory and edit just that one file:
#FIXME: Try to include in-place (apart from ncapi.h)
set( ncrystal_core_configured_includedir "${PROJECT_BINARY_DIR}/ncrystal_core_include_configured" )
file( MAKE_DIRECTORY "${ncrystal_core_configured_includedir}" )

#TODO: The trailing "/" on "ncrystal_core/include/" below is crucial. Is this behaviour portable?

file(
  COPY "ncrystal_core/include/"
  DESTINATION "${ncrystal_core_configured_includedir}"
  FILES_MATCHING
  PATTERN "*.h"
  PATTERN "*.hh"
)
if ( NOT EXISTS "${ncrystal_core_configured_includedir}/NCrystal/NCDefs.hh"
     OR NOT EXISTS "${ncrystal_core_configured_includedir}/NCrystal/internal/NCMath.hh" )
  message( FATAL_ERROR "Problems copying over include files" )
endif()
file( RENAME
  "${ncrystal_core_configured_includedir}/NCrystal/ncapi.h"
  "${ncrystal_core_configured_includedir}/NCrystal/ncapi.h.in"
  )
if ( NCRYSTAL_NAMESPACE )
  set( NCRYSTAL_CMAKE_HOOK_FOR_ADDING_NAMESPACE
    " -- CMake definitions begin -- */\n#define NCRYSTAL_NAMESPACE_PROTECTION ${NCRYSTAL_NAMESPACE}\n/* -- CMake definitions end --" )
else()
  set( NCRYSTAL_CMAKE_HOOK_FOR_ADDING_NAMESPACE
    " -- No definitions added by CMake -- " )
endif()
configure_file(
  "${ncrystal_core_configured_includedir}/NCrystal/ncapi.h.in"
  "${ncrystal_core_configured_includedir}/NCrystal/ncapi.h"
  @ONLY
)

file_glob( HDRS_NC "${ncrystal_core_configured_includedir}/NCrystal/*.hh" )
file_glob( HDRS_INTERNAL_NC "${ncrystal_core_configured_includedir}/NCrystal/internal/*.hh" )
file_glob( tmp "${ncrystal_core_configured_includedir}/NCrystal/*.h" )
list( APPEND HDRS_NC ${tmp} )
if ( NOT HDRS_NC OR NOT HDRS_INTERNAL_NC )
  message( FATAL_ERROR "Unexpected problem globbing for header files in ${ncrystal_core_configured_includedir}/" )
endif()

if ( NCRYSTAL_ENABLE_SETUPSH AND NOT NCRYSTAL_ENABLE_PYTHON )
  message( WARNING "NCRYSTAL_ENABLE_SETUPSH is not possible when NCRYSTAL_ENABLE_PYTHON is OFF (forcing NCRYSTAL_ENABLE_SETUPSH=OFF)." )
  set( NCRYSTAL_ENABLE_SETUPSH OFF )
endif()

include(CheckCXXCompilerFlag)

set( NC_STRICT_COMP_FLAGS -Wall -Wextra -pedantic -Werror )
if ( NCRYSTAL_BUILD_STRICT )
  string( REPLACE ";" " " tmp "${NC_STRICT_COMP_FLAGS}" )# check_cxx_compiler_flag needs single argument
  check_cxx_compiler_flag( "${tmp}" ncrystal_compiler_supports_strict_comp_flags )
endif()

set( STRICT_CPPSTD OFF )
set( STRICT_CSTD OFF )
if ( NCRYSTAL_BUILD_STRICT )
  #We also want to test the C-example with strict C standard (90, 99, 11). For
  #simplicity we simply map: C++11->C90, C++14->C99, C++17->C11. It is mainly
  #done to ensure unit test coverage.
  if ( NCRYSTAL_BUILD_STRICT STREQUAL "11" )
    set( STRICT_CSTD 90 )
    set( STRICT_CPPSTD 11 )
  elseif ( NCRYSTAL_BUILD_STRICT STREQUAL "14" )
    set( STRICT_CSTD 99 )
    set( STRICT_CPPSTD 14 )
  elseif ( NCRYSTAL_BUILD_STRICT STREQUAL "17" )
    set( STRICT_CSTD 11 )
    set( STRICT_CPPSTD 17 )
  elseif ( NCRYSTAL_BUILD_STRICT STREQUAL "20" )
    set( STRICT_CSTD 11 )#NB: Perhaps we will get C22 at some point?
    set( STRICT_CPPSTD 20 )
  endif()
endif()

set( ncrystal_extra_private_compile_options "" )
if ( "x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xIntelLLVM" )
  #Intel defaults to the equivalent of -ffast-math, revert back to proper math:
  set( ncrystal_extra_private_compile_options -fp-model=precise )
elseif( "${CMAKE_VERSION}" VERSION_LESS "3.20" AND "x${CMAKE_CXX_COMPILER_ID}" STREQUAL "xClang" )
  #Older cmake had the llvm-based intel compiler classified as "Clang". In
  #this case, we check whether or not the -fp-model=precise flag is supported
  #or not:
  check_cxx_compiler_flag( -fp-model=precise ncrystal_compiler_supports_fpmodelprecise )
  if ( ncrystal_compiler_supports_fpmodelprecise )
    set( ncrystal_extra_private_compile_options -fp-model=precise )
  endif()
endif()

function(set_target_common_props targetname)
  #Set private non-transitive properties. If strict builds are enabled, this can
  #enforce no warnings and compilation with a specific standards.
  if ( NCRYSTAL_BUILD_STRICT AND ncrystal_compiler_supports_strict_comp_flags )
    target_compile_options( ${targetname} PRIVATE ${NC_STRICT_COMP_FLAGS} )
  endif()
  if ( STRICT_CPPSTD )
    set_target_properties( ${targetname} PROPERTIES CXX_STANDARD ${STRICT_CPPSTD} CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
  endif()
  if ( STRICT_CSTD )
    set_target_properties( ${targetname} PROPERTIES C_STANDARD ${STRICT_CSTD} C_STANDARD_REQUIRED ON C_EXTENSIONS OFF)
  endif()
  #Always disallow M_PI and friends in our own code (they are not portable):
  target_compile_definitions( ${targetname} PRIVATE NCRYSTAL_NO_CMATH_CONSTANTS )
  if ( ncrystal_extra_private_compile_options )
    target_compile_options( ${targetname} PRIVATE ${ncrystal_extra_private_compile_options} )
  endif()
endfunction()

#Can we use -fno-math-errno when building the NCrystal library?
check_cxx_compiler_flag( -fno-math-errno ncrystal_compiler_supports_nomatherrno_flag )

set( ncrystal_threads_lib OFF )
if ( NOT NCRYSTAL_ENABLE_THREADS STREQUAL "OFF" )
  #set(THREADS_PREFER_PTHREAD_FLAG ON)
  find_package(Threads)
  if( Threads_FOUND )
    set( ncrystal_threads_lib "Threads::Threads" )
  else()
    if ( NOT NCRYSTAL_ENABLE_THREADS STREQUAL "IFAVAILABLE" )
      message( FATAL_ERROR "NCRYSTAL_ENABLE_THREADS set to ON but failed to enable thread support (set to IFAVAILABLE or OFF to proceed without it)." )
    endif()
  endif()
endif()

#NCrystal library and header files, including optional built-in modules if enabled:
add_library( NCrystal SHARED ${SRCS_NC} )
if ( ncrystal_compiler_supports_nomatherrno_flag )
  target_compile_options( NCrystal PRIVATE -fno-math-errno )
endif()

if ( ncrystal_threads_lib )
  message( STATUS "Thread support will be enabled.")
  set( _ncrystal_actual_enable_threads ON )
  target_link_libraries( NCrystal PRIVATE ${ncrystal_threads_lib} )
else()
  message( STATUS "Thread support will not be enabled.")
  set( _ncrystal_actual_enable_threads OFF )
  target_compile_definitions( NCrystal PRIVATE NCRYSTAL_DISABLE_THREADS )
endif()

if ( NCRYSTAL_ENABLE_SOVERSION AND NOT UNIX )
  message( WARNING "NCRYSTAL_ENABLE_SOVERSION ignored on non-unix platforms for now." )
  set( NCRYSTAL_ENABLE_SOVERSION OFF )
endif()

if ( DEFINED NCRYSTAL_OVERRIDE_SOVERSION )
  #For now this is a hidden feature, needed for unit tests in an environment
  #where NCrystal might already be installed:
  if ( NOT NCRYSTAL_ENABLE_SOVERSION )
    message( FATAL_ERROR "Do not use NCRYSTAL_OVERRIDE_SOVERSION with NCRYSTAL_ENABLE_SOVERSION=OFF" )
  endif()
  if ( NOT UNIX )
    message( FATAL_ERROR "NCRYSTAL_OVERRIDE_SOVERSION is only supported on unix builds for now." )
  endif()
  set( NCrystal_SOVERSION "${NCRYSTAL_OVERRIDE_SOVERSION}" )
else()
  #NOTE: Only the major number here, so we better make sure we take care of ABI
  #      compatibility within a major number in the future! We could have been
  #      more strict here, but then the assumed version pinning to major version
  #      numbers in conda and pypi would be a hassle instead:
  set( NCrystal_SOVERSION "${PROJECT_VERSION_MAJOR}" )
endif()
if ( NCRYSTAL_ENABLE_SOVERSION )
  set( _ncrystal_lib_version "${PROJECT_VERSION}")
  set_target_properties( NCrystal PROPERTIES VERSION "${_ncrystal_lib_version}" SOVERSION ${NCrystal_SOVERSION} )
else()
  #Make sure to export empty string to .in files:
  set( NCrystal_SOVERSION "" )
endif()

#Make sure client code will use at least c++11:
target_compile_features( NCrystal INTERFACE cxx_std_11 )

#Libname for ncrystal-config, _nclibpath, etc.:
set( _tmp_basiclibname "NCrystal" )
if ( NCRYSTAL_NAMESPACE )
  set( _tmp_basiclibname "${_tmp_basiclibname}-${NCRYSTAL_NAMESPACE}" )
endif()
set_target_properties( NCrystal PROPERTIES OUTPUT_NAME "${_tmp_basiclibname}" )
if ( NCRYSTAL_ENABLE_SOVERSION )
  #NB: We disabled this on windows for now.
  if( APPLE )
    #On osx, we get libNCrystal.3.8.0.dylib
    set( NCrystal_LIBNAME "${CMAKE_SHARED_LIBRARY_PREFIX}${_tmp_basiclibname}.${_ncrystal_lib_version}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
  else()
    #On linux, we get libNCrystal.so.3.8.0
    set( NCrystal_LIBNAME "${CMAKE_SHARED_LIBRARY_PREFIX}${_tmp_basiclibname}${CMAKE_SHARED_LIBRARY_SUFFIX}.${_ncrystal_lib_version}" )
  endif()
else()
    #No version, just libNCrystal<namespace>.so or libNCrystal<namespace>.dylib
  set( NCrystal_LIBNAME "${CMAKE_SHARED_LIBRARY_PREFIX}${_tmp_basiclibname}${CMAKE_SHARED_LIBRARY_SUFFIX}" )
endif()

#put libname in .txt file for downstream (e.g. conda) post-processing:
file(WRITE ${PROJECT_BINARY_DIR}/cfg_ncrystal_libname.txt "${NCrystal_LIBNAME}")

# Dynamic library loading:

# NCRYSTAL_ENABLE_DYNLOAD must be ON or OFF or IFAVAILABLE:
if ( NCRYSTAL_ENABLE_DYNLOAD STREQUAL "IFAVAILABLE" )
  if ( UNIX )
    set( _ncrystal_actual_disable_dynload OFF )
  else()
    set( _ncrystal_actual_disable_dynload ON )
  endif()
else()
  if ( NCRYSTAL_ENABLE_DYNLOAD STREQUAL ON )
    set( _ncrystal_actual_disable_dynload OFF )
  else()
    set( _ncrystal_actual_disable_dynload ON )
  endif()
endif()

if ( NOT _ncrystal_actual_disable_dynload AND CMAKE_DL_LIBS AND UNIX )
  #Cf. https://gitlab.kitware.com/cmake/cmake/-/merge_requests/1642 for why we
  #only do this on UNIX.
  target_link_libraries( NCrystal PRIVATE ${CMAKE_DL_LIBS} )
endif()
if ( _ncrystal_actual_disable_dynload )
  target_compile_definitions( NCrystal PRIVATE NCRYSTAL_DISABLE_DYNLOADER )
endif()

set_target_common_props( NCrystal )
target_link_libraries( NCrystal PRIVATE ncrystal_common )
target_include_directories(
  NCrystal PRIVATE "${PROJECT_SOURCE_DIR}/ncrystal_core/src"
  PUBLIC $<BUILD_INTERFACE:${ncrystal_core_configured_includedir}>
         $<INSTALL_INTERFACE:${NCrystal_INCDIR}>
)

#Make sure we link in math functions correctly (typically the linker needs libm on unix, but nothing on Windows).
set(TMP_TESTLIBMSRC "#include <math.h>\nint main(int argc,char** argv) { (void)argv;double a=(exp)(argc+1.0); return (int)(a*0.1); }\n")
set(TMP_TESTDIR ${PROJECT_BINARY_DIR}/test_libm)
file(WRITE ${TMP_TESTDIR}/test.c "${TMP_TESTLIBMSRC}")
try_compile(ALWAYS_HAS_MATH "${TMP_TESTDIR}" "${TMP_TESTDIR}/test.c")
if (NOT ALWAYS_HAS_MATH)
  set(TMP_TESTDIR ${PROJECT_BINARY_DIR}/test_libm2)
  file(WRITE ${TMP_TESTDIR}/test.c "${TMP_TESTLIBMSRC}")
  try_compile(MATH_NEEDS_LIBM "${TMP_TESTDIR}" "${TMP_TESTDIR}/test.c" LINK_LIBRARIES m)
  if (MATH_NEEDS_LIBM)
    target_link_libraries(NCrystal PRIVATE m)
  else()
    message( FATAL_ERROR "Could not figure out link flags needed to enable math functions" )
  endif()
endif()

#Embed requested external plugins directly into NCrystal:
function(parse_pluginentry entry)
  #Local entries are simply paths, while remote urls always contain at least one
  #semicolon and might look like:
  #
  # https://github.com/mctools/ncplugin-MyPlugin.git::develop <- repo url + '::' + git ref
  # https://github.com/mctools/ncplugin-MyPlugin.git          <- repo url (defaults to git ref 'main')
  # mctools:MyPlugin[::gitref]                                <- same in condensed form (no '/', a single ':')
  if(NOT entry MATCHES ":")
    set(pluginentry_islocal ON PARENT_SCOPE)
    return()
  endif()
  set(pluginentry_islocal OFF PARENT_SCOPE)
  string(REPLACE "::" ";" parts "${entry}")
  list(LENGTH parts nparts)
  if (nparts EQUAL 1)
    list(APPEND parts main)#default git ref tag
  elseif(NOT nparts EQUAL 2)
    message( FATAL_ERROR "Invalid syntax of entry in NCRYSTAL_BUILTIN_PLUGINS: ${entry}" )
  endif()
  list(GET parts 0 remoteurl)
  list(GET parts 1 remoteref)
  #condensed form has no slashes and a single ':'
  if(NOT remoteurl MATCHES "/")
    #condensed form has no slashes and a single ':'
    string(REPLACE ":" ";" urlparts "${remoteurl}")
    list(LENGTH urlparts nurlparts)
    if (nurlparts EQUAL 2)
      list(GET urlparts 0 remoteaccount)
      list(GET urlparts 1 remotereponameending)
      set(remoteurl "https://github.com/${remoteaccount}/ncplugin-${remotereponameending}.git")
    endif()
  endif()
  set(pluginentry_remoteurl ${remoteurl} PARENT_SCOPE)
  set(pluginentry_remoteref ${remoteref} PARENT_SCOPE)
endfunction()

if ( NCRYSTAL_ENABLE_DATA STREQUAL "EMBED" )
  #Embed data (needs to invoke python process to generate C++ code from .ncmat files)

  target_compile_definitions(NCrystal PRIVATE NCRYSTAL_STDCMAKECFG_EMBED_DATA_ON)

  #We must find python3 interpreter. The FindPython3 module is only available
  #from CMake 3.12 and is reported to not work well until CMake 3.15.5.
  message( STATUS "Looking for python interpreter (needed when NCRYSTAL_ENABLE_DATA=EMBED).")
  find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter)
  #Generate C++ code from the .ncmat files:
  execute_process(
    COMMAND "${Python3_EXECUTABLE}"
    "-BI"
    "${PROJECT_SOURCE_DIR}/NCrystal/_ncmat2cpp_impl.py"
    "-n" "NCrystal::AutoGenNCMAT::registerStdNCMAT"
    "--regfctname" "NCrystal::internal::registerEmbeddedNCMAT(const char*,const char*)"
    "-o" "${PROJECT_BINARY_DIR}/autogen_ncmat_data.cc" ${DATAFILES} RESULT_VARIABLE status
  )
  if(status AND NOT status EQUAL 0)
    message(FATAL_ERROR "Failure while trying to invoke ncmat2cpp (needed since NCRYSTAL_ENABLE_DATA=EMBED).")
  endif()
  target_sources(NCrystal PRIVATE "${PROJECT_BINARY_DIR}/autogen_ncmat_data.cc")#too late to just append to SRCS_NC
  message(STATUS "Generated autogen_ncmat_data.cc with embedded NCMAT data (will be compiled into the NCrystal library)." )
endif()

if ( NCRYSTAL_ENABLE_DATA STREQUAL "ON" )
  #Hardwiring NCRYSTAL_DATADIR in the binary, although handled with
  #NCRYSTAL_DATADIR env var in ncrystal_setup.sh. The environment variable makes
  #the installation relocatable (at least for users sourcing the installed
  #ncrystal_setup.sh):
  target_compile_definitions(NCrystal PRIVATE "NCRYSTAL_DATADIR=${CMAKE_INSTALL_PREFIX}/${NCrystal_DATAFILESDIR}")
  ncinstall(FILES ${DATAFILES} DESTINATION ${NCrystal_DATAFILESDIR})
endif()

ncinstall(TARGETS NCrystal EXPORT NCrystalTargets DESTINATION ${NCrystal_LIBDIR} )
ncinstall(FILES ${HDRS_NC} DESTINATION ${NCrystal_INCDIR}/NCrystal)
ncinstall(FILES ${HDRS_INTERNAL_NC} DESTINATION ${NCrystal_INCDIR}/NCrystal/internal)

#Examples:
if ( NCRYSTAL_ENABLE_EXAMPLES AND EXAMPLES_NC )
  foreach(ex ${EXAMPLES_NC})
    get_filename_component(exbn "${ex}" NAME_WE)
    add_executable(${exbn} "${ex}")
    set_target_common_props( ${exbn} )
    target_link_libraries(${exbn} NCrystal ncrystal_common )
    if ( NOT "x${_nc_executable_rpath}" STREQUAL "x" )
      set_target_properties(${exbn} PROPERTIES INSTALL_RPATH "${_nc_executable_rpath}" )
    endif()
    ncinstall(TARGETS ${exbn} DESTINATION ${NCrystal_BINDIR} )
  endforeach()
endif()

#Python interface:
if (NCRYSTAL_ENABLE_PYTHON)
  #NB: We don't actually require Python3 to be available, since we are just
  #copying over files. In principle a user can install python3 after installing
  #NCrystal.
  # --> Python module files:
  if ( NOT NCRYSTAL_SKIP_PYMODINST )
    ncinstall(FILES ${SRCS_NCPY} DESTINATION ${NCrystal_PYMODDIR})
    #autogenerated _nclibpath.py with relative location to library:
    file(WRITE "${PROJECT_BINARY_DIR}/_nclibpath.py.in"
      "#File autogenerated by NCrystal's CMakeLists.txt:\n"
      "liblocation='${NCrystal_relpath_PYMODDIR2LIBDIR}/${NCrystal_LIBNAME}'\n")
    if ( NCRYSTAL_NAMESPACE )
      file(APPEND "${PROJECT_BINARY_DIR}/_nclibpath.py.in"
        "ncrystal_namespace_protection='${NCRYSTAL_NAMESPACE}'\n")
    endif()
    configure_file("${PROJECT_BINARY_DIR}/_nclibpath.py.in" "${PROJECT_BINARY_DIR}/_nclibpath.py" COPYONLY)
    ncinstall( FILES "${PROJECT_BINARY_DIR}/_nclibpath.py" DESTINATION ${NCrystal_PYMODDIR})
  else()
    if ( NCRYSTAL_SETUPPY_MODE )
      #Static python files will be grapped from <src>/NCrystal by setup.py
      #itself, so we do nothing. However, _chooks.py also needs the information
      #about libname and namespace in _nclibpath.py:
      message( STATUS "Adding _nclibpath.py needed for setup.py installation.")
      file(WRITE "${PROJECT_BINARY_DIR}/_nclibpath.py.in"
        "#File autogenerated by NCrystal's CMakeLists.txt:\n"
        "import pathlib\n"
        "liblocation=(pathlib.Path(__file__).parent / 'ncrystal_pyinst_data' / 'lib' / '${NCrystal_LIBNAME}').absolute()\n")
      if ( NCRYSTAL_NAMESPACE )
        file(APPEND "${PROJECT_BINARY_DIR}/_nclibpath.py.in"
          "ncrystal_namespace_protection='${NCRYSTAL_NAMESPACE}'\n")
      endif()
      configure_file("${PROJECT_BINARY_DIR}/_nclibpath.py.in" "${PROJECT_BINARY_DIR}/_nclibpath.py" COPYONLY)
      ncinstall( FILES "${PROJECT_BINARY_DIR}/_nclibpath.py" DESTINATION "NCrystal" )
    else()
      #Similar but preparing a custom directory on which to run pip install during conda build:
      #scripts:
      file( MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/ncrystal_pypkg/NCrystal" )
      configure_file( "${PROJECT_SOURCE_DIR}/cmake/template_setup.py.in" "${PROJECT_BINARY_DIR}/ncrystal_pypkg/setup.py" @ONLY )
      file(COPY ${SRCS_NCPY} DESTINATION "${PROJECT_BINARY_DIR}/ncrystal_pypkg/NCrystal")
      message( STATUS "Will not install Python modules directly. Instead leaving pypkg skeleton in: ${PROJECT_BINARY_DIR}/ncrystal_pypkg/")
      message( STATUS "Note that the skeleton is incomplete (missing the ${PROJECT_BINARY_DIR}/ncrystal_pypkg/NCrystal/_nclibpath.py file).")
    endif()
  endif()

  # --> Python scripts:
  ncinstall(PROGRAMS ${SRCS_PYSCRIPTS} DESTINATION ${NCrystal_BINDIR})
  if ( NCRYSTAL_ENABLE_EXAMPLES )
    ncinstall(PROGRAMS "${PROJECT_SOURCE_DIR}/examples/ncrystal_example_py" DESTINATION ${NCrystal_BINDIR})
  endif()

endif()

#G4NCrystal
set( NCRYSTAL_GEANT4 OFF )
if ( NOT NCRYSTAL_ENABLE_GEANT4 STREQUAL "OFF" )
  find_package(Geant4)
  if( Geant4_FOUND )
    set( NCRYSTAL_GEANT4 ON )
  else()
    if ( NOT NCRYSTAL_ENABLE_GEANT4 STREQUAL "IFAVAILABLE" )
      message( FATAL_ERROR "NCRYSTAL_ENABLE_GEANT4 set to ON but failed to enable Geant4 support (set to IFAVAILABLE or OFF to proceed without Geant4 support)." )
    endif()
  endif()
endif()

if ( NCRYSTAL_GEANT4 )
  set( _tmp_basiclibnameg4 "G4NCrystal" )
  add_library(G4NCrystal SHARED ${SRCS_NCG4})#NB: No SOVERSION number since we plan to migrate this library out of the NCrystal repo.
  if ( NCRYSTAL_NAMESPACE )
    #Since we didn't implement the namespace for the G4NCrystal C++ namespace as
    #well, it is best to error out here. Once we have migrated to new and better
    #(and standalone) NCrystal-Geant4 hooks, they should certainly also have an
    #option for namespace protection there.
    message( WARNING "The NCRYSTAL_NAMESPACE feature is not actually fully implemented for the NCRYSTAL_GEANT4 extension (please get in touch if you have a different need)" )
    set( _tmp_basiclibnameg4 "${_tmp_basiclibnameg4}-${NCRYSTAL_NAMESPACE}" )
    set_target_properties( G4NCrystal PROPERTIES OUTPUT_NAME "${_tmp_basiclibnameg4}" )
  endif()
  set(G4NCrystal_LIBNAME "${CMAKE_SHARED_LIBRARY_PREFIX}${_tmp_basiclibnameg4}${CMAKE_SHARED_LIBRARY_SUFFIX}")
  #Make sure client code will use at least c++11:
  target_compile_features( G4NCrystal INTERFACE cxx_std_11 )
  set_target_common_props( G4NCrystal )
  target_include_directories(G4NCrystal
    PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/ncrystal_geant4/include>
    $<INSTALL_INTERFACE:${NCrystal_INCDIR}>
    )

  ncinstall(TARGETS G4NCrystal EXPORT G4NCrystalTargets DESTINATION ${NCrystal_LIBDIR} )

  ncinstall(FILES ${HDRS_NCG4} DESTINATION ${NCrystal_INCDIR}/G4NCrystal)
  #Transfer G4 flags uncovered by find_package(Geant4) call:
  #TODO: This should really use Geant4 targets once all known users move to
  #newer Geant4 release (G4 10.6?):
  target_compile_definitions(G4NCrystal PUBLIC ${Geant4_DEFINITIONS})
  target_link_libraries( G4NCrystal PUBLIC NCrystal ${Geant4_LIBRARIES} PRIVATE ncrystal_common )
  if ( NOT "x${_nc_libncg4_rpath}" STREQUAL "x" )
    set_target_properties( G4NCrystal PROPERTIES INSTALL_RPATH "${_nc_libncg4_rpath}" )
  endif()
  target_include_directories(G4NCrystal SYSTEM PUBLIC ${Geant4_INCLUDE_DIRS})
  set(Geant4_CXX_FLAGS_aslist ${Geant4_CXX_FLAGS})
  separate_arguments(Geant4_CXX_FLAGS_aslist)
  target_compile_options(G4NCrystal PUBLIC ${Geant4_CXX_FLAGS_aslist})
  #Check if compiler supports -Wno-overloaded-virtual. If so, we add it as a
  #public flag on the G4NCrystal target. This is because of Geant4 adds
  #-Wno-overloaded-virtual, and Geant4 v10.4 headers gives that warning with
  #gcc8/gcc9.
  check_cxx_compiler_flag( -Wno-overloaded-virtual ncrystal_compiler_supports_no_overloaded_virtual_flag )
  if ( ncrystal_compiler_supports_no_overloaded_virtual_flag )
    target_compile_options( G4NCrystal PUBLIC -Wno-overloaded-virtual )
  endif()

  #examples if needed:
  if ( NCRYSTAL_ENABLE_EXAMPLES AND EXAMPLES_NCG4 )
    foreach(ex ${EXAMPLES_NCG4})
      get_filename_component(exbn "${ex}" NAME_WE)
      add_executable(${exbn} "${ex}")
      set_target_common_props( ${exbn} )
      target_link_libraries(${exbn} G4NCrystal ncrystal_common )
      if ( NOT "x${_nc_executable_rpath}" STREQUAL "x" )
        set_target_properties(${exbn} PROPERTIES INSTALL_RPATH "${_nc_executable_rpath}" )
      endif()
      ncinstall(TARGETS ${exbn} DESTINATION ${NCrystal_BINDIR} )
    endforeach()
  endif()
endif()

#Package configuration files for downstream cmake projects:
ncinstall( EXPORT NCrystalTargets FILE NCrystalTargets.cmake NAMESPACE NCrystal:: DESTINATION ${NCrystal_CMAKEDIR} )
add_library(NCrystal::NCrystal ALIAS NCrystal)#always alias namespaces locally

if ( NCRYSTAL_GEANT4 )
  #Add G4NCrystal target in separate file:
  ncinstall( EXPORT G4NCrystalTargets FILE G4NCrystalTargets.cmake NAMESPACE NCrystal:: DESTINATION ${NCrystal_CMAKEDIR} )
  add_library(NCrystal::G4NCrystal ALIAS G4NCrystal)#always alias namespaces locally
endif()

if ( NCRYSTAL_ENABLE_SETUPSH )
  configure_file( "${PROJECT_SOURCE_DIR}/cmake/template_setup.sh.in" "${PROJECT_BINARY_DIR}/generated_setup.sh" @ONLY )
  configure_file( "${PROJECT_SOURCE_DIR}/cmake/template_unsetup.sh.in" "${PROJECT_BINARY_DIR}/generated_unsetup.sh" @ONLY )
  ncinstall( FILES "${PROJECT_BINARY_DIR}/generated_setup.sh" DESTINATION . RENAME ncrystal_setup.sh )
  ncinstall( FILES "${PROJECT_BINARY_DIR}/generated_unsetup.sh" DESTINATION . RENAME ncrystal_unsetup.sh )
endif()

if ( NCRYSTAL_ENABLE_TESTING )
  #NB: The enable_testing() call must occur before ANY add_subdirectory calls
  #(including those further down this file).
  if( "${CMAKE_VERSION}" VERSION_LESS "3.28" )
    #NB: Keep this minimum version value synchronised with the one on the top of
    #tests/CMakeLists.txt
    message( FATAL_ERROR "NCRYSTAL_ENABLE_TESTING requires cmake >= 3.28" )
  else()
    #We must export all symbols in libNCrystal.so, even those from internal
    #headers - since we also use internal files in our tests (for that reason we
    #also forced NCRYSTAL_SKIP_INSTALL above). We use
    #NCRYSTAL_PREVENT_WINDLLEXPORT to avoid explicitly exporting any symbols in
    #this case, since that triggers spurious warnings from VisualStudio.
    if ( WIN32 )
      set_target_properties(NCrystal PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
      target_compile_definitions(NCrystal PUBLIC NCRYSTAL_PREVENT_WINDLLEXPORT)
    endif()
    enable_testing()
    add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/tests")
  endif()
endif()

set(all_plugin_names_lc "stdscat;stdabs;stdncmat")
set(NCrystal_builtin_plugin_names "")
set(ncrystal_iplugin 0)
foreach(pluginentry ${NCRYSTAL_BUILTIN_PLUGINS})
  message( STATUS "Trying to add plugin: ${pluginentry}")
  parse_pluginentry("${pluginentry}")
  if (NOT pluginentry_islocal)
    #need unique id when using FetchContent:
    math(EXPR ncrystal_iplugin "${ncrystal_iplugin}+1")
    string(MD5 pluginid "${pluginentry}")
    set(pluginid "pl${pluginid}_${ncrystal_iplugin}")
    if (NOT COMMAND FetchContent_Declare)
      include(FetchContent)
    endif()
    message( STATUS "-- Trying to install plugin via remote git repo: ${pluginentry_remoteurl} (gitref ${pluginentry_remoteref})")
    FetchContent_Declare( ${pluginid}
      GIT_REPOSITORY "${pluginentry_remoteurl}"
      GIT_TAG        "${pluginentry_remoteref}"
      GIT_SHALLOW ON
      )
    FetchContent_GetProperties(${pluginid})
    if(NOT ${pluginid}_POPULATED)
      FetchContent_Populate(${pluginid})
      set(pluginlocalsrcdir "${${pluginid}_SOURCE_DIR}")
      message(STATUS "-- Fetched sources to ${pluginlocalsrcdir}")
    endif()
  else()
    message(STATUS "-- Trying to install plugin via local path: ${pluginentry}")
    if(NOT IS_ABSOLUTE ${pluginentry})
      message(FATAL_ERROR "Local path in NCRYSTAL_BUILTIN_PLUGINS is not an absolute paths: ${pluginentry}")
    endif()
    set(pluginlocalsrcdir "${pluginentry}")
  endif()
  if(NOT EXISTS "${pluginlocalsrcdir}/ncplugin_name.txt" OR NOT EXISTS "${pluginlocalsrcdir}/CMakeLists.txt")
    message(FATAL_ERROR "Entry in NCRYSTAL_BUILTIN_PLUGINS does not appear to contain proper sources (missing ncplugin_name.txt or CMakeLists.txt file): ${pluginentry}")
  endif()
  #Get plugin name from ncplugin_name.txt:
  file(STRINGS "${pluginlocalsrcdir}/ncplugin_name.txt" NCPlugin_NAME LIMIT_COUNT 1)
  string(STRIP "${NCPlugin_NAME}" NCPlugin_NAME)
  #Check unique-ness of plugin name (at least the builtin ones):
  string(TOLOWER "${NCPlugin_NAME}" NCPlugin_NAME_lowercase)
  if (NCPlugin_NAME_lowercase IN_LIST all_plugin_names_lc)
    message(FATAL_ERROR "Multiple plugins have the same name (must be unique and not clash with standard plugins): ${NCPlugin_NAME}.")
  endif()
  list(APPEND all_plugin_names_lc "${NCPlugin_NAME_lowercase}")
  #Add subdirectory after setting a few variables needed by the plugin's
  #CMakeLists.txt, and clearing a few ones we expect it to provide:
  set(NCPLUGIN_DEVMODE OFF)
  set(NCPLUGIN_ASBUILTIN ON)
  unset(NCPLUGIN_SRCFILES)
  unset(NCPLUGIN_DATAFILES)
  unset(NCPLUGIN_COMPILEDEFS)
  unset(NCPLUGIN_INCDIRS)
  set(plugin_binary_dir "${PROJECT_BINARY_DIR}/ncplugins/${NCPlugin_NAME}")
  add_subdirectory(${pluginlocalsrcdir} "${plugin_binary_dir}")
  if (NOT NCPLUGIN_SRCFILES OR NOT NCPlugin_NAME)
    message(FATAL_ERROR "Problem adding plugin (did not provide both NCPLUGIN_SRCFILES and NCPlugin_NAME variables): ${pluginentry}")
  endif()
  if (NCPLUGIN_COMPILEDEFS)
    set_property(SOURCE ${NCPLUGIN_SRCFILES} APPEND PROPERTY COMPILE_DEFINITIONS ${NCPLUGIN_COMPILEDEFS})
  endif()
  if (NCPLUGIN_INCDIRS)
    set_property(SOURCE ${NCPLUGIN_SRCFILES} APPEND PROPERTY INCLUDE_DIRECTORIES ${NCPLUGIN_INCDIRS})
  endif()
  if (NCPLUGIN_DATAFILES)
    foreach(df ${NCPLUGIN_DATAFILES})
      #Check that data file name follows convention: ncplugin-<pluginname>_*.ncmat
      get_filename_component(dfbn "${df}" NAME)
      if(NOT dfbn MATCHES "ncplugin-${NCPlugin_NAME}_.+\.ncmat")
        message(FATAL_ERROR "Problem adding plugin: name of datafile ${dfbn} does not have required form: ncplugin-${NCPlugin_NAME}_*.ncmat")
      endif()
      list(APPEND DATAFILES "${df}")
    endforeach()
  endif()
  list(APPEND NCrystal_builtin_plugin_names ${NCPlugin_NAME})
  target_sources(NCrystal PRIVATE ${NCPLUGIN_SRCFILES})
  message(STATUS "-- Configured plugin ${NCPlugin_NAME}")
endforeach()
if (NCrystal_builtin_plugin_names)
  target_compile_definitions(NCrystal PRIVATE NCRYSTAL_HAS_BUILTIN_PLUGINS)
  #Add generated .cc file which can load the plugins.
  set(ncplugcc "${PROJECT_BINARY_DIR}/autogen_ncplugins.cc.in")
  file(WRITE ${ncplugcc} "//Autogenerated file\n#include <iostream>\n#include \"NCrystal/NCPluginMgmt.hh\"\n")
  foreach(NCPlugin_NAME ${NCrystal_builtin_plugin_names})
    file(APPEND ${ncplugcc} "namespace NCrystalPlugin_${NCPlugin_NAME} { void registerPlugin(); }\n")
  endforeach()
  file(APPEND ${ncplugcc} "namespace NCRYSTAL_NAMESPACE {\n  void provideBuiltinPlugins() {\n")
  foreach(NCPlugin_NAME ${NCrystal_builtin_plugin_names})
    file(APPEND ${ncplugcc} "    Plugins::loadBuiltinPlugin(\"${NCPlugin_NAME}\",NCrystalPlugin_${NCPlugin_NAME}::registerPlugin);\n")
  endforeach()
  file(APPEND ${ncplugcc} "  }\n")
  file(APPEND ${ncplugcc} "}\n")
  configure_file("${ncplugcc}" "${PROJECT_BINARY_DIR}/autogen_ncplugins.cc" COPYONLY)
  target_sources(NCrystal PRIVATE "${PROJECT_BINARY_DIR}/autogen_ncplugins.cc")
  message( STATUS "Generated autogen_ncplugins.cc for enabling embedded plugins (will be compiled into the NCrystal library)." )
endif()

#NB: Important we export the config file AFTER all variables used have been set!
include(CMakePackageConfigHelpers)
write_basic_package_version_file( "${PROJECT_BINARY_DIR}/NCrystalConfigVersion.cmake"
                                  VERSION ${NCrystal_VERSION} COMPATIBILITY SameMajorVersion )

configure_file( "${PROJECT_SOURCE_DIR}/cmake/NCrystalConfig.cmake.in"
                "${PROJECT_BINARY_DIR}/NCrystalConfig.cmake" @ONLY )

ncinstall( FILES "${PROJECT_BINARY_DIR}/NCrystalConfigVersion.cmake" "${PROJECT_BINARY_DIR}/NCrystalConfig.cmake"
         DESTINATION ${NCrystal_CMAKEDIR} )

if (NCRYSTAL_ENABLE_PYTHON)
  #NB: Important we do this next call AFTER all variables used in ncrystal-config.in have been set!
  configure_file( "${PROJECT_SOURCE_DIR}/cmake/ncrystal-config.in" "${PROJECT_BINARY_DIR}/ncrystal-config" @ONLY )
  ncinstall( PROGRAMS "${PROJECT_BINARY_DIR}/ncrystal-config" DESTINATION ${NCrystal_BINDIR} )
endif()

if ( NOT NCRYSTAL_QUIET )
  foreach( optname ${_NCrystal_all_opts} )
    message( STATUS "NCrystal-cfg: ${optname}=${${optname}}" )
  endforeach()
endif()


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

#Tests require newer cmake (NB: Synchronize minimum value with that found in the
#check for NCRYSTAL_ENABLE_TESTING support in the ncrystal_core/CMakeLists.txt
#file):
cmake_minimum_required(VERSION 3.28...3.30)

project( NCrystalTests LANGUAGES CXX C )

if ( NOT TARGET NCrystal )
  message( FATAL_ERROR "Can only test as part of the NCrystal build" )
endif()
set( nctests_ncrystallib_target "NCrystal" )

enable_testing()

set(
  CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_LIST_DIR}/cmake"
  "${CMAKE_CURRENT_LIST_DIR}/../ncrystal_core/cmake/modules"
)
include( mctools_testutils )

#We want to symlink so that "import NCrystalDev" works instead of "import
#NCrystal" for the tests (to reduce conflicts with NCrystal installed in the
#environment). We copy instead as a fall-back (e.g. on Windows).

set( nctest_pypath "${PROJECT_BINARY_DIR}/ncrystal_tests_pypath" )
file( MAKE_DIRECTORY "${nctest_pypath}" )
file( CREATE_LINK
  "${NCrystalRoot_SOURCE_DIR}/ncrystal_python/src/NCrystal"
  "${nctest_pypath}/NCrystalDev"
  COPY_ON_ERROR
  SYMBOLIC
)

set( nctest_pyenvmod "" )

list(
  APPEND nctest_pyenvmod
  "PYTHONPATH=path_list_prepend:${nctest_pypath}"
  "PYTHONPATH=path_list_prepend:${PROJECT_SOURCE_DIR}/pypath"
)

list(
  APPEND nctest_pyenvmod
  "NCRYSTAL_LIB=set:$<TARGET_FILE:${nctests_ncrystallib_target}>"
)

if ( NOT NCRYSTAL_ENABLE_DATA STREQUAL "EMBED" )
  list( APPEND nctest_pyenvmod
        "NCRYSTAL_DATADIR=set:${NCrystalRoot_SOURCE_DIR}/data" )
endif()

list(
  APPEND nctest_pyenvmod
  "PYTHONIOENCODING=set:UTF-8"
  "PYTHONLEGACYWINDOWSSTDIO=set:UTF-8"
)

option(
  NCRYSTAL_TEST_ENABLE_LONG
  "Enable certain longer tests as part of CTest suite"
  OFF
)

if ( NCRYSTAL_TEST_ENABLE_LONG )
  set( nctest_ignore_pattern "" )
  message( STATUS "Enabling longer tests as well" )
else()
  set( nctest_ignore_pattern "^long_.*$" )
  message(
    STATUS
    "Disabling longer tests. Set NCRYSTAL_TEST_ENABLE_LONG=ON to enable."
  )
endif()

mctools_testutils_add_tests_pyscripts(
  "${CMAKE_CURRENT_LIST_DIR}/scripts" "${nctest_pyenvmod}" "${nctest_ignore_pattern}"
)

set( nctests_extra_link_libs "${nctests_ncrystallib_target}" )

#backtrace from FPE tests need -lexecinfo on FreeBSD:
find_package(Backtrace)
if( Backtrace_FOUND AND Backtrace_LIBRARY )
  #NB: From cmake 3.30, the FindBacktrace modules also provides a proper target
  #(Backtrace::Backtrace) which we should eventually use instead.
  list( APPEND nctests_extra_link_libs "${Backtrace_LIBRARY}" )
endif()

set( nctests_extra_inc_dirs PUBLIC "${CMAKE_CURRENT_LIST_DIR}/include" )
mctools_testutils_add_test_libs(
  "${CMAKE_CURRENT_LIST_DIR}/libs"
  "${nctests_extra_link_libs}"
  "${nctests_extra_inc_dirs}"
)
mctools_testutils_add_test_modules(
  "${CMAKE_CURRENT_LIST_DIR}/modules"
  "${nctests_extra_link_libs}"
  "${nctests_extra_inc_dirs}"
)

set(
  nctest_appenvmod
  "NCRYSTALTEST_DATADIR=set:${NCrystalRoot_SOURCE_DIR}/tests/data"
 )

if ( WIN32 )
  #Do not get "dll not found" errors for the NCrystal library. Note that PATH
  #plays the role of LD_LIBRARY_PATH on Windows (and there is seemingly not
  #really any rpath-like solution).
  list(
    APPEND nctest_appenvmod
    "PATH=path_list_prepend:$<TARGET_FILE_DIR:${nctests_ncrystallib_target}>"
  )
endif()

if ( NOT NCRYSTAL_ENABLE_DATA STREQUAL "EMBED" )
  list(
    APPEND nctest_appenvmod
    "NCRYSTAL_DATADIR=set:${NCrystalRoot_SOURCE_DIR}/data"
  )
endif()

mctools_testutils_add_tests_apps(
  "${CMAKE_CURRENT_LIST_DIR}/src"
  "${nctests_extra_link_libs}"
  "${nctests_extra_inc_dirs}"
  "${nctest_appenvmod}"
)
